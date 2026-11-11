// Experimental BTSplit: bitset filtering, support reachability, and clique coloring.
// Coloring is a standard clique upper bound, not a new BTSplit contribution.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Clock = std::chrono::steady_clock;
struct Expired {};
struct Bits {
    std::vector<uint64_t> w;
    Bits() {}
    explicit Bits(int n) : w((n+63)/64) {}
    void set(int i) { w[i>>6] |= uint64_t(1) << (i&63); }
    void reset(int i) { w[i>>6] &= ~(uint64_t(1) << (i&63)); }
    bool has(int i) const { return (w[i>>6] >> (i&63)) & 1; }
    int count() const { int n=0; for(auto x:w) n+=__builtin_popcountll(x); return n; }
    bool empty() const { for(auto x:w) if(x) return false; return true; }
    int pop() { for(size_t k=0;k<w.size();++k) if(w[k]) {
        int b=__builtin_ctzll(w[k]); w[k]&=w[k]-1; return int(k*64+b);
    } return -1; }
    void intersect(const Bits &b) { for(size_t k=0;k<w.size();++k) w[k]&=b.w[k]; }
    void unite(const Bits &b) { for(size_t k=0;k<w.size();++k) w[k]|=b.w[k]; }
    void subtract(const Bits &b) { for(size_t k=0;k<w.size();++k) w[k]&=~b.w[k]; }
    std::vector<int> vertices() const {
        std::vector<int> v; for(size_t k=0;k<w.size();++k) {
            uint64_t x=w[k]; while(x) { int b=__builtin_ctzll(x); v.push_back(int(k*64+b)); x&=x-1; }
        } return v;
    }
};
static Bits intersection(Bits a,const Bits &b) { a.intersect(b); return a; }
static int intersection_count(const Bits &a,const Bits &b) {
    int n=0; for(size_t k=0;k<a.w.size();++k) n+=__builtin_popcountll(a.w[k]&b.w[k]); return n;
}
struct Graph {
    int n=0; std::vector<Bits> out;
    Graph() {}
    explicit Graph(int size) : n(size),out(size,Bits(size)) {}
    bool edge(int u,int v) const { return out[u].has(v); }
};
struct Options {
    std::string g1,g2,candidates,output,batch,order="degree";
    double timeout=1800, heuristic_seconds=1.0;
    int restarts=256;
    uint64_t seed=1;
    bool reach=true,color=true,components=true;
};
struct Solver {
    Options opt;
    Clock::time_point started;
    Graph a,b;
    std::vector<std::pair<int,int>> pairs;
    std::vector<Bits> compat,link;
    std::vector<std::vector<int>> neighbors;
    std::vector<int> best,chosen;
    uint64_t nodes=0,color_calls=0,color_prunes=0,reach_calls=0,reach_removed=0,rng;
    int component_count=0,largest_component=0,greedy_trials=0;
    bool timed_out=false;
    int stage=0; // input, construction, heuristic, exact search, done
    double input_seconds=0,construction_seconds=0,heuristic_used=0;
    Solver(const Options &o,Clock::time_point t) : opt(o),started(t),rng(o.seed?o.seed:1) {}
    double elapsed() const { return std::chrono::duration<double>(Clock::now()-started).count(); }
    void check() const { if(opt.timeout>0 && elapsed()>=opt.timeout) throw Expired(); }
    uint64_t random() { rng^=rng<<13; rng^=rng>>7; rng^=rng<<17; return rng; }
    Graph read_sip(const std::string &path) {
        std::ifstream in(path); int n;
        if(!(in>>n)||n<0||n>250000) throw std::runtime_error("Invalid SIP dimension");
        if(uint64_t(n)*uint64_t((n+63)/64)*8 > uint64_t(4)*1024*1024*1024)
            throw std::runtime_error("Graph bit matrix exceeds 4 GiB");
        Graph g(n);
        for(int u=0;u<n;++u) {
            int d; if(!(in>>d)||d<0||d>n) throw std::runtime_error("Invalid SIP degree");
            for(int k=0;k<d;++k) { int v; if(!(in>>v)||v<0||v>=n) throw std::runtime_error("Invalid SIP neighbor"); g.out[u].set(v); }
            if((u&255)==0) check();
        }
        return g;
    }
    void load() {
        a=read_sip(opt.g1); b=read_sip(opt.g2);
        std::ifstream in(opt.candidates);
        if(!in) throw std::runtime_error("Cannot read candidates");
        std::string line;
        while(std::getline(in,line)) {
            if(line.find_first_not_of(" \t\r")==std::string::npos) continue;
            std::istringstream row(line); int u,v;
            if(!(row>>u>>v)||u<0||v<0||u>=a.n||v>=b.n) throw std::runtime_error("Invalid candidate");
            if(a.edge(u,u)==b.edge(v,v)) pairs.emplace_back(u,v);
        }
        std::sort(pairs.begin(),pairs.end());
        pairs.erase(std::unique(pairs.begin(),pairs.end()),pairs.end());
        input_seconds=elapsed();
        stage=1;
    }
    bool load_batch(std::istream &in) {
        int n,m; if(!(in>>n>>m)) return false;
        a=Graph(n); b=Graph(m); std::string s;
        for(int u=0;u<n;++u) { in>>s; if(int(s.size())!=n) throw std::runtime_error("Bad matrix"); for(int v=0;v<n;++v) if(s[v]=='1') a.out[u].set(v); }
        for(int u=0;u<m;++u) { in>>s; if(int(s.size())!=m) throw std::runtime_error("Bad matrix"); for(int v=0;v<m;++v) if(s[v]=='1') b.out[u].set(v); }
        for(int u=0;u<n;++u) { in>>s; if(int(s.size())!=m) throw std::runtime_error("Bad candidates"); for(int v=0;v<m;++v) if(s[v]=='1'&&a.edge(u,u)==b.edge(v,v)) pairs.emplace_back(u,v); }
        input_seconds=elapsed(); stage=1; return true;
    }
    void build() {
        const int n=pairs.size();
        if(uint64_t(n)*uint64_t((n+63)/64)*16 > uint64_t(4)*1024*1024*1024)
            throw std::runtime_error("Candidate bit matrices exceed 4 GiB");
        if(n) best={0};
        compat.assign(n,Bits(n)); link.assign(n,Bits(n)); neighbors.resize(n);
        for(int i=0;i<n;++i) {
            if((i&31)==0) check();
            auto [u,v]=pairs[i];
            for(int j=0;j<i;++j) {
                auto [x,y]=pairs[j];
                if(u==x||v==y) continue;
                bool f=a.edge(u,x), r=a.edge(x,u);
                if(f!=b.edge(v,y)||r!=b.edge(y,v)) continue;
                compat[i].set(j); compat[j].set(i);
                if(f||r) {
                    link[i].set(j); link[j].set(i);
                    neighbors[i].push_back(j); neighbors[j].push_back(i);
                }
            }
        }
        construction_seconds=elapsed()-input_seconds;
    }
    std::vector<Bits> components(const Bits &c) {
        check();
        Bits remaining=c; std::vector<std::pair<int,Bits>> collected;
        for(int v;(v=remaining.pop())>=0;) {
            check();
            Bits component(pairs.size()); component.set(v);
            std::vector<int> queue={v};
            uint64_t work=0;
            for(size_t k=0;k<queue.size();++k) for(int u:neighbors[queue[k]]) {
                if((++work&2047)==0) check();
                if(remaining.has(u)) { remaining.reset(u); component.set(u); queue.push_back(u); }
            }
            collected.emplace_back(int(queue.size()),std::move(component));
        }
        std::sort(collected.begin(),collected.end(),[](const auto &x,const auto &y){return x.first>y.first;});
        check();
        std::vector<Bits> result; result.reserve(collected.size());
        for(auto &item:collected) result.push_back(std::move(item.second));
        return result;
    }
    Bits reachable(const Bits &c,const Bits &f) {
        ++reach_calls;
        Bits result=intersection(c,f);
        std::vector<int> queue=result.vertices();
        uint64_t work=0;
        for(size_t k=0;k<queue.size();++k) for(int u:neighbors[queue[k]]) {
            if((++work&2047)==0) check();
            if(c.has(u)&&!result.has(u)) { result.set(u); queue.push_back(u); }
        }
        return result;
    }
    // Returns true only after constructing a complete proper coloring with <= limit colors.
    bool color_prunable(const Bits &c,int limit) {
        ++color_calls;
        check();
        Bits remaining=c;
        uint64_t work=0;
        for(int color=1;!remaining.empty();++color) {
            if((color&31)==0) check();
            if(color>limit) return false;
            Bits available=remaining;
            for(int v;(v=available.pop())>=0;) {
                if((++work&255)==0) check();
                remaining.reset(v);
                available.subtract(compat[v]);
            }
        }
        ++color_prunes; return true;
    }
    int select(const Bits &c,const Bits &f,bool randomized=false) {
        int winner=-1; double best_score=-1e100;
        for(int v:f.vertices()) {
            check();
            int degree=intersection_count(c,link[v]);
            int compatible=intersection_count(c,compat[v]);
            double score=0;
            if(opt.order=="degree") score=degree*1000.0+compatible;
            else if(opt.order=="compat") score=compatible*1000.0+degree;
            else {
                Bits next=intersection(c,compat[v]), nf=f; nf.unite(link[v]); nf.intersect(next);
                int potential=reachable(next,nf).count();
                score=potential*1000.0+compatible+degree*0.1;
            }
            if(randomized) score+=double(random()%1000000)/1000000.0*4000.0;
            if(score>best_score) { best_score=score; winner=v; }
        }
        return winner;
    }
    void improve(const std::vector<int> &solution,const char *phase) {
        if(solution.size()>best.size()) {
            best=solution;
            if(opt.batch.empty()) std::cerr<<"incumbent "<<best.size()<<" elapsed "<<std::setprecision(9)<<elapsed()<<" phase "<<phase<<"\n";
        }
    }
    void greedy(const Bits &component,int seed,bool randomized) {
        std::vector<int> solution={seed};
        Bits c=intersection(component,compat[seed]),f=intersection(c,link[seed]);
        while(!f.empty()) {
            check();
            c=reachable(c,f); f.intersect(c);
            if(solution.size()+c.count()<=best.size()) break;
            int v=select(c,f,randomized);
            if(v<0) break;
            solution.push_back(v); c.intersect(compat[v]); f.unite(link[v]); f.intersect(c);
        }
        improve(solution,"greedy");
    }
    void heuristics(const std::vector<Bits> &groups) {
        double began=elapsed();
        double allowance=std::min(opt.heuristic_seconds,opt.timeout>0?opt.timeout*0.1:opt.heuristic_seconds);
        for(const Bits &component:groups) {
            auto seeds=component.vertices();
            if(seeds.size()<=best.size()) continue;
            std::sort(seeds.begin(),seeds.end(),[&](int x,int y) {
                if(neighbors[x].size()!=neighbors[y].size()) return neighbors[x].size()>neighbors[y].size();
                return x<y;
            });
            int trials=std::min(opt.restarts,int(seeds.size())*2);
            if(pairs.size()<=40) trials=std::min(trials,4);
            for(int k=0;k<trials;++k) {
                if(k>0 && elapsed()-began>=allowance) break;
                int seed=k<int(seeds.size())?seeds[k]:seeds[random()%seeds.size()];
                greedy(component,seed,k>=4); ++greedy_trials;
            }
            if(elapsed()-began>=allowance) break;
        }
        heuristic_used=elapsed()-began;
    }
    void search(Bits c,Bits f) {
        improve(chosen,"search");
        while(!f.empty()) {
            ++nodes; check();
            if(opt.reach) {
                int before=c.count(); c=reachable(c,f); reach_removed+=before-c.count(); f.intersect(c);
            }
            int gap=int(best.size())-int(chosen.size());
            if(c.count()<=gap) return;
            if(opt.color && color_prunable(c,gap)) return;
            int v=select(c,f);
            if(v<0) return;
            Bits nc=intersection(c,compat[v]),nf=f; nf.unite(link[v]); nf.intersect(nc);
            chosen.push_back(v); search(std::move(nc),std::move(nf)); chosen.pop_back();
            c.reset(v); f.reset(v);
        }
    }
    void solve() {
        check();
        build();
        Bits all(pairs.size()); for(size_t i=0;i<pairs.size();++i) all.set(i);
        auto groups=components(all);
        component_count=groups.size();
        largest_component=groups.empty()?0:groups.front().count();
        construction_seconds=elapsed()-input_seconds;
        stage=2;
        heuristics(groups);
        stage=3;
        if(!opt.components) groups={all};
        for(Bits c:groups) {
            while(c.count()>int(best.size())) {
                check();
                // After root exclusions, disconnected components are independent subproblems.
                if(opt.components) {
                    auto split=components(c); Bits useful(pairs.size());
                    for(auto &g:split) if(g.count()>int(best.size())) useful.unite(g);
                    c=std::move(useful);
                    if(c.count()<=int(best.size())) break;
                }
                if(opt.color && color_prunable(c,best.size())) break;
                int v=select(c,c);
                chosen={v};
                Bits nc=intersection(c,compat[v]),nf=intersection(nc,link[v]);
                search(std::move(nc),std::move(nf));
                chosen.clear(); c.reset(v);
            }
        }
        check(); stage=4;
    }
    void emit() {
        if(stage==0) input_seconds=elapsed();
        if(stage==1) construction_seconds=elapsed()-input_seconds;
        if(stage==2) heuristic_used=elapsed()-input_seconds-construction_seconds;
        if(!opt.output.empty()) {
            std::ofstream out(opt.output);
            if(!out) throw std::runtime_error("Cannot write mapping");
            for(int v:best) out<<pairs[v].first<<'\t'<<pairs[v].second<<'\n';
        }
        if(opt.timeout>0 && elapsed()>=opt.timeout) timed_out=true;
        std::cout<<std::setprecision(12)<<"{\"status\":\""<<(timed_out?"timeout":"complete")
                 <<"\",\"solution_size\":"<<best.size()<<",\"nodes\":"<<nodes
                 <<",\"elapsed_seconds\":"<<elapsed()<<",\"input_seconds\":"<<input_seconds
                 <<",\"construction_seconds\":"<<construction_seconds<<",\"heuristic_seconds\":"<<heuristic_used
                 <<",\"greedy_trials\":"<<greedy_trials<<",\"color_calls\":"<<color_calls
                 <<",\"color_prunes\":"<<color_prunes<<",\"reach_calls\":"<<reach_calls
                 <<",\"reach_removed\":"<<reach_removed<<",\"components\":"<<component_count
                 <<",\"largest_component\":"<<largest_component<<",\"candidates\":"<<pairs.size()
                 <<",\"seed\":"<<opt.seed<<",\"pairs\":[";
        bool comma=false; for(int v:best) { if(comma) std::cout<<','; comma=true; std::cout<<'['<<pairs[v].first<<','<<pairs[v].second<<']'; }
        std::cout<<"]}\n"<<std::flush;
    }
};
int main(int argc,char **argv) {
    auto started=Clock::now();
    try {
        Options o;
        for(int i=1;i<argc;++i) {
            std::string key=argv[i];
            if(key=="--no-reach") o.reach=false;
            else if(key=="--no-color") o.color=false;
            else if(key=="--no-components") o.components=false;
            else {
                if(i+1==argc) throw std::runtime_error("Missing argument");
                std::string value=argv[++i];
                if(key=="--g1") o.g1=value;
                else if(key=="--g2") o.g2=value;
                else if(key=="--candidates") o.candidates=value;
                else if(key=="--output") o.output=value;
                else if(key=="--oracle-batch") o.batch=value;
                else if(key=="--timeout") o.timeout=std::stod(value);
                else if(key=="--seed") o.seed=std::stoull(value);
                else if(key=="--restarts") o.restarts=std::stoi(value);
                else if(key=="--heuristic-seconds") o.heuristic_seconds=std::stod(value);
                else if(key=="--order") o.order=value;
                else throw std::runtime_error("Unknown argument: "+key);
            }
        }
        if(!o.batch.empty()) {
            std::ifstream in(o.batch); if(!in) throw std::runtime_error("Cannot read batch");
            for(;;) {
                Solver solver(o,Clock::now());
                if(!solver.load_batch(in)) break;
                try { solver.solve(); } catch(Expired &) {solver.timed_out=true;}
                solver.emit();
            }
        } else {
            Solver solver(o,started);
            try {solver.load(); solver.solve();} catch(Expired &) {solver.timed_out=true;}
            solver.emit();
        }
        return 0;
    } catch(const std::exception &e) {
        std::cerr<<e.what()<<'\n'; std::cout<<"{\"status\":\"error\",\"solution_size\":0,\"nodes\":0,\"pairs\":[]}\n"; return 2;
    }
}
