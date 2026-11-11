// Persistent coverage cache for current BTSplit. K/H, Reach, coloring and
// weighted-degree search share one implementation in fresh and cached modes.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <psapi.h>
#endif

using Clock = std::chrono::steady_clock;
struct Stopped { std::string status; };
struct Bits {
    std::vector<uint64_t> w;
    Bits() = default;
    explicit Bits(int n):w((n+63)/64) {}
    void set(int i) { w[i>>6] |= uint64_t(1)<<(i&63); }
    void reset(int i) { w[i>>6] &= ~(uint64_t(1)<<(i&63)); }
    bool has(int i) const { return (w[i>>6]>>(i&63))&1; }
    bool empty() const { for(auto x:w) if(x) return false; return true; }
    int count() const { int n=0; for(auto x:w) n+=__builtin_popcountll(x); return n; }
    int pop() { for(size_t k=0;k<w.size();++k) if(w[k]) {
        int b=__builtin_ctzll(w[k]); w[k]&=w[k]-1; return int(k*64+b);
    } return -1; }
    void intersect(const Bits &b) { for(size_t k=0;k<w.size();++k) w[k]&=b.w[k]; }
    void unite(const Bits &b) { for(size_t k=0;k<w.size();++k) w[k]|=b.w[k]; }
    void subtract(const Bits &b) { for(size_t k=0;k<w.size();++k) w[k]&=~b.w[k]; }
    bool operator==(const Bits &b) const { return w==b.w; }
    bool numeric_less(const Bits &b) const {
        for(size_t k=w.size();k>0;--k) if(w[k-1]!=b.w[k-1]) return w[k-1]<b.w[k-1];
        return false;
    }
    std::vector<int> vertices() const {
        std::vector<int> result;
        for(size_t k=0;k<w.size();++k) { auto x=w[k]; while(x) {
            int b=__builtin_ctzll(x); result.push_back(int(64*k+b)); x&=x-1;
        }} return result;
    }
};
struct BitsHash { size_t operator()(const Bits &b) const {
    size_t h=1469598103934665603ULL;
    for(auto x:b.w) { h^=size_t(x); h*=1099511628211ULL; }
    return h;
}};
static Bits intersect(Bits a,const Bits &b) { a.intersect(b); return a; }
static int count_intersection(const Bits &a,const Bits &b) {
    int n=0; for(size_t k=0;k<a.w.size();++k) n+=__builtin_popcountll(a.w[k]&b.w[k]); return n;
}
static uint64_t rss_bytes() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX m{};
    if(GetProcessMemoryInfo(GetCurrentProcess(),reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&m),sizeof(m)))
        return uint64_t(m.WorkingSetSize);
#endif
    return 0;
}
struct Deadline {
    Clock::time_point start=Clock::now(),last_memory=start;
    double limit=1800;
    uint64_t memory_limit=uint64_t(8)*1024*1024*1024,peak_rss=0;
#ifdef BTSPLIT_AUDIT
    int fail_after_checks=-1; // Absent from production builds; deterministic fault injection.
#endif
    void begin(double seconds) { start=Clock::now(); limit=seconds; last_memory=start; }
    double elapsed() const { return std::chrono::duration<double>(Clock::now()-start).count(); }
    void check_memory() {
        auto value=rss_bytes(); peak_rss=std::max(peak_rss,value);
        if(memory_limit && value>memory_limit) throw Stopped{"resource_limit"};
    }
    void check() {
#ifdef BTSPLIT_AUDIT
        if(fail_after_checks>=0 && fail_after_checks--==0) throw Stopped{"timeout"};
#endif
        auto now=Clock::now();
        if(limit>=0 && std::chrono::duration<double>(now-start).count()>=limit) throw Stopped{"timeout"};
        if(std::chrono::duration<double>(now-last_memory).count()>=0.1) { check_memory(); last_memory=now; }
    }
};
struct Options {
    std::string g1,g2,universe,initial,updates,kh_case,mode="cached";
    double budget=1800,heuristic_seconds=1;
    int restarts=256;
    uint64_t seed=1,memory_limit=uint64_t(8)*1024*1024*1024;
    bool dump_cache=false;
};
struct Edit { bool insertion; int id; bool fresh; };
struct Graph {
    int n=0; std::vector<Bits> out;
    explicit Graph(int size=0):n(size),out(size,Bits(size)) {}
    bool edge(int u,int v) const { return out[u].has(v); }
};
struct Problem {
    int n=0;
    std::vector<std::pair<int,int>> pairs;
    std::vector<Bits> k,h;
    std::vector<std::vector<int>> neighbors;
    Bits initial,eligible;
    std::vector<Edit> edits;
    Deadline &deadline;
    explicit Problem(Deadline &d):deadline(d) {}
    Graph read_sip(const std::string &path) {
        std::ifstream in(path); int size;
        if(!(in>>size)||size<0||size>250000) throw std::runtime_error("Invalid SIP dimension: "+path);
        if(uint64_t(size)*((size+63)/64)*8>uint64_t(4)*1024*1024*1024)
            throw std::runtime_error("Input graph matrix exceeds 4 GiB");
        Graph g(size);
        for(int u=0;u<size;++u) {
            deadline.check(); int degree;
            if(!(in>>degree)||degree<0||degree>size) throw std::runtime_error("Invalid SIP degree");
            for(int j=0;j<degree;++j) { int v;
                if(!(in>>v)||v<0||v>=size) throw std::runtime_error("Invalid SIP endpoint");
                g.out[u].set(v);
            }
        }
        std::string extra; if(in>>extra) throw std::runtime_error("Trailing SIP data");
        return g;
    }
    std::vector<std::pair<int,int>> read_pairs(const std::string &path,int na,int nb) {
        std::ifstream in(path); if(!in) throw std::runtime_error("Cannot read pairs: "+path);
        std::vector<std::pair<int,int>> result; std::string line;
        while(std::getline(in,line)) {
            deadline.check(); if(line.find_first_not_of(" \t\r")==std::string::npos) continue;
            std::istringstream row(line); int u,v; std::string extra;
            if(!(row>>u>>v)||u<0||v<0||u>=na||v>=nb||(row>>extra)) throw std::runtime_error("Invalid pair row");
            result.emplace_back(u,v);
        }
        std::sort(result.begin(),result.end());
        result.erase(std::unique(result.begin(),result.end()),result.end()); return result;
    }
    void allocate() {
        n=int(pairs.size());
        if(uint64_t(n)*((n+63)/64)*16>uint64_t(4)*1024*1024*1024)
            throw std::runtime_error("Candidate matrices exceed 4 GiB");
        k.assign(n,Bits(n)); h.assign(n,Bits(n)); neighbors.resize(n);
        initial=Bits(n); eligible=Bits(n);
    }
    void build_neighbors() {
        for(int v=0;v<n;++v) { deadline.check(); neighbors[v]=h[v].vertices(); }
    }
    void load(const Options &o) {
        Graph a=read_sip(o.g1),b=read_sip(o.g2);
        pairs=read_pairs(o.universe,a.n,b.n); allocate();
        std::map<std::pair<int,int>,int> ids;
        for(int i=0;i<n;++i) { ids[pairs[i]]=i;
            if(a.edge(pairs[i].first,pairs[i].first)==b.edge(pairs[i].second,pairs[i].second)) eligible.set(i);
        }
        for(auto p:read_pairs(o.initial,a.n,b.n)) {
            auto it=ids.find(p); if(it==ids.end()) throw std::runtime_error("Initial pair outside universe"); initial.set(it->second);
        }
        for(int i=0;i<n;++i) {
            deadline.check(); if(!eligible.has(i)) continue;
            auto [u,v]=pairs[i];
            for(int j=0;j<i;++j) { if((j&1023)==0) deadline.check(); if(!eligible.has(j)) continue;
                auto [x,y]=pairs[j]; if(u==x||v==y) continue;
                bool f=a.edge(u,x),r=a.edge(x,u);
                if(f!=b.edge(v,y)||r!=b.edge(y,v)) continue;
                k[i].set(j); k[j].set(i);
                if(f||r) { h[i].set(j); h[j].set(i); }
            }
        }
        build_neighbors();
        std::ifstream in(o.updates); if(!in) throw std::runtime_error("Cannot read updates");
        std::string line;
        while(std::getline(in,line)) {
            if(line.find_first_not_of(" \t\r")==std::string::npos) continue;
            std::istringstream row(line); std::string op,extra; int u,v,fresh;
            if(!(row>>op>>u>>v>>fresh)||(row>>extra)||(fresh!=0&&fresh!=1)||(op!="insert"&&op!="delete"))
                throw std::runtime_error("Invalid update row");
            auto it=ids.find({u,v}); if(it==ids.end()) throw std::runtime_error("Edited pair outside universe");
            edits.push_back({op=="insert",it->second,bool(fresh)});
        }
        validate_edits();
    }
    void load_kh(const std::string &path) {
        std::ifstream in(path); int size;
        if(!(in>>size)||size<1||size>10000) throw std::runtime_error("Invalid KH dimension");
        for(int i=0;i<size;++i) pairs.emplace_back(i,i);
        allocate();
        for(auto *matrix:{&k,&h}) for(int i=0;i<n;++i) {
            std::string row; in>>row; if(int(row.size())!=n) throw std::runtime_error("Invalid KH row");
            for(int j=0;j<n;++j) { if(row[j]!='0'&&row[j]!='1') throw std::runtime_error("Invalid KH bit");
                if(row[j]=='1') (*matrix)[i].set(j);
            }
        }
        for(int i=0;i<n;++i) { eligible.set(i);
            if(k[i].has(i)||h[i].has(i)) throw std::runtime_error("KH has self edge");
            for(int j=0;j<n;++j) if(k[i].has(j)!=k[j].has(i)||h[i].has(j)!=h[j].has(i)||(h[i].has(j)&&!k[i].has(j)))
                throw std::runtime_error("KH must be symmetric with H subset K");
        }
        std::string row; in>>row; if(int(row.size())!=n) throw std::runtime_error("Invalid initial KH row");
        for(int i=0;i<n;++i) { if(row[i]!='0'&&row[i]!='1') throw std::runtime_error("Invalid initial bit"); if(row[i]=='1') initial.set(i); }
        int count; if(!(in>>count)||count<0) throw std::runtime_error("Invalid edit count");
        for(int i=0;i<count;++i) { std::string op; int id,fresh;
            if(!(in>>op>>id>>fresh)||id<0||id>=n||(op!="insert"&&op!="delete")||(fresh!=0&&fresh!=1))
                throw std::runtime_error("Invalid KH edit");
            edits.push_back({op=="insert",id,bool(fresh)});
        }
        build_neighbors(); validate_edits();
    }
    void validate_edits() {
        Bits active=initial;
        for(auto e:edits) { deadline.check(); if(active.has(e.id)==e.insertion) throw std::runtime_error("Edit does not toggle presence");
            if(e.insertion) active.set(e.id); else active.reset(e.id);
        }
    }
};
struct Region { Bits m,c; int upper; const char *kind; };
struct Stats {
    uint64_t nodes=0,visits=0,color_calls=0,color_prunes=0,cardinality_prunes=0;
    uint64_t reach_calls=0,reach_removed=0,root_skips=0,resumed=0,invalid_regions=0,restricted_regions=0,removed_solutions=0;
    uint64_t greedy_trials=0; bool early_return=false;
};
struct Engine {
    const Problem &p; const Options &opt; Deadline &deadline; bool capture;
    Bits active,best;
    std::unordered_set<Bits,BitsHash> pool;
    std::vector<Region> q;
    Stats stats; uint64_t rng;
    Engine(const Problem &problem,const Options &o,Deadline &d,bool save):
        p(problem),opt(o),deadline(d),capture(save),active(p.n),best(p.n),rng(o.seed?o.seed:1) {}
    void record(const Bits &m) {
        deadline.check(); ++stats.visits;
        if(capture) pool.insert(m);
        if(m.count()>best.count()) best=m;
    }
    void retain(const Bits &m,const Bits &c,int upper,const char *kind) {
        if(capture) q.push_back({m,c,upper,kind});
        deadline.check();
    }
    Bits boundary(const Bits &m,const Bits &c) {
        if(m.empty()) return c;
        Bits f(p.n); for(int v:m.vertices()) { deadline.check(); f.unite(p.h[v]); }
        f.intersect(c); return f;
    }
    std::vector<Bits> components(const Bits &c) {
        deadline.check(); Bits remaining=c; std::vector<Bits> result;
        for(int v;(v=remaining.pop())>=0;) {
            Bits group(p.n); group.set(v); std::vector<int> todo={v};
            for(size_t i=0;i<todo.size();++i) { deadline.check();
                for(int u:p.neighbors[todo[i]]) if(remaining.has(u)) { remaining.reset(u); group.set(u); todo.push_back(u); }
            }
            result.push_back(std::move(group));
        }
        std::stable_sort(result.begin(),result.end(),[](const Bits &a,const Bits &b){return a.count()>b.count();});
        return result;
    }
    Bits reachable(const Bits &c,const Bits &f) {
        ++stats.reach_calls; Bits result=intersect(c,f); auto todo=result.vertices();
        for(size_t i=0;i<todo.size();++i) { deadline.check();
            for(int u:p.neighbors[todo[i]]) if(c.has(u)&&!result.has(u)) {result.set(u); todo.push_back(u);}
        }
        return result;
    }
    // -1 means coloring was stopped early and does NOT provide a bound.
    int coloring_if_prunable(const Bits &c,int limit) {
        ++stats.color_calls; Bits remaining=c; int colors=0;
        while(!remaining.empty()) {
            deadline.check(); if(++colors>limit) return -1;
            Bits available=remaining;
            for(int v;(v=available.pop())>=0;) { deadline.check(); remaining.reset(v); available.subtract(p.k[v]); }
        }
        ++stats.color_prunes; return colors;
    }
    uint64_t random() { rng^=rng<<13; rng^=rng>>7; rng^=rng<<17; return rng; }
    int select(const Bits &c,const Bits &f,bool randomized=false) {
        int winner=-1; double winning=-1;
        for(int v:f.vertices()) {
            deadline.check();
            double score=1000.0*count_intersection(c,p.h[v])+count_intersection(c,p.k[v]);
            if(randomized) score+=double(random()%1000000)/1000000.0*4000.0;
            if(score>winning) { winning=score; winner=v; }
        }
        return winner;
    }
    void greedy(const Bits &group,int seed,bool randomized) {
        Bits m(p.n); m.set(seed); record(m);
        Bits c=intersect(group,p.k[seed]),f=intersect(c,p.h[seed]);
        while(!f.empty()) {
            deadline.check(); c=reachable(c,f); f.intersect(c);
            if(m.count()+c.count()<=best.count()) break;
            int v=select(c,f,randomized); if(v<0) break;
            m.set(v); record(m); c.intersect(p.k[v]); f.unite(p.h[v]); f.intersect(c);
        }
    }
    void heuristics(const std::vector<Bits> &groups) {
        if(opt.restarts<=0||opt.heuristic_seconds<=0) return;
        double began=deadline.elapsed();
        double allowance=std::min(opt.heuristic_seconds,deadline.limit>=0?deadline.limit*0.1:opt.heuristic_seconds);
        for(const auto &group:groups) {
            auto seeds=group.vertices(); if(int(seeds.size())<=best.count()) continue;
            std::sort(seeds.begin(),seeds.end(),[&](int x,int y) {
                int a=count_intersection(active,p.h[x]),b=count_intersection(active,p.h[y]);
                return a!=b?a>b:x<y;
            });
            int trials=std::min(opt.restarts,int(seeds.size())*2);
            if(active.count()<=40) trials=std::min(trials,4);
            for(int i=0;i<trials;++i) {
                if(i>0&&deadline.elapsed()-began>=allowance) break;
                int seed=i<int(seeds.size())?seeds[i]:seeds[random()%seeds.size()];
                greedy(group,seed,i>=4); ++stats.greedy_trials;
            }
            if(deadline.elapsed()-began>=allowance) break;
        }
    }
    void search(const Bits &m,Bits c) {
        record(m); Bits f=boundary(m,c); int selected=m.count();
        while(!f.empty()) {
            deadline.check(); ++stats.nodes;
            if(selected) {
                int before=c.count(); c=reachable(c,f); stats.reach_removed+=before-c.count(); f.intersect(c);
            } else {
                // Residual root components may split after an exclusion. Every
                // skipped component remains in Q, including never-visited roots.
                Bits useful(p.n);
                for(const auto &group:components(c)) {
                    if(group.count()<=best.count()) { ++stats.root_skips; retain(m,group,group.count(),"root_size"); }
                    else useful.unite(group);
                }
                c=std::move(useful); f=c;
            }
            if(c.empty()||f.empty()) return;
            int gap=best.count()-selected;
            if(c.count()<=gap) { ++stats.cardinality_prunes; retain(m,c,selected+c.count(),"cardinality"); return; }
            int colors=coloring_if_prunable(c,gap);
            if(colors>=0) { retain(m,c,selected+colors,"color"); return; }
            int v=select(c,f); if(v<0) return;
            Bits next_m=m; next_m.set(v);
            search(next_m,intersect(c,p.k[v]));
            c.reset(v); f.reset(v);
        }
    }
    void initialize(const Bits &raw_active) {
        active=intersect(raw_active,p.eligible); best=Bits(p.n); pool.clear(); q.clear(); stats=Stats{};
        rng=opt.seed?opt.seed:1; record(Bits(p.n));
        auto groups=components(active); heuristics(groups);
        for(const auto &group:groups) {
            deadline.check();
            if(group.count()<=best.count()) { ++stats.root_skips; retain(Bits(p.n),group,group.count(),"root_size"); }
            else search(Bits(p.n),group);
        }
        deadline.check();
    }
    void insert(int v) {
        stats=Stats{};
        if(!p.eligible.has(v)) { deadline.check(); return; }
        if(active.has(v)) throw std::runtime_error("Inserting active candidate");
        active.set(v); deadline.check(); deadline.check_memory();
        Bits m(p.n); m.set(v);
        search(m,intersect(active,p.k[v])); deadline.check();
    }
    void erase(int v) {
        stats=Stats{};
        if(!p.eligible.has(v)) { deadline.check(); return; }
        if(!active.has(v)) throw std::runtime_error("Deleting inactive candidate");
        Bits old_best=best; active.reset(v);
        // A timeout during maintenance must still expose a feasible lower bound.
        if(old_best.has(v)) best=Bits(p.n);
        deadline.check(); deadline.check_memory();
        for(auto it=pool.begin();it!=pool.end();) {
            deadline.check();
            if(it->has(v)) { it=pool.erase(it); ++stats.removed_solutions; } else ++it;
        }
        std::vector<Region> surviving; surviving.reserve(q.size());
        for(auto &region:q) {
            deadline.check();
            if(region.m.has(v)) { ++stats.invalid_regions; continue; }
            if(region.c.has(v)) { region.c.reset(v); ++stats.restricted_regions; }
            surviving.push_back(std::move(region)); // U unchanged and still safe.
        }
        q=std::move(surviving);
        if(!old_best.has(v)) { stats.early_return=true; deadline.check(); return; }
        best=Bits(p.n);
        for(const auto &m:pool) { deadline.check();
            if(m.count()>best.count()||(m.count()==best.count()&&m.numeric_less(best))) best=m;
        }
        // New frontier produced by search has U <= the then-current incumbent.
        // It need not be reopened again while this incumbent only increases.
        auto old_q=std::move(q); q.clear();
        for(auto &region:old_q) {
            deadline.check();
            if(region.upper>best.count()) { ++stats.resumed; search(region.m,std::move(region.c)); }
            else q.push_back(std::move(region));
        }
        deadline.check();
    }
};
static void emit_ids(const Bits &b) {
    std::cout<<'['; bool comma=false;
    for(int v:b.vertices()) { if(comma) std::cout<<','; comma=true; std::cout<<v; } std::cout<<']';
}
static void emit_record(const Problem &p,const Options &o,const Engine &e,int index,const Edit *edit,
                        const std::string &status,double seconds,double cumulative,int raw_active,uint64_t rss) {
    const auto &s=e.stats;
    std::cout<<std::setprecision(12)<<"{\"type\":\"result\",\"mode\":\""<<o.mode<<"\",\"update_index\":"<<index
        <<",\"operation\":\""<<(edit?(edit->insertion?"insert":"delete"):"initial")<<"\",\"status\":\""<<status
        <<"\",\"seconds\":"<<seconds<<",\"cached_cumulative_seconds\":"<<cumulative
        <<",\"active_count\":"<<raw_active<<",\"effective_active_count\":"<<e.active.count()
        <<",\"solution_size\":"<<e.best.count()<<",\"nodes\":"<<s.nodes<<",\"visits\":"<<s.visits
        <<",\"color_calls\":"<<s.color_calls<<",\"color_prunes\":"<<s.color_prunes
        <<",\"cardinality_prunes\":"<<s.cardinality_prunes<<",\"reach_calls\":"<<s.reach_calls
        <<",\"reach_removed\":"<<s.reach_removed<<",\"root_skips\":"<<s.root_skips
        <<",\"resumed_regions\":"<<s.resumed<<",\"invalid_regions\":"<<s.invalid_regions
        <<",\"restricted_regions\":"<<s.restricted_regions<<",\"removed_solutions\":"<<s.removed_solutions
        <<",\"greedy_trials\":"<<s.greedy_trials<<",\"early_return\":"<<(s.early_return?"true":"false")
        <<",\"pool_size\":"<<e.pool.size()<<",\"frontier_size\":"<<e.q.size()<<",\"peak_rss_bytes\":"<<rss;
    if(edit) std::cout<<",\"edited_id\":"<<edit->id<<",\"fresh_checkpoint\":"<<(edit->fresh?"true":"false");
    std::cout<<",\"solution_ids\":"; emit_ids(e.best);
    std::cout<<",\"pairs\":["; bool comma=false;
    for(int v:e.best.vertices()) { if(comma) std::cout<<','; comma=true; std::cout<<'['<<p.pairs[v].first<<','<<p.pairs[v].second<<']'; }
    std::cout<<']';
    if(o.dump_cache&&e.capture&&status=="complete") {
        std::cout<<",\"pool\":["; std::vector<Bits> ordered(e.pool.begin(),e.pool.end());
        std::sort(ordered.begin(),ordered.end(),[](const Bits &a,const Bits &b){return a.numeric_less(b);});
        comma=false; for(const auto &m:ordered) { if(comma) std::cout<<','; comma=true; emit_ids(m); }
        std::cout<<"],\"frontier\":["; comma=false;
        for(const auto &r:e.q) { if(comma) std::cout<<','; comma=true;
            std::cout<<"{\"m\":"; emit_ids(r.m); std::cout<<",\"c\":"; emit_ids(r.c);
            std::cout<<",\"upper\":"<<r.upper<<",\"kind\":\""<<r.kind<<"\"}";
        }
        std::cout<<']';
    }
    std::cout<<"}\n"<<std::flush;
}
int main(int argc,char **argv) {
    auto process_start=Clock::now();
    try {
        Options o;
        for(int i=1;i<argc;++i) {
            std::string key=argv[i];
            if(key=="--dump-cache") { o.dump_cache=true; continue; }
            if(i+1==argc) throw std::runtime_error("Missing option value");
            std::string value=argv[++i];
            if(key=="--g1") o.g1=value; else if(key=="--g2") o.g2=value;
            else if(key=="--universe") o.universe=value; else if(key=="--initial") o.initial=value;
            else if(key=="--updates") o.updates=value; else if(key=="--mode") o.mode=value;
            else if(key=="--kh-case") o.kh_case=value; else if(key=="--budget") o.budget=std::stod(value);
            else if(key=="--heuristic-seconds") o.heuristic_seconds=std::stod(value);
            else if(key=="--restarts") o.restarts=std::stoi(value); else if(key=="--seed") o.seed=std::stoull(value);
            else if(key=="--rss-limit-bytes") o.memory_limit=std::stoull(value);
            else throw std::runtime_error("Unknown option: "+key);
        }
        if(o.mode!="cached"&&o.mode!="fresh") throw std::runtime_error("Mode must be cached or fresh");
        if(o.budget<=0||o.restarts<0||o.heuristic_seconds<0) throw std::runtime_error("Invalid numeric option");
        Deadline deadline; deadline.memory_limit=o.memory_limit; deadline.begin(o.budget);
        Problem p(deadline);
        if(o.kh_case.empty()) p.load(o); else p.load_kh(o.kh_case);
        deadline.check(); double preprocess=deadline.elapsed();
        std::cout<<std::setprecision(12)<<"{\"type\":\"metadata\",\"mode\":\""<<o.mode
            <<"\",\"universe_size\":"<<p.n<<",\"eligible_candidates\":"<<p.eligible.count()
            <<",\"edits\":"<<p.edits.size()<<",\"preprocessing_seconds\":"<<preprocess
            <<",\"budget_seconds\":"<<o.budget<<",\"seed\":"<<o.seed<<",\"restarts\":"<<o.restarts
            <<",\"heuristic_seconds\":"<<o.heuristic_seconds<<",\"score\":\"1000H+K\",\"cache_schema\":1}\n"<<std::flush;
        Bits raw_active=p.initial; double cumulative=preprocess; int processed=0,certified=0; bool failed=false;
        if(o.mode=="cached") {
            deadline.begin(std::max(0.0,o.budget-cumulative));
            Engine engine(p,o,deadline,true);
            for(size_t index=0;index<=p.edits.size();++index) {
                std::string status="complete"; const Edit *edit=index?&p.edits[index-1]:nullptr;
                if(index) deadline.begin(std::max(0.0,o.budget-cumulative));
                try {
                    if(!edit) { deadline.check_memory(); engine.initialize(raw_active); }
                    else if(edit->insertion) { raw_active.set(edit->id); engine.insert(edit->id); }
                    else { raw_active.reset(edit->id); engine.erase(edit->id); }
                    deadline.check(); deadline.check_memory();
                } catch(const Stopped &stop) { status=stop.status; failed=true; }
                double seconds=deadline.elapsed(); cumulative+=seconds;
                emit_record(p,o,engine,int(index),edit,status,seconds,cumulative,raw_active.count(),deadline.peak_rss);
                if(index) { ++processed; if(status=="complete") ++certified; }
                if(failed) break;
            }
        } else {
            for(size_t i=0;i<p.edits.size();++i) {
                const auto &edit=p.edits[i]; if(edit.insertion) raw_active.set(edit.id); else raw_active.reset(edit.id);
                if(!edit.fresh) continue;
                deadline.begin(o.budget); Engine engine(p,o,deadline,false); std::string status="complete";
                try { deadline.check_memory(); engine.initialize(raw_active); deadline.check(); deadline.check_memory(); }
                catch(const Stopped &stop) { status=stop.status; failed=true; }
                emit_record(p,o,engine,int(i+1),&edit,status,deadline.elapsed(),0,raw_active.count(),deadline.peak_rss);
                ++processed; if(status=="complete") ++certified;
                // A failed fresh point does not invalidate subsequent independent points.
            }
        }
        std::cout<<std::setprecision(12)<<"{\"type\":\"summary\",\"mode\":\""<<o.mode<<"\",\"status\":\""
            <<(failed?"incomplete":"complete")<<"\",\"processed_update_records\":"<<processed
            <<",\"certified_update_records\":"<<certified
            <<",\"process_seconds\":"<<std::chrono::duration<double>(Clock::now()-process_start).count()
            <<",\"peak_rss_bytes\":"<<deadline.peak_rss<<"}\n"<<std::flush;
        return failed?3:0;
    } catch(const Stopped &stop) {
        std::cout<<"{\"type\":\"fatal\",\"status\":\""<<stop.status<<"\",\"phase\":\"preprocessing\"}\n"<<std::flush; return 3;
    } catch(const std::exception &error) {
        std::cerr<<error.what()<<'\n';
        std::cout<<"{\"type\":\"fatal\",\"status\":\"error\"}\n"<<std::flush; return 2;
    }
}
