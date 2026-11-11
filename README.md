# BTSplit

An exact solver for Constraint-Aware Maximum Common Induced Subgraph (CA-MCIS)
on directed graphs, with support for candidate insertions and deletions.

- `src/btsplit.cpp`: static solver.
- `src/btsplit_dynamic.cpp`: cached and fresh solving under candidate updates.

## Build

Requires a GCC-compatible C++17 compiler and GNU Make. On Linux/macOS, use a
GCC-compatible toolchain and run `make`. On Windows with MinGW, run
`mingw32-make`. Windows executables have the `.exe` suffix.

## Static solving

```sh
./btsplit --g1 graph1.sip --g2 graph2.sip --candidates pairs.tsv --timeout 1800 --output solution.tsv
```

Graphs use directed SIP format: the first line is the vertex count; each
following line contains the out-degree and zero-based outgoing neighbor IDs
of one vertex. Candidate files contain two zero-based vertex IDs per line,
without a header. The solver writes a JSON result to standard output and an
optional matching to `--output`. `complete` certifies an optimum; unfinished
results provide a feasible lower bound.

## Candidate updates

```sh
./btsplit_dynamic --g1 graph1.sip --g2 graph2.sip --universe universe.tsv --initial initial.tsv --updates updates.tsv --mode cached --budget 1800
```

`universe.tsv` contains every candidate that may appear; `initial.tsv` specifies
the initial active subset. Each update row is `delete|insert u v fresh`, where
`fresh` is a 0/1 checkpoint flag. Use `--mode fresh` for independent solves at
flagged checkpoints. Results are written as JSON lines to standard output.

The dynamic solver's internal memory monitor uses Windows PSAPI; on other
systems, apply an external memory limit when needed.
