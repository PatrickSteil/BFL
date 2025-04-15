
# Bloom Filter Labeling (BFL)
Author: Patrick Steil

This repository implements a reachability-labeling algorithm for efficient reachability queries on directed acyclic graphs (DAGs). 

## Build Instructions
To build the project, run the provided compile.sh script:
```bash
./compile.sh
```

Alternatively, to build only the release version:
```bash
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

The project supports two build configurations:
- Release: Optimized for performance.
- Debug: Includes debug symbols and error checks.

To build and run the release versions:
```bash
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
./BFL -h
```

## Example Execution
Here is an example of running the BFL algorithm on a graph read from a DIMACS file:
```bash
./BFL -i ../example/swiss.dimacs -s -b
```

Example Output:
```bash
Reading graph from dimacs ... done [12069ms]
Graph Statistics:
  Number of vertices:            5861608
  Number of edges:               12515291
  Min degree:                    0
  Max degree:                    31
  Average degree:                2.13513
  Number of isolated vertices:   1618
Build Index ... done [1143ms]
Memory Consumption:
  Labels memory: 447.205 mb
  Discovery/Finish times memory: 89.441 mb
  Total memory: 536.646 mb
Label-based DFS-Rec Benchmark Results:
  Number of queries:            10000
  Total time (ms):              5651
  Time per query (ms):          0.5651
  Number of positive queries:   3520
Label-based DFS-Iter Benchmark Results:
  Number of queries:            10000
  Total time (ms):              5110
  Time per query (ms):          0.511
  Number of positive queries:   3520
```

## Reference
The BFL algorithm is based on:
- [**"Reachability Querying: Can It Be Even Faster?"** – Su et al.](https://ieeexplore.ieee.org/document/7750623)