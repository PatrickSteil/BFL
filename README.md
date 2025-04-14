
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
Reading graph from dimacs ... done [19242ms]
Graph Statistics:
  Number of vertices:            5861608
  Number of edges:               12515291
  Min degree:                    0
  Max degree:                    31
  Average degree:                2.13513
  Number of isolated vertices:   1618
DFS Visit Times ... done [690ms]
Merge Vertices ... done [34ms]
Build Index ... done [1734ms]
Memory Consumption:
  Labels memory: 178.88208 mb
  Discovery/Finish times memory: 44.72052 mb
  Total memory: 223.60260 mb
Simple BFS Benchmark Results:
  Number of queries:            1000
  Total time (ms):              127856
  Time per query (ms):          127.85600
  Number of positive queries:   361
Label-based BFS Benchmark Results:
  Number of queries:            1000
  Total time (ms):              36064
  Time per query (ms):          36.06400
  Number of positive queries:   361
Label-based DFS Benchmark Results:
  Number of queries:            1000
  Total time (ms):              56853
  Time per query (ms):          56.85300
  Number of positive queries:   361
```

## Reference
The BFL algorithm is based on:
- [**"Reachability Querying: Can It Be Even Faster?"** – Su et al.](https://ieeexplore.ieee.org/document/7750623)