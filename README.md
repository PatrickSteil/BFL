
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
./BFL -i ../example/kvv.dimacs -s -b
```

Example Output:
```bash
Reading graph from dimacs ... done [588ms]
Graph Statistics:
  Number of vertices:  300223
  Number of edges: 708101
  Min degree:  0
  Max degree:  13
  Average degree:  2.35858
  Number of isolated vertices: 97
DFS Visit Times ... done [25ms]
Merge Vertices ... done [1ms]
Build Index ... done [65ms]
Benchmark Results:
  Number of queries:  10000
  Total time (ms):  42127
  Time per query (ms):  4.21270
  Number of positive queries: 4487
```

## Reference
The BFL algorithm is based on:
- [**"Reachability Querying: Can It Be Even Faster?"** – Su et al.](https://ieeexplore.ieee.org/document/7750623)