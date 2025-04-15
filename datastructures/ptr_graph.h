#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "types.h"
// A simple Graph class that stores, for each vertex,
// pointers to contiguous arrays for outgoing and incoming vertex ids.
struct Graph {
  std::size_t numVertices;
  // For each direction (FWD and BWD), we store a vector of pointers
  // to contiguous arrays of neighbor ids.
  std::array<std::vector<int*>, 2> edges;
  // For each direction, we store the degree (i.e. number of neighbors) for each
  // vertex.
  std::array<std::vector<std::size_t>, 2> degrees;

  std::vector<std::uint16_t> lastSeen;
  std::uint16_t timestamp;

  Graph(std::size_t n)
      : numVertices(n),
        edges({std::vector<int*>(n, nullptr), std::vector<int*>(n, nullptr)}),
        degrees(
            {std::vector<std::size_t>(n, 0), std::vector<std::size_t>(n, 0)}),
        lastSeen(n, 0),
        timestamp(0) {}

  ~Graph() {
    for (const auto& vec : edges) {
      for (auto ptr : vec) {
        delete[] ptr;
      }
    }
  }

  void reset() {
    ++timestamp;

    if (timestamp == 0) {
      std::fill(lastSeen.begin(), lastSeen.end(), 0);
      ++timestamp;
    }
  }

  // Check whether v is a valid vertex.
  bool isVertex(int v) const {
    return (v >= 0) && (static_cast<std::size_t>(v) < numVertices);
  }

  // Added member function showStats() prints out:
  //   Graph Statistics:
  //     Number of vertices:            [numVertices]
  //     Number of edges:               [totalEdges]
  //     Min degree:                    [minDegree]
  //     Max degree:                    [maxDegree]
  //     Average degree:                [avgDegree with 5 decimal precision]
  //     Number of isolated vertices:   [numIsolated]
  void showStats() const {
    std::size_t totalEdges = 0;
    std::size_t minDegree = numVertices > 0 ? static_cast<std::size_t>(-1) : 0;
    std::size_t maxDegree = 0;
    std::size_t numIsolated = 0;

    // We'll use the forward degrees for our computations.
    const auto& deg = degrees[FWD];
    for (std::size_t i = 0; i < numVertices; ++i) {
      std::size_t d = deg[i];
      totalEdges += d;
      if (d < minDegree) {
        minDegree = d;
      }
      if (d > maxDegree) {
        maxDegree = d;
      }
      if (d == 0) {
        numIsolated++;
      }
    }
    double avgDegree =
        numVertices ? static_cast<double>(totalEdges) / numVertices : 0.0;

    std::cout << "Graph Statistics:\n";
    std::cout << "  Number of vertices:            " << numVertices << "\n";
    std::cout << "  Number of edges:               " << totalEdges << "\n";
    std::cout << "  Min degree:                    " << minDegree << "\n";
    std::cout << "  Max degree:                    " << maxDegree << "\n";
    std::cout << "  Average degree:                " << std::fixed
              << std::setprecision(5) << avgDegree << "\n";
    std::cout << "  Number of isolated vertices:   " << numIsolated << "\n";
  }
};

//
// Reads a DIMACS graph file.
// The expected format is:
//   c <any comment text>
//   p edge <num_vertices> <num_edges>
//   e <u> <v>
// where vertices u and v are 1-indexed and each edge is defined by a line
// starting with 'e'. This function builds both the outgoing neighbor lists (for
// FWD) and the incoming neighbor lists (for BWD).
//

Graph readDimacsGraphFromFile(const std::string filename) {
  std::ifstream inFile(filename);
  std::string line;
  int n = 0, m = 0;

  // Read header lines: skip comment lines and process the problem line.
  while (std::getline(inFile, line)) {
    if (line.empty()) continue;
    // In DIMACS sp files, comments might start with 'c'
    if (line[0] == 'c') continue;
    if (line[0] == 'p') {
      std::istringstream iss(line);
      std::string pStr, typeStr;
      iss >> pStr >> typeStr >> n >> m;
      break;
    }
  }

  Graph g(n);

  // Temporary storage for accumulating outgoing and incoming edges.
  std::vector<std::vector<int>> tempOut(n);
  std::vector<std::vector<int>> tempIn(n);

  // Process arc lines.
  while (std::getline(inFile, line)) {
    if (line.empty() || line[0] == 'c') continue;
    // Expecting lines starting with 'a'
    if (line[0] == 'a') {
      std::istringstream iss(line);
      char a;
      int u, v;
      if (!(iss >> a >> u >> v)) {
        std::cerr << "Error reading arc from line: " << line << "\n";
        continue;
      }
      // Convert from 1-indexed to 0-indexed.
      u--;
      v--;
      if (!g.isVertex(u) || !g.isVertex(v)) continue;
      tempOut[u].push_back(v);
      tempIn[v].push_back(u);
    }
  }
  inFile.close();

  // Allocate contiguous arrays for each vertex and copy the temporary storage.
  for (std::size_t i = 0; i < g.numVertices; ++i) {
    // Outgoing neighbors.
    g.degrees[FWD][i] = tempOut[i].size();
    if (g.degrees[FWD][i] > 0) {
      g.edges[FWD][i] = new int[g.degrees[FWD][i]];
      std::copy(tempOut[i].begin(), tempOut[i].end(), g.edges[FWD][i]);
    }
    // Incoming neighbors.
    g.degrees[BWD][i] = tempIn[i].size();
    if (g.degrees[BWD][i] > 0) {
      g.edges[BWD][i] = new int[g.degrees[BWD][i]];
      std::copy(tempIn[i].begin(), tempIn[i].end(), g.edges[BWD][i]);
    }
  }

  return g;
}