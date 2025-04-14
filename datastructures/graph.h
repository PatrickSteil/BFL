/*
 * Licensed under MIT License.
 * Author: Patrick Steil
 */

#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "types.h"

struct Graph {
  std::vector<std::vector<Vertex>> edges;

  Graph(const std::size_t numVertices = 0) : edges(numVertices) {}

  bool isVertex(const Vertex v) const { return v < edges.size(); }

  std::size_t numVertices() const { return edges.size(); }

  std::size_t numEdges() const {
    return std::accumulate(edges.begin(), edges.end(), std::size_t(0),
                           [](std::size_t sum, const std::vector<Vertex>& v) {
                             return sum + v.size();
                           });
  }

  std::size_t degree(const Vertex v) const {
    assert(isVertex(v));
    return edges[v].size();
  }

  void addEdge(const Vertex from, const Vertex to) {
    assert(isVertex(from));
    assert(isVertex(to));
    edges[from].push_back(to);
  }

  template <typename FUNC>
  void relaxEdges(const Vertex v, const FUNC&& function) const {
    for (auto w : edges[v]) {
      function(v, w);
    }
  }

  void showStats() const {
    std::size_t n = numVertices();
    std::size_t m = numEdges();

    std::size_t minDegree = (n > 0 ? edges[0].size() : 0);
    std::size_t maxDegree = 0;
    std::size_t sumDegree = 0;
    std::size_t isolatedCount = 0;

    for (const auto& neighbors : edges) {
      std::size_t deg = neighbors.size();
      minDegree = std::min(minDegree, deg);
      maxDegree = std::max(maxDegree, deg);
      sumDegree += deg;
      if (deg == 0) {
        ++isolatedCount;
      }
    }

    double averageDegree = (n > 0 ? static_cast<double>(sumDegree) / n : 0.0);

    std::cout << "Graph Statistics:" << std::endl;
    std::cout << "  Number of vertices:            " << n << std::endl;
    std::cout << "  Number of edges:               " << m << std::endl;
    std::cout << "  Min degree:                    " << minDegree << std::endl;
    std::cout << "  Max degree:                    " << maxDegree << std::endl;
    std::cout << "  Average degree:                " << std::fixed
              << std::setprecision(5) << averageDegree << std::endl;
    std::cout << "  Number of isolated vertices:   " << isolatedCount
              << std::endl;
  }

  Graph reversed() const {
    Graph rev(numVertices());
    for (std::size_t u = 0; u < edges.size(); ++u) {
      for (std::uint32_t v : edges[u]) {
        rev.edges[v].push_back(static_cast<Vertex>(u));
      }
    }
    return rev;
  }

  static Graph readFromDimacs(std::istream& in) {
    std::string line;
    Graph g;

    while (std::getline(in, line)) {
      if (line.empty()) continue;
      if (line[0] == 'c') continue;

      std::istringstream iss(line);
      std::string token;
      iss >> token;

      if (token == "p") {
        std::string graphType;
        std::size_t numVertices = 0, numEdges = 0;
        if (!(iss >> graphType >> numVertices >> numEdges)) {
          throw std::runtime_error("Invalid problem line in DIMACS file.");
        }
        g = Graph(numVertices);
      } else if (token == "a" || token == "e") {
        Vertex u, v;
        if (!(iss >> u >> v)) {
          throw std::runtime_error("Invalid edge line in DIMACS file.");
        }

        if (u == 0 || v == 0) {
          throw std::runtime_error(
              "DIMACS file vertices must be 1-indexed (nonzero).");
        }
        g.addEdge(u - 1, v - 1);
      }
    }
    return g;
  }

  static Graph readFromDimacsFile(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open file: " + filename);
    return readFromDimacs(in);
  }

  static Graph readFromEdgeList(std::istream& in) {
    std::string line;
    std::vector<std::pair<Vertex, Vertex>> edges;

    edges.reserve(100000);

    Vertex maxV = 0;

    while (std::getline(in, line)) {
      if (line.empty()) continue;
      if (line[0] == '#') continue;

      std::istringstream iss(line);
      Vertex from, to;
      iss >> from >> to;
      edges.emplace_back(from, to);

      maxV = std::max({from, to, maxV});
    }

    Graph g(maxV + 1);

    for (const auto& e : edges) {
      g.addEdge(e.first, e.second);
    }
    return g;
  }

  static Graph readFromEdgeListFile(const std::string& filename) {
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open file: " + filename);
    return readFromEdgeList(in);
  }
};
