/*
 * Licensed under MIT License.
 * Author: Patrick Steil
 */

#pragma once

#ifdef __GNUC__
#define PREFETCH(addr) __builtin_prefetch(addr)
#else
#define PREFETCH(addr)
#endif

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "status_log.h"
#include "types.h"
#include "utils.h"

struct DynGraph {
  std::vector<std::vector<Vertex>> edges;

  DynGraph(const std::size_t numVertices = 0) : edges(numVertices) {}

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
      if (function(v, w)) break;
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

    std::cout << "DynGraph Statistics:" << std::endl;
    std::cout << "  Number of vertices:            " << n << std::endl;
    std::cout << "  Number of edges:               " << m << std::endl;
    std::cout << "  Min degree:                    " << minDegree << std::endl;
    std::cout << "  Max degree:                    " << maxDegree << std::endl;
    std::cout << "  Average degree:                " << std::fixed
              << std::setprecision(5) << averageDegree << std::endl;
    std::cout << "  Number of isolated vertices:   " << isolatedCount
              << std::endl;
  }

  DynGraph reversed() const {
    DynGraph rev(numVertices());
    for (std::size_t u = 0; u < edges.size(); ++u) {
      for (std::uint32_t v : edges[u]) {
        rev.edges[v].push_back(static_cast<Vertex>(u));
      }
    }
    return rev;
  }

  std::vector<Vertex> topologicalOrder() const {
    std::vector<bool> visited(numVertices(), false);
    std::vector<Vertex> order;
    order.reserve(numVertices());

    std::function<void(Vertex)> dfs = [&](Vertex v) {
      visited[v] = true;
      for (const Vertex nbr : edges[v]) {
        if (!visited[nbr]) {
          dfs(nbr);
        }
      }
      order.push_back(v);
    };

    for (Vertex v = 0; v < numVertices(); ++v) {
      if (!visited[v]) {
        dfs(v);
      }
    }

    std::reverse(order.begin(), order.end());
    return order;
  }

  static DynGraph readFromDimacs(std::istream& in) {
    std::string line;
    DynGraph g;

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
        g = DynGraph(numVertices);
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

  static DynGraph readFromDimacsFile(const std::string& filename) {
    StatusLog log("Reading graph from dimacs");
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open file: " + filename);
    return readFromDimacs(in);
  }

  static DynGraph readFromEdgeList(std::istream& in) {
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

    DynGraph g(maxV + 1);

    for (const auto& e : edges) {
      g.addEdge(e.first, e.second);
    }
    return g;
  }

  static DynGraph readFromEdgeListFile(const std::string& filename) {
    StatusLog log("Reading graph from edge list");
    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open file: " + filename);
    return readFromEdgeList(in);
  }

  void printAsGReach() const {
    std::cout << "graph_for_greach\n";
    std::cout << numVertices() << "\n";

    for (Vertex v = 0; v < numVertices(); ++v) {
      std::cout << v << ": ";
      for (Vertex w : edges[v]) std::cout << w << " ";
      std::cout << "#\n";
    }
  }
};

std::vector<std::pair<Vertex, Vertex>> getTopologicallySortedEdges(
    const DynGraph& g) {
  StatusLog log("Compute Topo-Ordering and sort edges");
  std::vector<std::pair<Vertex, Vertex>> sortedEdges;
  sortedEdges.reserve(g.numEdges());

  std::vector<Vertex> topoOrder = g.topologicalOrder();

  std::vector<std::size_t> position(g.numVertices());
  for (std::size_t i = 0; i < topoOrder.size(); ++i) {
    position[topoOrder[i]] = i;
  }

  for (Vertex u = 0; u < g.numVertices(); ++u) {
    for (const Vertex v : g.edges[u]) {
      sortedEdges.emplace_back(u, v);
    }
  }
  std::sort(sortedEdges.begin(), sortedEdges.end(),
            [&position](const std::pair<Vertex, Vertex>& left,
                        const std::pair<Vertex, Vertex>& right) {
              return std::tie(position[left.first], position[left.second]) <
                     std::tie(position[right.first], position[right.second]);
            });

  return sortedEdges;
}

struct Edge {
  Vertex from;
  Vertex to;

  Edge() = default;
  Edge(Vertex from, Vertex to) : from(from), to(to) {}

  auto operator<=>(const Edge& other) const = default;
};

struct Graph {
  std::vector<std::size_t> adjArray;
  std::vector<Vertex> toVertex;

  Graph() : adjArray(1), toVertex() {};
  Graph(std::size_t numVertices, std::size_t numEdges)
      : adjArray(numVertices + 1, 0), toVertex(numEdges, 0) {};

  Graph(const Graph& other)
      : adjArray(other.adjArray), toVertex(other.toVertex) {};

  Graph(Graph&& other) noexcept
      : adjArray(std::move(other.adjArray)),
        toVertex(std::move(other.toVertex)) {}

  Graph& operator=(Graph&& other) noexcept {
    if (this != &other) {
      adjArray = std::move(other.adjArray);
      toVertex = std::move(other.toVertex);
    }
    return *this;
  }

  bool isVertex(const Vertex v) const { return v < numVertices(); }

  std::size_t numVertices() const { return adjArray.size() - 1; }
  std::size_t numEdges() const { return toVertex.size(); }

  void print() const {
    std::cout << "NumVertices: " << numVertices() << std::endl;
    std::cout << "NumEdges: " << numEdges() << std::endl;

    for (Vertex v = 0; v < numVertices(); ++v) {
      std::cout << "Edges from " << v << std::endl;

      for (std::size_t i = beginEdge(v); i < endEdge(v); ++i) {
        std::cout << toVertex[i] << " ";
      }
      std::cout << std::endl;
    }
  }

  template <typename FUNC>
  void doForAllEdges(FUNC&& function) const {
    for (Vertex v = 0; v < numVertices(); ++v) {
      for (std::size_t i = beginEdge(v); i < endEdge(v); ++i) {
        if (i + 4 < endEdge(v)) {
          PREFETCH(&toVertex[i + 4]);
        }

        function(v, toVertex[i]);
      }
    }
  }

  template <typename FUNC>
  void relaxEdges(const Vertex from, FUNC&& function) const {
    for (std::size_t i = beginEdge(from); i < endEdge(from); ++i) {
      if (i + 4 < endEdge(from)) {
        PREFETCH(&toVertex[i + 4]);
      }

      if (function(from, toVertex[i])) break;
    }
  }

  std::size_t degree(const Vertex v) const {
    assert(isVertex(v));
    return endEdge(v) - beginEdge(v);
  }

  std::size_t beginEdge(const Vertex v) const {
    assert(isVertex(v));
    assert(v < adjArray.size());
    return adjArray[v];
  }

  std::size_t endEdge(const Vertex v) const {
    assert(isVertex(v));
    assert(v + 1 < adjArray.size());
    return adjArray[v + 1];
  }

  void clear() {
    adjArray.clear();
    toVertex.clear();
  }

  void sortByRank(const std::vector<std::size_t>& rank) {
    for (Vertex v = 0; v < adjArray.size() - 1; ++v) {
      std::sort(toVertex.begin() + adjArray[v],
                toVertex.begin() + adjArray[v + 1],
                [&](const Vertex left, const Vertex right) {
                  return rank[left] < rank[right];
                });
    }
  }

  void buildFromEdgeList(std::vector<std::pair<Vertex, Vertex>>& edgeList,
                         std::size_t numVertices) {
    std::vector<std::size_t> newAdjArray(numVertices + 1, 0);
    std::vector<Vertex> newToVertex;

    sortAndRemoveDuplicates(edgeList);

    assert(std::is_sorted(edgeList.begin(), edgeList.end(),
                          [](const auto& a, const auto& b) {
                            return (a.first < b.first) ||
                                   (a.first == b.first && a.second < b.second);
                          }));

    for (const auto& [u, v] : edgeList) {
      assert(u < newAdjArray.size());
      assert(v < newAdjArray.size());
      assert(u < numVertices);
      assert(v < numVertices);

      ++newAdjArray[u + 1];
    }

    for (std::size_t i = 1; i < newAdjArray.size(); ++i) {
      newAdjArray[i] += newAdjArray[i - 1];
    }

    assert(newAdjArray.back() == edgeList.size());

    newToVertex.resize(edgeList.size());
    std::vector<std::size_t> offset = newAdjArray;
    for (const auto& [u, v] : edgeList) {
      assert(u < offset.size());
      assert(offset[u] < newAdjArray[u + 1]);
      assert(offset[u] < newToVertex.size());
      newToVertex[offset[u]++] = v;
    }

    adjArray = std::move(newAdjArray);
    toVertex = std::move(newToVertex);
  }

  void readFromEdgeList(const std::string& fileName) {
    clear();

    std::ifstream file(fileName);
    if (!file.is_open()) {
      throw std::runtime_error("Cannot open file: " + fileName);
    }

    std::vector<std::pair<Vertex, Vertex>> edges;
    Vertex u, v, maxV;
    std::string line;
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      if (!(iss >> u >> v)) {
        continue;
      }
      edges.emplace_back(u - 1, v - 1);

      maxV = std::max({maxV, u, v});
    }
    file.close();

    buildFromEdgeList(edges, maxV);
  }

  void readDimacs(const std::string& fileName) {
    StatusLog log("Reading graph from dimacs");
    clear();

    std::ifstream file(fileName);
    if (!file.is_open()) {
      throw std::runtime_error("Cannot open file: " + fileName);
    }

    std::string line;
    std::vector<std::pair<Vertex, Vertex>> edges;
    Vertex numVertices = 0, numEdges = 0;

    while (std::getline(file, line)) {
      if (line.empty() || line[0] == 'c') {
        continue;
      }

      if (line[0] == 'p') {
        std::istringstream iss(line);
        std::string tmp;
        if (iss >> tmp >> tmp >> numVertices >> numEdges) {
          adjArray.assign(numVertices + 1, std::size_t(0));
          toVertex.assign(numEdges, Vertex(0));
          edges.reserve(numEdges);
        }
      } else if (line[0] == 'a') {
        std::istringstream iss(line);
        char a;
        Vertex u, v;
        if (iss >> a >> u >> v) {
          edges.emplace_back(u - 1, v - 1);
        }
      }
    }

    file.close();
    std::sort(edges.begin(), edges.end(),
              [](const auto& left, const auto& right) {
                return std::tie(left.first, left.second) <
                       std::tie(right.first, right.second);
              });

    for (const auto& [u, v] : edges) {
      ++adjArray[u + 1];
    }

    for (std::size_t i = 1; i < adjArray.size(); ++i) {
      adjArray[i] += adjArray[i - 1];
    }

    adjArray.back() = edges.size();

    toVertex.resize(edges.size());
    std::vector<std::size_t> offset = adjArray;

    for (const auto& [u, v] : edges) {
      toVertex[offset[u]++] = v;
    }
  }

  void toDimacs(const std::string& fileName) const {
    std::ofstream file(fileName);
    if (!file.is_open()) {
      throw std::runtime_error("Cannot open file: " + fileName);
    }

    file << "p sp " << numVertices() << " " << numEdges() << "\n";

    for (Vertex v = 0; v < numVertices(); ++v) {
      for (std::size_t i = beginEdge(v); i < endEdge(v); ++i) {
        file << "a " << (v + 1) << " " << (toVertex[i] + 1) << "\n";
      }
    }

    file.close();
  }

  void readSnap(const std::string& fileName) {
    StatusLog log("Reading graph from .snap format");
    clear();

    std::ifstream file(fileName);
    if (!file.is_open()) {
      throw std::runtime_error("Cannot open file: " + fileName);
    }

    std::vector<std::pair<Vertex, Vertex>> edges;
    Vertex maxVertex = 0;

    std::string line;
    while (std::getline(file, line)) {
      if (line.empty() || line[0] == '#') {
        continue;
      }

      std::istringstream iss(line);
      Vertex u, v;
      if (!(iss >> u >> v)) {
        throw std::runtime_error("Invalid line format in .snap file: " + line);
      }

      edges.emplace_back(u, v);
      maxVertex = std::max(maxVertex, v);
    }

    file.close();

    adjArray.resize(maxVertex + 2, 0);

    std::sort(edges.begin(), edges.end(),
              [](const auto& left, const auto& right) {
                return std::tie(left.first, left.second) <
                       std::tie(right.first, right.second);
              });

    auto it = std::unique(edges.begin(), edges.end());
    edges.erase(it, edges.end());

    for (const auto& [u, v] : edges) {
      ++adjArray[u + 1];
    }

    for (std::size_t i = 1; i < adjArray.size(); ++i) {
      adjArray[i] += adjArray[i - 1];
    }

    adjArray.back() = edges.size();

    toVertex.resize(edges.size());
    std::vector<std::size_t> offset = adjArray;

    for (const auto& [u, v] : edges) {
      toVertex[offset[u]++] = v;
    }
  }

  bool rankIsPermutation(const std::vector<std::size_t>& rank) {
    const std::size_t n = rank.size();
    std::vector<bool> seen(n, false);
    for (std::size_t i = 0; i < n; ++i) {
      const std::size_t r = rank[i];

      if (!(r < n && !seen[r])) return false;
      seen[r] = true;
    }

    return std::all_of(seen.begin(), seen.end(), [](bool val) { return val; });
  }

  void reorderByRank(const std::vector<std::size_t>& rank) {
    assert(rankIsPermutation(rank));
    assert(rank.size() == numVertices());

    std::vector<std::size_t> newAdjArray(numVertices() + 1, 0);
    std::vector<Vertex> newToVertex(numEdges(), 0);

    for (std::size_t v = 0; v < numVertices(); ++v) {
      newAdjArray[rank[v] + 1] += degree(v);
    }

    for (std::size_t v = 0; v < numVertices(); ++v) {
      newAdjArray[v + 1] += newAdjArray[v];
    }

    std::vector<std::size_t> placedEdges(numVertices(), 0);

    for (std::size_t v = 0; v < numVertices(); ++v) {
      relaxEdges(v, [&](const Vertex /* from */, const Vertex to) {
        newToVertex[newAdjArray[rank[v]] + placedEdges[rank[v]]] = rank[to];
        placedEdges[rank[v]]++;

        return false;
      });
    }

    for (std::size_t v = 0; v < numVertices(); ++v) {
      std::sort(newToVertex.begin() + newAdjArray[v],
                newToVertex.begin() + newAdjArray[v + 1]);
    }

    adjArray = std::move(newAdjArray);
    toVertex = std::move(newToVertex);
  }

  Graph reverseGraph() const {
    Graph reversed;
    reversed.adjArray = adjArray;
    reversed.toVertex = toVertex;
    reversed.flip();
    return reversed;
  }

  void flip() {
    std::vector<std::size_t> flippedAdjArray(numVertices() + 1, 0);
    std::vector<Vertex> flippedToVertex(numEdges(), noVertex);

    for (Vertex fromV(0); fromV < numVertices(); ++fromV) {
      for (std::size_t i = adjArray[fromV]; i < adjArray[fromV + 1]; ++i) {
        flippedAdjArray[toVertex[i] + 1]++;
      }
    }

    for (Vertex v = 1; v <= numVertices(); ++v) {
      flippedAdjArray[v] += flippedAdjArray[v - 1];
    }

    std::vector<std::size_t> offset = flippedAdjArray;

    for (Vertex fromV(0); fromV < numVertices(); ++fromV) {
      for (std::size_t i = adjArray[fromV]; i < adjArray[fromV + 1]; ++i) {
        Vertex toV = toVertex[i];
        flippedToVertex[offset[toV]++] = fromV;
      }
    }

    adjArray = std::move(flippedAdjArray);
    toVertex = std::move(flippedToVertex);
  }

  void showStats() const {
    if (numEdges() == 0) {
      std::cout << "Graph is empty.\n";
      return;
    }

    std::size_t numIsolatedVertices = 0;
    std::size_t minDegree = std::numeric_limits<std::size_t>::max();
    std::size_t maxDegree = 0;
    std::size_t totalDegree = 0;

    for (Vertex v = 0; v < numVertices(); ++v) {
      std::size_t deg = degree(v);
      numIsolatedVertices += (deg == 0);

      minDegree = std::min(minDegree, deg);
      maxDegree = std::max(maxDegree, deg);
      totalDegree += deg;
    }

    double avgDegree = static_cast<double>(totalDegree) / numVertices();

    std::cout << "Graph Statistics:\n";
    std::cout << "  Number of vertices:            " << numVertices() << "\n";
    std::cout << "  Number of edges:               " << numEdges() << "\n";
    std::cout << "  Min degree:                    " << minDegree << "\n";
    std::cout << "  Max degree:                    " << maxDegree << "\n";
    std::cout << "  Average degree:                " << avgDegree << "\n";
    std::cout << "  Number of isolated vertices:   " << numIsolatedVertices
              << "\n";
  }

  template <typename FUNC>
  void removeEdges(const FUNC&& predicate) {
    std::vector<Vertex> newToVertex;
    std::vector<std::size_t> newAdjArray(numVertices() + 1, 0);

    for (Vertex v = 0; v < numVertices(); ++v) {
      std::size_t begin = beginEdge(v);
      std::size_t end = endEdge(v);

      for (std::size_t i = begin; i < end; ++i) {
        if (i + 4 < end) {
          PREFETCH(&toVertex[i + 4]);
        }

        if (!predicate(v, toVertex[i])) {
          newToVertex.push_back(toVertex[i]);
        }
      }
      newAdjArray[v + 1] = newToVertex.size();
    }

    toVertex = std::move(newToVertex);
    adjArray = std::move(newAdjArray);
  }

  std::vector<Vertex> removeVertices(
      const std::vector<std::uint8_t>& partition,
      const std::vector<Vertex>& representation) {
    assert(partition.size() == numVertices());
    assert(representation.size() == numVertices());

    auto keepVertex = [&](const Vertex v) -> bool {
      return (partition[v] == 3 || representation[v] == v);
    };

    const std::size_t oldNumVertices = numVertices();

    std::vector<Vertex> oldToNew(oldNumVertices, static_cast<Vertex>(-1));
    std::size_t newNumVertices = 0;
    for (Vertex u = 0; u < oldNumVertices; ++u) {
      if (keepVertex(u)) oldToNew[u] = newNumVertices++;
    }

    std::vector<std::size_t> newAdjArray(newNumVertices + 1, 0);
    std::size_t newNumEdges = 0;

    for (Vertex u = 0; u < oldNumVertices; ++u) {
      if (!keepVertex(u)) continue;
      newAdjArray[oldToNew[u]] = newNumEdges;
      for (std::size_t j = beginEdge(u); j < endEdge(u); ++j) {
        Vertex w = toVertex[j];
        if (keepVertex(w)) ++newNumEdges;
      }
    }
    newAdjArray[newNumVertices] = newNumEdges;

    std::vector<Vertex> newToVertex(newNumEdges, 0);
    std::size_t edgeIndex = 0;
    for (Vertex u = 0; u < oldNumVertices; ++u) {
      if (!keepVertex(u)) continue;
      for (std::size_t j = beginEdge(u); j < endEdge(u); ++j) {
        Vertex w = toVertex[j];
        if (keepVertex(w)) {
          newToVertex[edgeIndex++] = oldToNew[w];
        }
      }
    }
    assert(edgeIndex == newNumEdges);

    adjArray = std::move(newAdjArray);
    toVertex = std::move(newToVertex);

    return oldToNew;
  }
};