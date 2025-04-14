/*
 * Licensed under MIT License.
 * Author: Patrick Steil
 */

#pragma once

#include <bitset>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#include "generation_checker.h"
#include "graph.h"
#include "status_log.h"
#include "types.h"

template <int S_MAX = 64>
struct BFL {
  std::array<const Graph*, 2> graphs;
  std::array<std::vector<std::bitset<S_MAX>>, 2> labels;

  std::vector<std::uint32_t> discoveredTime;
  std::vector<std::uint32_t> finishedTime;

  std::vector<Vertex> postOrderTraversalOrder;

  GenerationChecker<> visited;
  std::vector<Vertex> stack;
  std::size_t index;

  std::vector<std::uint32_t> hash;
  std::vector<Vertex> g;

  BFL(const Graph& fwdGraph, const Graph& bwdGraph)
      : graphs{&fwdGraph, &bwdGraph},
        labels{std::vector<std::bitset<S_MAX>>(fwdGraph.numVertices(), 0),
               std::vector<std::bitset<S_MAX>>(fwdGraph.numVertices(), 0)},
        discoveredTime(fwdGraph.numVertices(), 0),
        finishedTime(fwdGraph.numVertices(), 0),
        postOrderTraversalOrder(),
        visited(fwdGraph.numVertices()),
        stack(fwdGraph.numVertices()),
        index(0),
        hash(fwdGraph.numVertices()),
        g(fwdGraph.numVertices()) {
    for (Vertex v = 0; v < fwdGraph.numVertices(); ++v) {
      hash[v] = hashVertex(v);
    }
  }

  // Wrapper to reset index before DFS usage
  void resetStack() {
    index = 0;
    visited.reset();
  }

  void buildIndex() {
    StatusLog log("Build Index");
    visited.reset();

    std::function<void(std::uint32_t, DIRECTION)> compute =
        [&](const std::uint32_t v, const DIRECTION dir = FWD) -> void {
      visited.mark(v);
      labels[dir][v][hash[v]] = 1;

      graphs[FWD]->relaxEdges(v, [&](const Vertex, const Vertex w) {
        if (!visited.isMarked(w)) compute(w, dir);

        labels[dir][v] |= labels[dir][w];
      });
    };

    for (std::uint32_t v = 0; v < graphs[FWD]->numVertices(); ++v) {
      if (visited.isMarked(v)) continue;
      compute(v, FWD);
    }

    visited.reset();
    for (std::uint32_t v = 0; v < graphs[FWD]->numVertices(); ++v) {
      if (visited.isMarked(v)) continue;
      compute(v, BWD);
    }
  }
  std::uint16_t hashVertex(const Vertex v) const {
    return std::hash<std::uint32_t>{}(v) % S_MAX;
  }

  void mergeVertices(const int d = 1600) {
    StatusLog log("Merge Vertices");

    std::size_t n = graphs[FWD]->numVertices();
    std::size_t sliceSize = n / d;
    std::size_t index = 0;
    for (std::uint16_t slice = 0; slice < d; ++slice) {
      if (index >= n) break;

      std::uint32_t rep = postOrderTraversalOrder[index];

      std::size_t currentSliceSize = sliceSize;
      if (slice == d - 1) {
        currentSliceSize = n - index;
      }

      for (std::size_t j = 0; j < currentSliceSize; ++j, ++index) {
        g[postOrderTraversalOrder[index]] = rep;
        hash[postOrderTraversalOrder[index]] = hash[rep];
      }
    }
  }

  void dfsVisit(const Vertex v, std::uint32_t& timeCounter) {
    visited.mark(v);
    discoveredTime[v] = timeCounter++;

    graphs[FWD]->relaxEdges(v, [&](const Vertex, const Vertex w) {
      if (!visited.isMarked(w)) {
        dfsVisit(w, timeCounter);
      }
    });
    postOrderTraversalOrder.push_back(v);
    finishedTime[v] = timeCounter++;
  }

  void computeTimesAndOrder() {
    StatusLog log("DFS Visit Times");
    postOrderTraversalOrder.reserve(graphs[FWD]->numVertices());

    std::uint32_t timeCounter = 0;

    for (Vertex v = 0; v < graphs[FWD]->numVertices(); v++) {
      if (graphs[BWD]->degree(v) == 0) {
        dfsVisit(v, timeCounter);
      }
    }
  }

  bool queryIter(const Vertex from, const Vertex to) {
    assert(graphs[FWD]->isVertex(from));
    assert(graphs[FWD]->isVertex(to));

    if (from == to) [[unlikely]]
      return true;

    visited.reset();
    resetStack();

    stack[index++] = from;
    visited.mark(from);

    const auto& targetFWDLabels = labels[FWD][to];
    const auto& targetBWDLabels = labels[BWD][to];
    const auto targetDiscovered = discoveredTime[to];
    const auto targetFinished = finishedTime[to];

    while (index > 0) {
      assert(index < stack.size());
      Vertex curr = stack[--index];

      if (discoveredTime[curr] <= targetDiscovered &&
          targetFinished <= finishedTime[curr]) {
        return true;
      }

      if ((targetFWDLabels & ~labels[FWD][curr]).any() ||
          (labels[BWD][curr] & ~targetBWDLabels).any()) {
        continue;
      }

      for (const Vertex w : graphs[FWD]->edges[curr]) {
        if (!visited.isMarked(w)) {
          visited.mark(w);
          stack[index++] = w;
        }
      }
    }
    return false;
  }

  bool query(const Vertex from, const Vertex to) {
    assert(graphs[FWD]->isVertex(from));
    assert(graphs[FWD]->isVertex(to));

    if (from == to) [[unlikely]]
      return true;
    visited.reset();

    visited.mark(from);

    return queryRec(from, to);
  }

  bool queryRec(const Vertex from, const Vertex to) {
    if (discoveredTime[from] <= discoveredTime[to] &&
        finishedTime[to] <= finishedTime[from]) {
      return true;
    }

    if (((labels[FWD][to] & labels[FWD][from]) != labels[FWD][to]) ||
        ((labels[BWD][from] & labels[BWD][to]) != labels[BWD][from])) {
      return false;
    }

    for (const Vertex w : graphs[FWD]->edges[from]) {
      if (!visited.isMarked(w)) {
        visited.mark(w);
        if (queryRec(w, to)) {
          return true;
        }
      }
    }

    return false;
  }

  void run_benchmark(const std::size_t numberOfQueries = 10000) {
    std::vector<std::pair<Vertex, Vertex>> queries;
    queries.reserve(numberOfQueries);

    std::size_t n = graphs[FWD]->numVertices();
    if (n == 0) {
      std::cout << "Graph is empty.\n";
      return;
    }

    std::mt19937 gen(42);
    std::uniform_int_distribution<Vertex> dis(0, static_cast<Vertex>(n - 1));

    for (std::size_t i = 0; i < numberOfQueries; i++) {
      Vertex u = dis(gen);
      Vertex v = dis(gen);
      queries.emplace_back(u, v);
    }

    // Run the queries and measure the overall time in milliseconds.
    auto startTime = std::chrono::high_resolution_clock::now();

    std::size_t positiveCount = 0;
    for (const auto& pr : queries) {
      positiveCount += query(pr.first, pr.second);
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    auto elapsedNanos = std::chrono::duration_cast<std::chrono::milliseconds>(
                            endTime - startTime)
                            .count();

    double timePerQuery = static_cast<double>(elapsedNanos) / numberOfQueries;

    std::cout << "Benchmark Results:\n";
    std::cout << "  Number of queries:            " << numberOfQueries << "\n";
    std::cout << "  Total time (ms):              " << elapsedNanos << "\n";
    std::cout << "  Time per query (ms):          " << timePerQuery << "\n";
    std::cout << "  Number of positive queries:   " << positiveCount << "\n";
  }
};