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

#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "generation_checker.h"
#include "graph.h"
#include "status_log.h"
#include "types.h"

template <int S_MAX = 64>
struct BFL {
  static constexpr int NUM_WORDS = (S_MAX) / 64;
  struct Label {
    std::array<std::array<std::uint64_t, NUM_WORDS>, 2> bits;
    std::pair<std::uint32_t, std::uint32_t> times = {0, 0};
    std::pair<std::uint16_t, std::uint32_t> degrees = {0, 0};

    Label() {
      bits[FWD].fill(0);
      bits[BWD].fill(0);
    }

    void set(const DIRECTION dir, std::size_t pos) {
      std::size_t word = pos / 64;
      std::size_t offset = pos % 64;
      bits[FWD][word] |= (1ULL << offset);
    }
  };

  std::array<const Graph*, 2> graphs;
  std::vector<Label> labels;

  GenerationChecker<> visited;

  std::vector<Vertex> stack;
  std::size_t timestamp;
  std::size_t index;

  BFL(const Graph& fwdGraph, const Graph& bwdGraph)
      : graphs{&fwdGraph, &bwdGraph},
        labels(fwdGraph.numVertices(), Label()),
        visited(fwdGraph.numVertices()),
        stack(fwdGraph.numVertices(), 0),
        timestamp(0),
        index(0) {
    for (Vertex v = 0; v < fwdGraph.numVertices(); ++v) {
      labels[v].degrees.first = fwdGraph.degree(v);
      labels[v].degrees.second = bwdGraph.degree(v);
    }
  }

  void printMemoryConsumption() const {
    std::size_t numVertices = graphs[FWD]->numVertices();

    std::size_t labelsMemory = 2 * numVertices * sizeof(Label);
    // std::size_t timesMemory = numVertices * sizeof(std::pair<std::uint32_t,
    // std::uint32_t>); std::size_t totalMemory = labelsMemory + timesMemory;

    std::cout << "Memory Consumption:\n";
    // std::cout << "  Labels memory: "
    //           << static_cast<double>(labelsMemory) / (1024.0 * 1024.0)
    //           << " mb\n";
    // std::cout << "  Discovery/Finish times memory: "
    //           << static_cast<double>(timesMemory) / (1024.0 * 1024.0)
    //           << " mb\n";
    std::cout << "  Total memory: "
              << static_cast<double>(labelsMemory) / (1024.0 * 1024.0)
              << " mb\n";
  }

  std::size_t nextHash() {
    static std::size_t c = 0, r = rand();
    if (c >= (std::size_t)labels.size() / (S_MAX * 100)) {
      c = 0;
      r = rand();
    }
    c++;
    return r;
  }

  void buildIndex() {
    StatusLog log("Build Index");

    timestamp = 0;

    auto setDiscovery = [&](Vertex v) { labels[v].times.first = timestamp++; };
    auto setFinish = [&](Vertex v) { labels[v].times.second = timestamp++; };

    std::function<void(Vertex, DIRECTION)> compute =
        [&](const Vertex v, const DIRECTION dir = FWD) -> void {
      visited.mark(v);
      if (dir == FWD) {
        setDiscovery(v);
      }

      labels[v].set(dir, nextHash() % S_MAX);

      for (std::size_t start = graphs[dir]->beginEdge(v),
                       end = graphs[dir]->endEdge(v);
           start < end; ++start) {
        const Vertex w = graphs[dir]->toVertex[start];
        if (!visited.isMarked(w)) compute(w, dir);
        for (int i = 0; i < NUM_WORDS; ++i) {
          labels[v].bits[dir][i] |= labels[w].bits[dir][i];
        }
      }
      if (dir == FWD) {
        setFinish(v);
      }
    };

    visited.reset();
    for (std::uint32_t v = 0; v < graphs[FWD]->numVertices(); ++v) {
      if (visited.isMarked(v)) continue;
      compute(v, FWD);
    }

    visited.reset();
    for (std::uint32_t v = 0; v < graphs[BWD]->numVertices(); ++v) {
      if (visited.isMarked(v)) continue;
      compute(v, BWD);
    }
  }

  void exportData(const std::string& filename) const {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
      std::cerr << "Error: Could not open file " << filename
                << " for writing.\n";
      return;
    }

    outFile
        << "Vertex,DiscoveredTime,FinishedTime,ForwardLabel,BackwardLabel\n";

    const std::size_t numVertices = graphs[FWD]->numVertices();
    for (Vertex v = 0; v < numVertices; ++v) {
      outFile << v << "," << labels[v].times.first << ","
              << labels[v].times.second << ",";
      for (auto word : labels[v].bits[FWD]) outFile << word << " ";
      outFile << ",";
      for (auto word : labels[v].bits[BWD]) outFile << word << " ";
      outFile << "\n";
    }

    outFile.close();
    std::cout << "Data successfully exported to " << filename << "\n";
  }

  // Query stuff

  template <typename QueryFunction>
  void run_benchmark(const std::string& name, QueryFunction&& queryFn,
                     std::size_t numberOfQueries = 10000) {
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

    auto startTime = std::chrono::high_resolution_clock::now();

    std::size_t positiveCount = 0;
    for (const auto& pr : queries) {
      positiveCount += queryFn(pr.first, pr.second);
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    auto elapsedMillis = std::chrono::duration_cast<std::chrono::milliseconds>(
                             endTime - startTime)
                             .count();

    double timePerQuery = static_cast<double>(elapsedMillis) / numberOfQueries;

    std::cout << name << " Benchmark Results:\n";
    std::cout << "  Number of queries:            " << numberOfQueries << "\n";
    std::cout << "  Total time (ms):              " << elapsedMillis << "\n";
    std::cout << "  Time per query (ms):          " << timePerQuery << "\n";
    std::cout << "  Number of positive queries:   " << positiveCount << "\n";
  }

  // void run_dfs_benchmark(std::size_t numberOfQueries = 10000) {
  //   run_benchmark(
  //       "Simple DFS",
  //       [this](Vertex from, Vertex to) { return dfsQuery(from, to); },
  //       numberOfQueries);
  // }

  // void run_bfs_benchmark(std::size_t numberOfQueries = 10000) {
  //   run_benchmark(
  //       "Simple BFS",
  //       [this](Vertex from, Vertex to) { return bfsQuery(from, to); },
  //       numberOfQueries);
  // }

  // void run_label_bfs_benchmark(std::size_t numberOfQueries = 10000) {
  //   run_benchmark(
  //       "Label-based BFS",
  //       [this](Vertex from, Vertex to) { return bfsQueryPruned(from, to); },
  //       numberOfQueries);
  // }

  void run_label_dfs_benchmark_rec(std::size_t numberOfQueries = 10000) {
    run_benchmark(
        "Label-based DFS-Rec",
        [this](Vertex from, Vertex to) { return dfsRecPruned(from, to); },
        numberOfQueries);
  }

  void run_label_dfs_benchmark_iter(std::size_t numberOfQueries = 10000) {
    run_benchmark(
        "Label-based DFS-Iter",
        [this](Vertex from, Vertex to) { return dfsIterativePruned(from, to); },
        numberOfQueries);
  }

  bool dfsRecPruned(const Vertex from, const Vertex to) {
    assert(graphs[FWD]->isVertex(from));
    assert(graphs[FWD]->isVertex(to));

    if (from == to) [[unlikely]]
      return true;
    visited.reset();

    visited.mark(from);

    return queryRec(from, to);
  }

  bool queryRec(const Vertex from, const Vertex to) {
    if (labels[from].times.second < labels[to].times.second) {
      return false;
    } else if (labels[from].times.first <= labels[to].times.first) {
      return true;
    }

    if ((labels[from].degrees.first == 0) || (labels[to].degrees.second == 0))
      return false;

    for (int i = 0; i < NUM_WORDS; ++i) {
      if ((labels[to].bits[FWD][i] & labels[from].bits[FWD][i]) !=
          labels[to].bits[FWD][i])
        return false;
    }

    for (int i = 0; i < NUM_WORDS; ++i) {
      if ((labels[from].bits[BWD][i] & labels[to].bits[BWD][i]) !=
          labels[from].bits[BWD][i])
        return false;
    }

    for (std::size_t start = graphs[FWD]->beginEdge(from),
                     end = graphs[FWD]->endEdge(from);
         start < end; ++start) {
      PREFETCH(&(graphs[FWD]->toVertex[start + 4]));

      assert(start < graphs[FWD]->toVertex.size());

      const Vertex w = graphs[FWD]->toVertex[start];
      assert(graphs[FWD]->isVertex(w));

      if (!visited.isMarked(w)) {
        visited.mark(w);
        if (queryRec(w, to)) [[unlikely]] {
          return true;
        }
      }
    }

    return false;
  }

  bool dfsIterativePruned(const Vertex from, const Vertex to) {
    assert(graphs[FWD]->isVertex(from));
    assert(graphs[FWD]->isVertex(to));

    if (from == to) [[unlikely]]
      return true;

    visited.reset();
    visited.mark(from);

    index = 0;

    stack[index++] = from;

    while (index > 0) {
      const Vertex u = stack[--index];

      if (labels[u].times.second < labels[to].times.second) {
        continue;
      } else if (labels[u].times.first <= labels[to].times.first) {
        return true;
      }

      if ((labels[u].degrees.first == 0) || (labels[to].degrees.second == 0))
        continue;

      bool pruned = false;

      for (int i = 0; i < NUM_WORDS; ++i) {
        pruned |= ((labels[to].bits[FWD][i] & labels[u].bits[FWD][i]) !=
                   labels[to].bits[FWD][i]);
      }
      if (pruned) continue;

      for (int i = 0; i < NUM_WORDS; ++i) {
        pruned |= ((labels[u].bits[BWD][i] & labels[to].bits[BWD][i]) !=
                   labels[u].bits[BWD][i]);
      }

      if (pruned) continue;

      std::size_t start = graphs[FWD]->beginEdge(u);
      std::size_t end = graphs[FWD]->endEdge(u);

      for (; start < end; ++start) {
        PREFETCH(&(graphs[FWD]->toVertex[start + 4]));

        assert(start < graphs[FWD]->toVertex.size());
        const Vertex w = graphs[FWD]->toVertex[start];
        assert(graphs[FWD]->isVertex(w));

        if (!visited.isMarked(w)) {
          visited.mark(w);

          stack[index++] = w;
        }
      }
    }

    return false;
  }
};
