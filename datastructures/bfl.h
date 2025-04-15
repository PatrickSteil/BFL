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

template <int K = 5>
struct BFL {
  static constexpr int NUM_WORDS = K;
  struct Label {
    std::array<std::array<std::uint32_t, NUM_WORDS>, 2> bits;
    std::array<std::uint32_t, 2> times = {0, 0};
    std::array<std::uint16_t, 2> degrees = {0, 0};
    std::array<Vertex*, 2> edges = {nullptr, nullptr};
    std::uint32_t lastSeen{0};

    Label() {
      bits[FWD].fill(0);
      bits[BWD].fill(0);
    }

    void set(const DIRECTION dir, std::size_t pos) {
      std::size_t word = pos / 32;
      std::size_t offset = pos % 32;
      bits[dir][word] |= (1ULL << offset);
    }

    bool isMarked(std::uint32_t timestamp) { return lastSeen == timestamp; }

    void mark(std::uint32_t timestamp) { lastSeen = timestamp; }
  };

  std::vector<Label> labels;

  std::vector<Vertex> stack;
  std::size_t timestamp;
  std::size_t index;
  std::uint32_t timer;

  BFL(Graph& fwdGraph, Graph& bwdGraph)
      : labels(fwdGraph.numVertices(), Label()),
        stack(fwdGraph.numVertices(), 0),
        timestamp(0),
        index(0),
        timer(0) {
    for (Vertex v = 0; v < fwdGraph.numVertices(); ++v) {
      labels[v].degrees[FWD] = fwdGraph.degree(v);
      labels[v].degrees[BWD] = bwdGraph.degree(v);

      labels[v].edges[FWD] = &(fwdGraph.toVertex[fwdGraph.beginEdge(v)]);
      labels[v].edges[BWD] = &(bwdGraph.toVertex[bwdGraph.beginEdge(v)]);
    }
  }

  void printMemoryConsumption() const {
    std::size_t numVertices = labels.size();

    std::size_t labelsMemory = 2 * numVertices * sizeof(Label::bits);
    std::size_t timesMemory = 2 * numVertices * sizeof(Label::times);
    std::size_t totalMemory = labelsMemory + timesMemory;

    std::cout << "Memory Consumption:\n";
    std::cout << "  Labels memory: "
              << static_cast<double>(labelsMemory) / (1024.0 * 1024.0)
              << " mb\n";
    std::cout << "  Discovery/Finish times memory: "
              << static_cast<double>(timesMemory) / (1024.0 * 1024.0)
              << " mb\n";
    std::cout << "  Total memory: "
              << static_cast<double>(totalMemory) / (1024.0 * 1024.0)
              << " mb\n";
  }

  std::size_t nextHash() {
    static std::size_t c = 0, r = rand();
    if (c >=
        (std::size_t)labels.size() / (K * (sizeof(std::uint32_t) << 3) * 10)) {
      c = 0;
      r = rand();
    }
    c++;
    return r;
  }

  void resetTimer() {
    ++timer;

    if (timer == 0) {
      timer = 1;

      for (auto& n : labels) n.mark(0);
    }
  }

  void buildIndex() {
    StatusLog log("Build Index");

    timestamp = 0;

    auto setDiscovery = [&](Vertex v) { labels[v].times[FWD] = timestamp++; };
    auto setFinish = [&](Vertex v) { labels[v].times[BWD] = timestamp++; };

    std::function<void(Vertex, DIRECTION)> compute =
        [&](const Vertex v, const DIRECTION dir = FWD) -> void {
      labels[v].mark(timer);
      if (dir == FWD) {
        setDiscovery(v);
      }

      labels[v].set(dir, nextHash() % (K * (sizeof(std::uint32_t) << 3)));

      for (int i = 0; i < labels[v].degrees[dir]; ++i) {
        const Vertex w = labels[v].edges[dir][i];
        if (!labels[w].isMarked(timer)) compute(w, dir);

        for (int i = 0; i < NUM_WORDS; ++i) {
          labels[v].bits[dir][i] |= labels[w].bits[dir][i];
        }
      }
      if (dir == FWD) {
        setFinish(v);
      }
    };

    resetTimer();
    for (int v = 0; v < labels.size(); ++v) {
      if (labels[v].isMarked(timer)) continue;
      compute(v, FWD);
    }

    resetTimer();
    for (int v = 0; v < labels.size(); ++v) {
      if (labels[v].isMarked(timer)) continue;
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

    const std::size_t numVertices = labels.size();
    for (Vertex v = 0; v < numVertices; ++v) {
      outFile << v << "," << labels[v].times[FWD] << "," << labels[v].times[BWD]
              << ",";
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

    std::size_t n = labels.size();
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
    assert(from < labels.size());
    assert(to < labels.size());

    if (from == to) [[unlikely]]
      return true;
    resetTimer();
    labels[from].mark(timer);

    return queryRec(labels[from], labels[to]);
  }

  bool queryRec(Label& from, Label& to) {
    if (from.times[BWD] < to.times[BWD]) {
      return false;
    } else if (from.times[FWD] <= to.times[FWD]) {
      return true;
    }

    if ((from.degrees[FWD] == 0) || (to.degrees[BWD] == 0)) return false;

    for (int i = 0; i < NUM_WORDS; ++i) {
      if ((to.bits[FWD][i] & from.bits[FWD][i]) != to.bits[FWD][i])
        return false;
    }

    for (int i = 0; i < NUM_WORDS; ++i) {
      if ((from.bits[BWD][i] & to.bits[BWD][i]) != from.bits[BWD][i])
        return false;
    }

    for (int i = 0; i < from.degrees[FWD]; ++i) {
      const Vertex w = from.edges[FWD][i];
      if (!labels[w].isMarked(timer)) {
        labels[w].mark(timer);
        if (queryRec(labels[w], to)) [[unlikely]] {
          return true;
        }
      }
    }

    return false;
  }

  bool dfsIterativePruned(const Vertex from, const Vertex to) {
    assert(from < labels.size());
    assert(to < labels.size());

    if (from == to) [[unlikely]]
      return true;

    resetTimer();
    labels[from].mark(timer);

    index = 0;

    stack[index++] = from;

    while (index > 0) {
      const Vertex u = stack[--index];

      if (labels[u].times[BWD] < labels[to].times[BWD]) {
        continue;
      } else if (labels[u].times[FWD] <= labels[to].times[FWD]) [[unlikely]] {
        return true;
      }

      if ((labels[u].degrees[FWD] == 0) || (labels[to].degrees[BWD] == 0))
        continue;

      bool pruned = false;

      for (int i = 0; i < NUM_WORDS; ++i) {
        pruned |= ((labels[to].bits[FWD][i] & labels[u].bits[FWD][i]) !=
                   labels[to].bits[FWD][i]);
      }
      if (pruned) [[unlikely]]
        continue;

      for (int i = 0; i < NUM_WORDS; ++i) {
        pruned |= ((labels[u].bits[BWD][i] & labels[to].bits[BWD][i]) !=
                   labels[u].bits[BWD][i]);
      }

      if (pruned) [[unlikely]]
        continue;

      for (int i = 0; i < labels[u].degrees[FWD]; ++i) {
        if (4 < labels[u].degrees[FWD]) {
          PREFETCH(&labels[labels[u].edges[FWD][i + 4]]);
        }
        const Vertex w = labels[u].edges[FWD][i];
        if (!labels[w].isMarked(timer)) {
          labels[w].mark(timer);

          stack[index++] = w;
        }
      }
    }

    return false;
  }
};
