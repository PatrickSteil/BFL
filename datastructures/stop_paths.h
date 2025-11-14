/*
 * Licensed under MIT License.
 * Author: Patrick Steil
 */

#pragma once

#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "bfl.h"
#include "status_log.h"
#include "types.h"

struct Paths {
  std::vector<std::vector<Vertex>> paths;
  mutable std::mt19937 rng{std::random_device{}()};

  Paths(const std::size_t numPaths = 0) : paths(numPaths) {}

  const std::size_t numPaths() const { return paths.size(); }

  const std::vector<Vertex>& getPath(const std::size_t pathId) const {
    assert(pathId < numPaths());

    return paths[pathId];
  }

  void loadPathsFromFile(const std::string& filename) {
    StatusLog log("Reading paths from file");

    paths.clear();

    std::ifstream in(filename);
    if (!in) throw std::runtime_error("Cannot open file: " + filename);

    std::string line;

    if (!std::getline(in, line))
      throw std::runtime_error("File is empty: " + filename);

    {
      std::istringstream iss(line);
      char P;
      std::size_t numPaths = 0;

      if (!(iss >> P >> numPaths) || P != 'P') {
        throw std::runtime_error("Expected first line: 'P <numPaths>'");
      }

      paths.resize(numPaths);
    }

    std::size_t pathIndex = 0;

    while (std::getline(in, line)) {
      if (line.empty()) continue;

      if (pathIndex >= paths.size()) {
        throw std::runtime_error(
            "More path lines in file than declared in header.");
      }

      std::istringstream iss(line);
      Vertex v;

      while (iss >> v) {
        paths[pathIndex].push_back(v - 1);
      }

      ++pathIndex;
    }

    if (pathIndex != paths.size()) {
      std::ostringstream err;
      err << "File declared " << paths.size() << " paths but only contained "
          << pathIndex << ".";
      throw std::runtime_error(err.str());
    }
  }

  void showStats() const {
    const std::size_t numPaths = paths.size();
    if (numPaths == 0) {
      std::cout << "No paths loaded." << std::endl;
      return;
    }

    std::size_t totalVertices = 0;
    std::size_t minLen = std::numeric_limits<std::size_t>::max();
    std::size_t maxLen = 0;

    for (const auto& p : paths) {
      const std::size_t len = p.size();
      totalVertices += len;
      if (len < minLen) minLen = len;
      if (len > maxLen) maxLen = len;
    }

    const double avgLen = static_cast<double>(totalVertices) / numPaths;

    std::cout << "Paths:             " << numPaths << "\n";
    std::cout << "Total vertices:    " << totalVertices << "\n";
    std::cout << "Average length:    " << avgLen << "\n";
    std::cout << "Min length:        " << minLen << "\n";
    std::cout << "Max length:        " << maxLen << "\n";
  }

  template <int K>
  int query(std::size_t fromPathId, std::size_t toPathId, BFL<K>& bfl) const {
    assert(fromPathId < paths.size() && toPathId < paths.size());
    const auto& fromPath = paths[fromPathId];
    const auto& toPath = paths[toPathId];

    if (fromPath.empty() || toPath.empty()) [[unlikely]]
      return -1;

    /* std::uniform_int_distribution<std::size_t> dist(0, fromPath.size() - 1);
     */
    /* Vertex source = fromPath[dist(rng)]; */
    Vertex source = fromPath[fromPath.size() / 2];

    int best = -1;
    int L = 0;
    int R = (int)toPath.size() - 1;

    while (L <= R) {
      int mid = (L + R) / 2;
      Vertex v = toPath[mid];

      if (bfl.dfsIterativePruned(source, v)) {
        best = mid;
        R = mid - 1;
      } else {
        L = mid + 1;
      }
    }

    return best;
  }

  template <int K>
  void run_benchmark(const std::string& name, BFL<K>& bfl,
                     std::size_t numberOfQueries = 10000) const {
    std::size_t n = paths.size();
    if (n == 0) {
      std::cout << "No paths loaded, benchmark aborted.\n";
      return;
    }

    std::mt19937 gen(42);
    std::uniform_int_distribution<std::size_t> dis(0, n - 1);

    std::vector<std::pair<std::size_t, std::size_t>> queries;
    queries.reserve(numberOfQueries);

    for (std::size_t i = 0; i < numberOfQueries; i++) {
      std::size_t fromId = dis(gen);
      std::size_t toId = dis(gen);
      queries.emplace_back(fromId, toId);
    }

    auto start = std::chrono::high_resolution_clock::now();

    std::size_t positiveCount = 0;
    for (const auto& pr : queries) {
      positiveCount += (query(pr.first, pr.second, bfl) != -1);
    }

    auto end = std::chrono::high_resolution_clock::now();

    auto elapsedMs =
        std::chrono::duration<double, std::milli>(end - start).count();
    auto elapsedNs =
        std::chrono::duration<double, std::nano>(end - start).count();

    double timePerQueryMs = elapsedMs / numberOfQueries;
    double timePerQueryUs = (elapsedMs * 1000.0) / numberOfQueries;
    double timePerQueryNs = elapsedNs / numberOfQueries;

    std::cout << name << " Benchmark Results:\n";
    std::cout << "  Number of queries:          " << numberOfQueries << "\n";
    std::cout << "  Total time:                 " << elapsedMs << " ms\n";
    std::cout << "  Time per query:             " << timePerQueryMs
              << " ms  |  " << timePerQueryUs << " µs  |  " << timePerQueryNs
              << " ns\n";
    std::cout << "  Positive query results:     " << positiveCount << "\n";
  }
};
