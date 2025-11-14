/*
 * Licensed under MIT License.
 * Author: Patrick Steil
 */

#include <stdlib.h>

#include "bfl.h"
#include "cmdparser.hpp"
#include "graph.h"
#include "stop_paths.h"

void configure_parser(cli::Parser& parser) {
  parser.set_required<std::string>("i", "input_graph", "Input graph file.");
  parser.set_optional<bool>("s", "show_stats", false,
                            "Show statistics about the computed hub labels.");
  parser.set_optional<bool>("b", "run_benchmark", false,
                            "Runs 100.000 random queries.");
  parser.set_optional<std::string>(
      "o", "output_file", "",
      "Output filename to write labels and times into.");
  parser.set_optional<std::string>("p", "path_file", "", "Path input file.");
};

int main(int argc, char* argv[]) {
  cli::Parser parser(argc, argv, "BFL - Bloom Filter Labeling.");
  configure_parser(parser);
  parser.run_and_exit_if_error();

  const std::string inputFileName = parser.get<std::string>("i");
  const std::string outputFileName = parser.get<std::string>("o");
  const bool showstats = parser.get<bool>("s");
  const bool benchmark = parser.get<bool>("b");
  const std::string pathFileName = parser.get<std::string>("p");
  const std::size_t numQueries = 10000;

  Graph fwdGraph;
  fwdGraph.readDimacs(inputFileName);

  if (showstats) fwdGraph.showStats();

  Graph bwdGraph = fwdGraph.reverseGraph();

  BFL<8> bfl(fwdGraph, bwdGraph);

  bfl.buildIndex();

  if (showstats) bfl.printMemoryConsumption();

  if (outputFileName != "") bfl.exportData(outputFileName);

  if (benchmark) {
    // bfl.run_dfs_benchmark(numQueries);
    // bfl.run_bfs_benchmark(numQueries);
    // bfl.run_label_bfs_benchmark(numQueries);
    bfl.run_label_dfs_benchmark_rec(numQueries);
    bfl.run_label_dfs_benchmark_iter(numQueries);
  }

  if (pathFileName != "") {
    Paths paths;
    paths.loadPathsFromFile(pathFileName);
    if (showstats) paths.showStats();

    if (benchmark) {
      paths.run_benchmark("Path Queries", bfl, numQueries);
    }
  }

  return 0;
}
