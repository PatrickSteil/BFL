/*
 * Licensed under MIT License.
 * Author: Patrick Steil
 */

#include <stdlib.h>

#include "bfl.h"
#include "cmdparser.hpp"
#include "graph.h"

void configure_parser(cli::Parser &parser) {
  parser.set_required<std::string>("i", "input_graph", "Input graph file.");
  parser.set_optional<bool>("s", "show_stats", false,
                            "Show statistics about the computed hub labels.");
  parser.set_optional<bool>("b", "run_benchmark", false,
                            "Runs 100.000 random queries.");
  parser.set_optional<std::string>(
      "o", "output_file", "",
      "Output filename to write labels and times into.");
};

int main(int argc, char *argv[]) {
  cli::Parser parser(argc, argv, "BFL - Bloom Filter Labeling.");
  configure_parser(parser);
  parser.run_and_exit_if_error();

  const std::string inputFileName = parser.get<std::string>("i");
  const std::string outputFileName = parser.get<std::string>("o");
  const bool showstats = parser.get<bool>("s");
  const bool benchmark = parser.get<bool>("b");

  // Graph fwdGraph = Graph::readFromEdgeListFile(inputFileName);
  Graph fwdGraph;
  fwdGraph.readDimacs(inputFileName);

  if (showstats) fwdGraph.showStats();

  Graph bwdGraph = fwdGraph.reverseGraph();

  BFL<6> bfl(fwdGraph, bwdGraph);

  bfl.buildIndex();

  if (showstats) bfl.printMemoryConsumption();

  if (outputFileName != "") bfl.exportData(outputFileName);

  if (benchmark) {
    // bfl.run_dfs_benchmark(10000);
    // bfl.run_bfs_benchmark(10000);
    // bfl.run_label_bfs_benchmark(10000);
    bfl.run_label_dfs_benchmark_rec(10000);
    bfl.run_label_dfs_benchmark_iter(10000);
  }

  return 0;
}
