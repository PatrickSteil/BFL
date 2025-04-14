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
};

int main(int argc, char *argv[]) {
  cli::Parser parser(argc, argv, "BFL - Bloom Filter Labeling.");
  configure_parser(parser);
  parser.run_and_exit_if_error();

  const std::string inputFileName = parser.get<std::string>("i");
  const bool showstats = parser.get<bool>("s");
  const bool benchmark = parser.get<bool>("b");

  // Graph fwdGraph = Graph::readFromEdgeListFile(inputFileName);
  Graph fwdGraph = Graph::readFromDimacsFile(inputFileName);

  if (showstats) fwdGraph.showStats();

  Graph bwdGraph = fwdGraph.reversed();

  const int s = 256;
  const int d = s * 100;
  BFL<s> bfl(fwdGraph, bwdGraph);

  bfl.computeTimesAndOrder();
  bfl.mergeVertices(d);
  bfl.buildIndex();

  if (benchmark) bfl.run_benchmark();

  return 0;
}
