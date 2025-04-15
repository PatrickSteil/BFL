/*
 * Licensed under MIT License.
 * Author: Patrick Steil
 */

#pragma once

#include <cstdint>
#include <limits>
#include <vector>

typedef int Vertex;
typedef std::size_t Index;

constexpr int noVertex = int(-1);
constexpr std::size_t noIndex = std::size_t(-1);

enum DIRECTION : bool { FWD, BWD };
