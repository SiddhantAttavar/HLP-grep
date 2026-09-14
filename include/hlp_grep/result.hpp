/**
 * @file result.hpp
 * @brief A single approximate match returned by Solver::query().
 */
#pragma once
#include <cstddef>

/**
 * @brief A single approximate match returned by Solver::query().
 */
struct Result {
	std::size_t id; ///< Index of the matched sequence in the dictionary (as passed to the constructor).
	int dist;       ///< Edit distance between the matched sequence and the query string.
};
