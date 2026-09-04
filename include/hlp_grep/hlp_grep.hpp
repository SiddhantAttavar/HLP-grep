/**
 * @file hlp_grep.hpp
 * @brief Edit distance search queries on a DNA sequence dictionary using
 *        heavy-light decomposition on pangenome paths.
 *
 * HLP-grep indexes a dictionary of DNA sequences (represented as paths in a
 * pangenome) and answers edit distance search queries: given a query string
 * and a threshold k, it finds all dictionary sequences within edit distance
 * k of the query string.
 */
#pragma once
#include <vector>
#include <string>

/**
 * @brief A single approximate match returned by Solver::query().
 */
struct Result {
	std::size_t id;   ///< Index of the matched sequence in the dictionary (as passed to the constructor).
	int dist; ///< Edit distance between the matched sequence and the query string.
};

/**
 * @brief Index and query engine for edit distance search over a sequence dictionary.
 *
 * Usage: construct a concrete Solver implementation with the dictionary, then
 * call query() to retrieve all dictionary sequences within a given edit
 * distance threshold of a query string.
 */
class Solver {
public:
	virtual ~Solver() = default;

	/**
	 * @brief Finds all dictionary sequences within edit distance k of the query.
	 *
	 * @param query The query DNA sequence.
	 * @param k     Maximum allowed edit distance (threshold).
	 * @return A vector of Result, one per matched dictionary sequence,
	 *         each containing the sequence's index in the dictionary and its
	 *         edit distance to the query.
	 */
	virtual std::vector<Result> query(const std::string &query, int k) = 0;
};
