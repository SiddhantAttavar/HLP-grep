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
#include <hlp_grep/cost_model.hpp>
#include <hlp_grep/poa_graph.hpp>
#include <hlp_grep/result.hpp>

#include <string>
#include <stdexcept>
#include <utility>
#include <vector>

/**
 * @brief Edit distance search engine over a sequence dictionary using
 *        heavy-light decomposition on pangenome paths.
 *
 * Indexes the dictionary sequences (represented as paths in a pangenome) and
 * answers edit distance search queries.
 *
 * Usage: construct the solver with the dictionary, then call query() to
 * retrieve all dictionary sequences within a given edit distance threshold
 * of a query string.
 */
class Solver {
public:
	/**
	 * @brief Constructs the solver over the given sequence dictionary.
	 *
	 * @param dict Dictionary of DNA sequences to search. The position of each
	 *             sequence in this vector defines the `id` reported in Result.
	 * @param cost Cost model defining the costs of the basic edit operations.
	 *             Defaults to the unit-cost model (`CostModel{}`).
	 */
	explicit Solver(std::vector<std::string> dict, CostModel cost = {})
	    : dict(std::move(dict)), cost(cost),
	      graph(dict, cost) {
		compressed_paths.reserve(dict.size());
		for (std::size_t i = 0; i < dict.size(); ++i) {
			compressed_paths.push_back(graph.compressed_path(i));
		}
	}

	/**
	 * @brief Finds all dictionary sequences within edit distance k of the query.
	 *
	 * @param query The query DNA sequence.
	 * @param k     Maximum allowed edit distance (threshold).
	 * @return A vector of Result, one per matched dictionary sequence,
	 *         each containing the sequence's index in the dictionary and its
	 *         edit distance to the query. Results are sorted by `id` in
	 *         ascending order.
	 */
	std::vector<Result> query(const std::string &query, int k) const {
		throw std::logic_error("Solver::query() not implemented yet");
	}

private:
	std::vector<std::string> dict; ///< Dictionary of DNA sequences to search.
	CostModel cost;                ///< Cost model used for the edit distance computations.
	POAGraph graph;                ///< Pangenome (POA) graph built from the dictionary.
	/// Compressed (heavy-chain / light-step) representation of each dict path.
	std::vector<POAGraph::CompressedPath> compressed_paths;
};
