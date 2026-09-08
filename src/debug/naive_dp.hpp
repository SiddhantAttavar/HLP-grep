/**
 * @file naive_dp.hpp
 * @brief Naive Solver implementation: full Levenshtein DP of the query
 *        against every dictionary sequence. For debugging and verification.
 */
#pragma once

#include <hlp_grep/hlp_grep.hpp>

#include <algorithm>
#include <numeric>

/**
 * @brief Computes the edit distance between two strings under the given
 *        cost model.
 *
 * Classic dynamic programming algorithm, O(|a| * |b|) time and
 * O(min(|a|, |b|)) space. With the default unit-cost model this is the
 * Levenshtein distance.
 */
inline int edit_distance(const std::string &a, const std::string &b,
                         const CostModel &cost = {}) {
	std::vector<int> prev(b.size() + 1), curr(b.size() + 1);
	for (std::size_t j = 0; j <= b.size(); ++j)
		prev[j] = static_cast<int>(j) * cost.ins;

	for (std::size_t i = 1; i <= a.size(); ++i) {
		curr[0] = static_cast<int>(i) * cost.del;
		for (std::size_t j = 1; j <= b.size(); ++j)
			curr[j] = std::min({prev[j] + cost.del, curr[j - 1] + cost.ins,
			                    prev[j - 1] + cost.match(a[i - 1], b[j - 1])});
		std::swap(prev, curr);
	}
	return prev[b.size()];
}

/**
 * @brief Naive solver: checks the query against every dictionary
 *        sequence with the classic edit distance dynamic program.
 */
class NaiveSolver : public Solver {
public:
	/**
	 * @brief Constructs the solver over the given sequence dictionary.
	 *
	 * @param dict Dictionary of DNA sequences to search. The position of each
	 *             sequence in this vector defines the `id` reported in Result.
	 * @param cost Cost model defining the costs of the basic edit operations.
	 *             Defaults to the unit-cost model (`CostModel{}`).
	 */
	explicit NaiveSolver(std::vector<std::string> dict, CostModel cost = {})
	    : Solver(std::move(dict), cost) {}

	std::vector<Result> query(const std::string &query, int k) override {
		std::vector<Result> results;
		for (std::size_t i = 0; i < dict.size(); ++i) {
			const int dist = edit_distance(query, dict[i], cost);
			if (dist <= k)
				results.push_back({i, dist});
		}
		return results;
	}
};
