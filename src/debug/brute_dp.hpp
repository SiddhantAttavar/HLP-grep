/**
 * @file brute_dp.hpp
 * @brief Brute-force Solver implementation: full Levenshtein DP of the query
 *        against every dictionary sequence. For debugging and verification.
 */
#pragma once

#include <hlp_grep/hlp_grep.hpp>

#include <algorithm>
#include <numeric>

/**
 * @brief Brute-force solver: checks the query against every dictionary
 *        sequence with the classic edit distance dynamic program.
 */
class BruteSolver : public Solver {
public:
	/**
	 * @brief Constructs the solver over the given sequence dictionary.
	 *
	 * @param dict Dictionary of DNA sequences to search. The position of each
	 *             sequence in this vector defines the `id` reported in Result.
	 */
	explicit BruteSolver(const std::vector<std::string> &dict) : dict_(dict) {}

	std::vector<Result> query(const std::string &query, int k) override {
		std::vector<Result> results;
		for (std::size_t i = 0; i < dict_.size(); ++i) {
			const int dist = edit_distance(query, dict_[i]);
			if (dist <= k)
				results.push_back({i, dist});
		}
		return results;
	}

private:
	std::vector<std::string> dict_;

	/**
	 * @brief Computes the Levenshtein edit distance between two strings.
	 *
	 * Classic dynamic programming algorithm, O(|a| * |b|) time and
	 * O(min(|a|, |b|)) space.
	 */
	static int edit_distance(const std::string &a, const std::string &b) {
		if (a.size() < b.size())
			return edit_distance(b, a);

		std::vector<int> prev(b.size() + 1), curr(b.size() + 1);
		std::iota(prev.begin(), prev.end(), 0);

		for (std::size_t i = 1; i <= a.size(); ++i) {
			curr[0] = static_cast<int>(i);
			for (std::size_t j = 1; j <= b.size(); ++j)
				curr[j] = std::min({prev[j] + 1, curr[j - 1] + 1,
				                    prev[j - 1] + (a[i - 1] != b[j - 1])});
			std::swap(prev, curr);
		}
		return prev[b.size()];
	}
};
