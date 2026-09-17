/**
 * @file naive_solver.cpp
 * @brief Parses testcase files and generates solution files (`.sol`) with the
 *        built-in NaiveSolver reference implementation.
 *
 * For each testcase file `<name>.txt`, a solution file `<name>.sol` is
 * written next to it, in the format described in tests/testcases/README.md:
 * for each query, three lines: the number of matches, their 1-based
 * dictionary indices in ascending order, and their edit distances.
 *
 * Usage: solve <testcase-file-or-dir>...
 */

#include <hlp_grep/hlp_grep.hpp>

#include "testcase.hpp"

#include <algorithm>
#include <climits>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>


using namespace hlp_grep;

namespace {

/**
 * @brief Minimum insertion/deletion and substitution costs of the model.
 *
 * Enumerates the full 256-byte domain (out-of-alphabet characters fall back
 * Characters outside the alphabet fall back to index 0 (i.e. alphabet[0]),
 * so scanning the testcase alphabet covers every table entry the DP can
 * ever use. The caller computes this once per query so no per-solver state
 * is needed.
 */
std::pair<int, int> model_mins(const CostModel &cost,
                               const std::string &alphabet) {
	int c_min = INT_MAX, g_min = INT_MAX;
	for (const char ch : alphabet) {
		c_min = std::min({c_min, cost.ins(ch), cost.del(ch)});
		for (const char dh : alphabet)
			g_min = std::min(g_min, cost.consume(ch, dh));
	}
	return {c_min, g_min};
}

/**
 * @brief Computes the edit distance between two strings under the given
 *        cost model, capped at the threshold k.
 *
 * Returns the exact distance when it is <= k and k + 1 otherwise (the
 * .sol format only records matches, so the cap is invisible downstream).
 * When every operation costs >= 0 and every insertion/deletion costs >= 1,
 * any alignment of cost <= k stays within |i - j| <= k / c_min of the main
 * diagonal, so only that Ukkonen band is computed: O((|a| + |b|) * W) time
 * and O(|b|) space with W = k / c_min. Other models fall back to the
 * classic full O(|a| * |b|) dynamic program, which stays exact.
 * With the default unit-cost model this is the Levenshtein distance.
 */
int edit_distance(const std::string &a, const std::string &b,
                  const CostModel &cost, int k, int c_min, int g_min) {
	if (k < 0)
		return k + 1;
	if (g_min < 0 || c_min < 1) {
		std::vector<int> prev(b.size() + 1), curr(b.size() + 1);
		prev[0] = 0;
		for (std::size_t j = 1; j <= b.size(); ++j)
			prev[j] = prev[j - 1] + cost.ins(b[j - 1]);

		for (std::size_t i = 1; i <= a.size(); ++i) {
			curr[0] = prev[0] + cost.del(a[i - 1]);
			for (std::size_t j = 1; j <= b.size(); ++j)
				curr[j] = std::min(
				    {prev[j] + cost.del(a[i - 1]),
				     curr[j - 1] + cost.ins(b[j - 1]),
				     prev[j - 1] + cost.consume(a[i - 1], b[j - 1])});
			std::swap(prev, curr);
		}
		return prev[b.size()];
	}

	const long n = static_cast<long>(a.size());
	const long m = static_cast<long>(b.size());
	const int cap = (k >= INT_MAX - 1) ? INT_MAX : k + 1;
	const long width = static_cast<long>(k) / c_min;
	if (std::labs(n - m) > width)
		return cap;

	// Out-of-band cells cost > k (every path to (i, j) needs >= |i - j|
	// indel steps at >= c_min each), so cap doubles as their value: with
	// nonnegative steps a capped cell can never seed a <= k path that the
	// true optimum would miss.
	std::vector<int> prev(m + 1, cap), curr(m + 1, cap);
	const long row0_hi = std::min(m, width);
	prev[0] = 0;
	for (long j = 1; j <= row0_hi; ++j) {
		const long v = static_cast<long>(prev[j - 1]) + cost.ins(b[j - 1]);
		prev[j] = v > cap ? cap : static_cast<int>(v);
	}

	for (long i = 1; i <= n; ++i) {
		const long lo = std::max(0L, i - width);
		const long hi = std::min(m, i + width);
		const long plo = std::max(0L, i - 1 - width);
		const long phi = std::min(m, i - 1 + width);
		for (long j = lo; j <= hi; ++j) {
			long best = cap;
			if (j >= plo && j <= phi)
				best = std::min(
				    best, static_cast<long>(prev[j]) + cost.del(a[i - 1]));
			if (j - 1 >= plo && j - 1 <= phi)
				best = std::min(best, static_cast<long>(prev[j - 1]) +
				                           cost.consume(a[i - 1], b[j - 1]));
			if (j > lo)
				best = std::min(
				    best, static_cast<long>(curr[j - 1]) + cost.ins(b[j - 1]));
			curr[j] = best > cap ? cap : static_cast<int>(best);
		}
		std::swap(prev, curr);
	}
	return prev[m];
}

/**
 * @brief Reference solver: checks the query against every dictionary
 *        sequence with the classic edit distance dynamic program.
 *
 * Used to generate the solution files against which the library Solver is
 * validated by the test_hlp_grep executable.
 */
class NaiveSolver {
public:
	/**
	 * @brief Constructs the solver over the given sequence dictionary.
	 *
	 * @param dict Dictionary of DNA sequences to search. The position of each
	 *             sequence in this vector defines the `id` reported in Result.
	 * @param cost Cost model defining the costs of the basic edit operations;
	 *             defaults to the unit-cost model. The referenced model must
	 *             outlive the solver.
	 * @param alphabet Characters of the testcase alphabet, in matrix order;
	 *             band parameters are derived from these characters only.
	 */
	explicit NaiveSolver(std::vector<std::string> dict,
	                     const CostModel &cost = DEFAULT_COST_MODEL,
	                     std::string alphabet = "AGCT")
	    : dict(std::move(dict)), cost(cost), alphabet(std::move(alphabet)) {}

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
		// Band parameters from the testcase alphabet, computed once per
		// query (not per pair).
		const auto [c_min, g_min] = model_mins(cost, alphabet);
		std::vector<Result> results;
		for (std::size_t i = 0; i < dict.size(); ++i) {
			const int dist =
			    edit_distance(query, dict[i], cost, k, c_min, g_min);
			if (dist <= k)
				results.push_back({i, dist});
		}
		return results;
	}

private:
	std::vector<std::string> dict; ///< Dictionary of DNA sequences to search.
	/// Cost model used for the edit distance computations; must outlive
	/// the solver.
	const CostModel &cost;
	/// Testcase alphabet band parameters are derived from.
	std::string alphabet;
};

void solve_testcase(const fs::path &file, const Testcase &tc) {
	const NaiveSolver solver(tc.dict, *tc.cost, tc.alphabet);

	const fs::path out_path = solution_path(file);
	std::ofstream out(out_path);
	if (!out) {
		std::cerr << "cannot write " << out_path << '\n';
		std::exit(1);
	}

	for (const auto &[k, query] : tc.queries) {
		const auto results = solver.query(query, k);

		out << results.size() << '\n';
		for (std::size_t i = 0; i < results.size(); ++i)
			out << (i ? " " : "") << results[i].id + 1; // 1-based indices
		out << '\n';
		for (std::size_t i = 0; i < results.size(); ++i)
			out << (i ? " " : "") << results[i].dist;
		out << '\n';
	}
	std::cout << "SOLVED " << file << " -> " << out_path << '\n';
}

} // namespace

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cerr << "usage: " << argv[0] << " <testcase-file-or-dir>...\n";
		return 1;
	}

	const auto files = collect_files(argc, argv);
	if (files.empty()) {
		std::cout << "no testcase files found, nothing to do\n";
		return 0;
	}

	for (const auto &file : files) {
		const Testcase tc = parse_testcase(file);
		solve_testcase(file, tc);
	}

	std::cout << files.size() << " solution file(s) generated\n";
	return 0;
}
