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
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>


using namespace hlp_grep;

namespace {

/**
 * @brief Computes the edit distance between two strings under the given
 *        cost model.
 *
 * Classic dynamic programming algorithm, O(|a| * |b|) time and O(|b|) space.
 * With the default unit-cost model this is the Levenshtein distance.
 */
int edit_distance(const std::string &a, const std::string &b,
                  const CostModel &cost) {
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
	 */
	explicit NaiveSolver(std::vector<std::string> dict,
	                     const CostModel &cost = DEFAULT_COST_MODEL)
	    : dict(std::move(dict)), cost(cost) {}

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
		std::vector<Result> results;
		for (std::size_t i = 0; i < dict.size(); ++i) {
			const int dist = edit_distance(query, dict[i], cost);
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
};

void solve_testcase(const fs::path &file, const Testcase &tc) {
	const NaiveSolver solver(tc.dict, *tc.cost);

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
