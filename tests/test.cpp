/**
 * @file test.cpp
 * @brief Parses testcase files and validates the library Solver against the
 *        corresponding solution files (`.sol`).
 *
 * For each testcase file `<name>.txt`, the solution file `<name>.sol` must
 * exist (generate it with the `solve` executable). The results reported by
 * the library Solver for every query are compared against the expected
 * matches from the solution file: match count, matched dictionary ids, and
 * their edit distances.
 *
 * Usage: test_hlp_grep <testcase-file-or-dir>...
 */

#include <hlp_grep/hlp_grep.hpp>

#include "testcase.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

namespace {

/** Expected results for a single query, as read from a solution file. */
struct Expected {
	std::vector<std::size_t> ids;  ///< 0-based dictionary indices, ascending.
	std::vector<int> dists;        ///< Edit distances, parallel to ids.
};

/**
 * @brief Parses a solution file: three lines per query — match count,
 *        1-based indices, and edit distances.
 */
std::vector<Expected> parse_solution(const fs::path &file, std::size_t n_queries) {
	std::ifstream in(file);
	if (!in) {
		std::cerr << "cannot open " << file << " (generate it with `solve`)\n";
		std::exit(1);
	}

	std::vector<Expected> expected;
	expected.reserve(n_queries);
	for (std::size_t q = 0; q < n_queries; ++q) {
		std::string line;
		std::size_t l = 0;
		{
			if (!std::getline(in, line))
				fail(file, static_cast<int>(3 * q + 1),
				     "unexpected EOF while reading match count");
			try {
				l = std::stoul(line);
			} catch (...) {
				fail(file, static_cast<int>(3 * q + 1),
				     "match count is not an integer: " + line);
			}
		}

		Expected e;
		e.ids.reserve(l);
		e.dists.reserve(l);
		{
			if (!std::getline(in, line))
				fail(file, static_cast<int>(3 * q + 2),
				     "unexpected EOF while reading indices");
			std::istringstream iss(line);
			for (std::size_t i = 0; i < l; ++i) {
				std::size_t a = 0;
				if (!(iss >> a))
					fail(file, static_cast<int>(3 * q + 2),
					     "expected " + std::to_string(l) + " indices");
				if (a < 1)
					fail(file, static_cast<int>(3 * q + 2),
					     "index out of range: " + std::to_string(a));
				e.ids.push_back(a - 1); // 1-based -> 0-based
			}
		}
		{
			if (!std::getline(in, line))
				fail(file, static_cast<int>(3 * q + 3),
				     "unexpected EOF while reading distances");
			std::istringstream iss(line);
			for (std::size_t i = 0; i < l; ++i) {
				int d = 0;
				if (!(iss >> d))
					fail(file, static_cast<int>(3 * q + 3),
					     "expected " + std::to_string(l) + " distances");
				e.dists.push_back(d);
			}
		}
		if (!std::is_sorted(e.ids.begin(), e.ids.end()))
			fail(file, static_cast<int>(3 * q + 2), "indices are not sorted");
		expected.push_back(std::move(e));
	}
	return expected;
}

void test_testcase(const fs::path &file, const Testcase &tc,
                   const std::vector<Expected> &expected) {
	const Solver solver(tc.dict, tc.cost);

	for (std::size_t q = 0; q < tc.queries.size(); ++q) {
		const auto &[k, query] = tc.queries[q];
		const auto results = solver.query(query, k);

		if (results.size() != expected[q].ids.size()) {
			std::cerr << file << ": query " << q << " (\"" << query << "\", k=" << k
			          << "): expected " << expected[q].ids.size() << " matches, got "
			          << results.size() << '\n';
			std::exit(1);
		}
		for (std::size_t i = 0; i < results.size(); ++i) {
			if (results[i].id != expected[q].ids[i] ||
			    results[i].dist != expected[q].dists[i]) {
				std::cerr << file << ": query " << q << " (\"" << query << "\", k=" << k
				          << "): mismatch at position " << i << ": expected (id="
				          << expected[q].ids[i] << ", dist=" << expected[q].dists[i]
				          << "), got (id=" << results[i].id << ", dist=" << results[i].dist
				          << ")\n";
				std::exit(1);
			}
		}
	}
	std::cout << "PASS " << file << " (" << tc.dict.size() << " sequences, "
	          << tc.queries.size() << " queries)\n";
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
		const auto expected =
		    parse_solution(solution_path(file), tc.queries.size());
		test_testcase(file, tc, expected);
	}

	std::cout << files.size() << " testcase file(s) passed\n";
	return 0;
}
