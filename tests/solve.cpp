/**
 * @file solve.cpp
 * @brief Parses testcase files and generates solution files (`.sol`) with the
 *        given solver.
 *
 * For each testcase file `<name>.txt`, a solution file `<name>.sol` is
 * written next to it, in the format described in tests/testcases/README.md:
 * for each query, three lines: the number of matches, their 1-based
 * dictionary indices in ascending order, and their edit distances.
 *
 * Usage: solve <solver-name> <testcase-file-or-dir>...
 */

#include <hlp_grep/hlp_grep.hpp>

#include "testcase.hpp"

#include <fstream>
#include <iostream>

namespace {

void solve_testcase(Solver &solver, const fs::path &file, const Testcase &tc) {
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
	if (argc < 3) {
		std::cerr << "usage: " << argv[0]
		          << " <solver-name> <testcase-file-or-dir>...\n";
		return 1;
	}

	const std::string solver_name = argv[1];
	const auto files = collect_files(argc, argv, 2);
	if (files.empty()) {
		std::cout << "no testcase files found, nothing to do\n";
		return 0;
	}

	for (const auto &file : files) {
		const Testcase tc = parse_testcase(file);
		const auto solver = make_solver(solver_name, tc);
		solve_testcase(*solver, file, tc);
	}

	std::cout << files.size() << " solution file(s) generated\n";
	return 0;
}
