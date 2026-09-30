/**
 * @file load.cpp
 * @brief Loads a saved POAGraph file and answers a testcase's queries.
 *
 * For the testcase file `<name>.txt` and the graph file written by `save`,
 * a solution file `<name>.sol` is written (or the path given with `--out`),
 * in the same format as `solve`: for each query, three lines holding the
 * number of matches, their 1-based dictionary indices in ascending order,
 * and their edit distances.
 *
 * The dictionary build step is skipped: the solver is reconstructed from
 * the saved graph via Solver::load_graph() instead of the dictionary
 * sequences. The graph file must have been saved from the testcase's
 * dictionary (its stored cost model must also match the testcase's).
 *
 * Usage: load <testcase-file> <graph-file> [--out <sol-file>]
 */

#include <hlp_grep/hlp_grep.hpp>

#include "testcase.hpp"

#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace hlp_grep;

namespace {

[[noreturn]] void usage(const char *argv0) {
	std::cerr << "usage: " << argv0
	          << " <testcase-file> <graph-file> [--out <sol-file>]\n";
	std::exit(1);
}

Solver load_solver(const Testcase &tc, const std::string &graph_file) {
	try {
		return Solver::load_graph(graph_file, tc.cost);
	} catch (const std::exception &e) {
		std::cerr << "load failed: " << e.what() << '\n';
		std::exit(1);
	}
}

} // namespace

int main(int argc, char **argv) {
	std::string testcase;
	std::string graph_file;
	std::string out_path;
	for (int i = 1; i < argc; ++i) {
		const std::string arg = argv[i];
		if (arg == "--out") {
			if (i + 1 >= argc)
				usage(argv[0]);
			out_path = argv[++i];
		} else if (!arg.empty() && arg[0] != '-' && testcase.empty()) {
			testcase = arg;
		} else if (!arg.empty() && arg[0] != '-' && graph_file.empty()) {
			graph_file = arg;
		} else {
			std::cerr << "unexpected argument: " << arg << '\n';
			usage(argv[0]);
		}
	}
	if (testcase.empty()) {
		std::cerr << "missing testcase file\n";
		usage(argv[0]);
	}
	if (graph_file.empty()) {
		std::cerr << "missing graph file\n";
		usage(argv[0]);
	}
	const fs::path file(testcase);
	if (!fs::is_regular_file(file)) {
		std::cerr << "not a testcase file: " << file << '\n';
		return 1;
	}
	if (out_path.empty())
		out_path = solution_path(file).string();

	const Testcase tc = parse_testcase(file);
	// Skips the dictionary build: the solver is reconstructed from the
	// saved graph, and only the testcase queries are answered.
	Solver solver = load_solver(tc, graph_file);

	std::ofstream out(out_path);
	if (!out) {
		std::cerr << "cannot write " << out_path << '\n';
		return 1;
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
	return 0;
}
