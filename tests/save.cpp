/**
 * @file save.cpp
 * @brief Builds the solver for a testcase dictionary and saves its POAGraph.
 *
 * For the testcase file `<name>.txt`, a graph file `<name>.hlpg` is written
 * (or the path given with `--out`), holding the complete finalized graph:
 * node labels, edges, paths, topological order, and heavy-chain metadata.
 * The `load` executable reconstructs a solver from that file without
 * rebuilding the dictionary.
 *
 * Usage: save <testcase-file> [--out <graph-file>] [--no-compact]
 */

#include <hlp_grep/hlp_grep.hpp>

#include "testcase.hpp"

#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>

using namespace hlp_grep;

namespace {

[[noreturn]] void usage(const char *argv0) {
	std::cerr << "usage: " << argv0
	          << " <testcase-file> [--out <graph-file>] [--no-compact]\n";
	std::exit(1);
}

} // namespace

int main(int argc, char **argv) {
	std::string testcase;
	std::string out_path;
	bool compact_nodes = true;
	for (int i = 1; i < argc; ++i) {
		const std::string arg = argv[i];
		if (arg == "--out") {
			if (i + 1 >= argc)
				usage(argv[0]);
			out_path = argv[++i];
		} else if (arg == "--no-compact") {
			compact_nodes = false; // POAGraph run compaction off
		} else if (!arg.empty() && arg[0] != '-' && testcase.empty()) {
			testcase = arg;
		} else {
			std::cerr << "unexpected argument: " << arg << '\n';
			usage(argv[0]);
		}
	}
	if (testcase.empty()) {
		std::cerr << "missing testcase file\n";
		usage(argv[0]);
	}
	const fs::path file(testcase);
	if (!fs::is_regular_file(file)) {
		std::cerr << "not a testcase file: " << file << '\n';
		return 1;
	}
	if (out_path.empty()) {
		fs::path out = file;
		out.replace_extension(".hlpg");
		out_path = out.string();
	}

	const Testcase tc = parse_testcase(file);
	try {
		const Solver solver(tc.dict, tc.cost, compact_nodes);
		solver.save_graph(out_path);
	} catch (const std::exception &e) {
		std::cerr << "save failed: " << e.what() << '\n';
		return 1;
	}
	std::cout << "SAVED " << file << " -> " << out_path
	          << " (dict=" << tc.dict.size() << ")\n";
	return 0;
}
