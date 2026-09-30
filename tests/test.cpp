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
#include <filesystem>
#include <iostream>
#include <sstream>
#include <vector>


using namespace hlp_grep;

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

void test_query_prefix_sharing() {
	const std::vector<std::string> dict{
	    "ACGT", "ACGA", "A", "ACGC", "ACGTA", "ACTT", "TCGT", "ACGT"};
	const Solver solver(dict);
	const auto results = solver.query("ACGT", 1);
	if (!std::is_sorted(results.begin(), results.end(),
	                    [](const Result &a, const Result &b) {
		                    return a.id < b.id;
	                    })) {
		std::cerr << "prefix-sharing: result ids are not sorted\n";
		std::exit(1);
	}
	const std::vector<Result> expected{
	    {0, 0}, {1, 1}, {3, 1}, {4, 1}, {5, 1}, {6, 1}, {7, 0}};
	if (results.size() != expected.size()) {
		std::cerr << "prefix-sharing: expected " << expected.size()
		          << " matches, got " << results.size() << '\n';
		std::exit(1);
	}
	for (std::size_t i = 0; i < expected.size(); ++i)
		if (results[i].id != expected[i].id ||
		    results[i].dist != expected[i].dist) {
			std::cerr << "prefix-sharing: wrong result at position " << i << '\n';
			std::exit(1);
		}
	std::error_code ec;
	const fs::path cache_file =
	    fs::temp_directory_path(ec) / "hlp-grep-roundtrip.hlpg";
	if (ec) {
		std::cerr << "persistence: no temp directory: " << ec.message()
		          << '\n';
		std::exit(1);
	}
	try {
		solver.save_graph(cache_file.string());
	} catch (const std::exception &e) {
		std::cerr << "persistence: save failed: " << e.what() << '\n';
		std::exit(1);
	}
	const Solver loaded = Solver::load_graph(cache_file.string());
	const auto loaded_results = loaded.query("ACGT", 1);
	if (loaded_results.size() != results.size()) {
		std::cerr << "persistence: loaded solver returned different result count\n";
		std::exit(1);
	}
	for (std::size_t i = 0; i < results.size(); ++i)
		if (loaded_results[i].id != results[i].id ||
		    loaded_results[i].dist != results[i].dist) {
			std::cerr << "persistence: loaded solver result differs\n";
			std::exit(1);
		}
	std::filesystem::remove(cache_file);

	const Solver with_empty({"", "AC"});
	const auto empty_results = with_empty.query("A", 1);
	if (empty_results.size() != 2 || empty_results[0].id != 0 ||
	    empty_results[0].dist != 1 || empty_results[1].id != 1 ||
	    empty_results[1].dist != 1) {
		std::cerr << "prefix-sharing: empty-path scoring failed\n";
		std::exit(1);
	}

	const CostModel weighted(2, 2, 0, 1);
	const Solver weighted_solver({"AA"}, weighted, false);
	const auto over_threshold = weighted_solver.query("AAAA", 2);
	if (!over_threshold.empty()) {
		std::cerr << "prefix-sharing: trailing clipped cells changed the query-end score\n";
		std::exit(1);
	}
	const auto within_threshold = weighted_solver.query("AAA", 2);
	if (within_threshold.size() != 1 || within_threshold[0].id != 0 ||
	    within_threshold[0].dist != 2) {
		std::cerr << "prefix-sharing: two-sided row clipping lost an endpoint match\n";
		std::exit(1);
	}
	std::cout << "PASS query-prefix-sharing\n";
}

void test_lifter_decompose() {
	const std::vector<std::string> dict(4, "ACGTACGTACGT");
	const POAGraph graph(dict, DEFAULT_COST_MODEL, false);
	std::vector<POAGraph::CompressedPath> paths;
	for (std::size_t i = 0; i < dict.size(); ++i)
		paths.push_back(graph.compressed_path(i));
	if (paths[0].steps.size() != 1 ||
	    paths[0].steps[0].type != POAGraph::EdgeType::HEAVY ||
	    paths[0].steps[0].length != 11) {
		std::cerr << "decompose: expected one 11-edge heavy run before split\n";
		std::exit(1);
	}

	BinaryLifter lifter(graph);
	lifter.decompose(paths);
	for (const auto &path : paths) {
		if (path.steps.size() != 3 ||
		    path.steps[0].length != 1 || path.steps[1].length != 2 ||
		    path.steps[2].length != 8) {
			std::cerr << "decompose: length-11 chain was not split as 1+2+8\n";
			std::exit(1);
		}
	}
	std::cout << "PASS lifter-decompose\n";
}

void test_clip_row_two_sided() {
	const POAGraph graph({"A"});
	BinaryLifter lifter(graph, 1);
	std::vector<POAGraph::CompressedPath> paths{graph.compressed_path(0)};
	lifter.decompose(paths);
	lifter.build("ACGT", 1);

	auto range = lifter.window(paths[0].start);
	const auto original_range = range;
	std::vector<int> row(range.second - range.first, 2);
	if (row.size() < 3) {
		std::cerr << "clip-row: test window is too narrow\n";
		std::exit(1);
	}
	row[1] = 1;
	lifter.clip_row(row, range);
	if (row != std::vector<int>{1} || range.first != original_range.first + 1 ||
	    range.second != original_range.first + 2) {
		std::cerr << "clip-row: did not trim both sides of the row\n";
		std::exit(1);
	}
	std::cout << "PASS clip-row-two-sided\n";
}

void test_step_bounded() {
	// Bounded step() must match the full-window computation: stepping
	// with a clipped subrange must give the same clipped result as
	// stepping with the full row. Covers the single-character fast path
	// (no compaction) and the generic sweep (compacted multi-char node).
	for (const bool compact : {false, true}) {
		const POAGraph graph({"ACGTACGT", "ACGTACGA", "ACGTTCGT"},
		                     DEFAULT_COST_MODEL, compact);
		BinaryLifter lifter(graph, 1);
		std::vector<POAGraph::CompressedPath> paths{
		    graph.compressed_path(0)};
		lifter.decompose(paths);
		lifter.build("ACGTACGT", 2);
		bool covered_single = false, covered_multi = false;
		for (std::size_t s = 0; s < graph.num_sequences(); ++s) {
			const auto &path = graph.path(s);
			for (std::size_t i = 1; i < path.size(); ++i) {
				const auto u = path[i - 1];
				const auto v = path[i];
				const bool single = graph.seq(v).size() == 1;
				if (single ? covered_single : covered_multi)
					continue;
				auto full_range = lifter.window(u);
				const std::size_t width =
				    full_range.second - full_range.first;
				if (width < 5)
					continue;
				// Full row with above-threshold margins around a
				// feasible core.
				std::vector<int> full(width, DistMatrix::INF);
				for (std::size_t j = 2; j + 2 < width; ++j)
					full[j] = static_cast<int>(j % 3);
				auto range_a = full_range;
				auto row_a = full;
				lifter.step(u, v, row_a, range_a);
				auto range_b = full_range;
				auto row_b = full;
				lifter.clip_row(row_b, range_b);
				if (row_b.empty())
					continue;
				lifter.step(u, v, row_b, range_b);
				if (row_a != row_b || range_a != range_b) {
					std::cerr << "step-bounded: clipped input "
					             "changed the result (compact="
					          << compact << " single=" << single
					          << ")\n";
					std::exit(1);
				}
				// Fully above-threshold input must stay dead.
				auto range_c = full_range;
				std::vector<int> row_c(width, DistMatrix::INF);
				lifter.step(u, v, row_c, range_c);
				if (!row_c.empty()) {
					std::cerr << "step-bounded: dead input "
					             "produced a live row\n";
					std::exit(1);
				}
				if (single)
					covered_single = true;
				else
					covered_multi = true;
			}
		}
		if (!covered_single) {
			std::cerr << "step-bounded: no single-char edge covered\n";
			std::exit(1);
		}
		if (compact && !covered_multi) {
			std::cerr << "step-bounded: no multi-char edge covered\n";
			std::exit(1);
		}
	}
	std::cout << "PASS step-bounded\n";
}

} // namespace

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cerr << "usage: " << argv[0] << " <testcase-file-or-dir>...\n";
		return 1;
	}

	const auto files = collect_files(argc, argv);
	test_query_prefix_sharing();
	test_lifter_decompose();
	test_clip_row_two_sided();
	test_step_bounded();
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
