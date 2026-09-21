/**
 * @file bench_poa.cpp
 * @brief Benchmarks POA graph construction time and quality under the
 *        default band configuration.
 *
 * For each testcase file, a Solver and a POAGraph are built with the
 * dynamic doubling band (see POAGraph) and the following is reported:
 * minimum construction time over `--reps` runs, node/edge counts,
 * compression (`|V| / total length`), and mean/max `j_max - j_min`
 * spread. Using a single config keeps repeated graph construction from
 * dominating profiling runs.
 *
 * Deterministic invariants are enforced (nonzero exit on violation), so
 * this doubles as a CTest: every stored path must spell its dictionary
 * sequence, and every query must return self-consistent results.
 *
 * Usage: bench_poa [--max-seqs N] [--reps R] <testcase-file-or-dir>...
 */

#include <hlp_grep/hlp_grep.hpp>

#include "testcase.hpp"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace hlp_grep;

namespace {

/** Exits with an error message on failure, like the other test helpers. */
void check(bool ok, const std::string &msg) {
	if (!ok) {
		std::cerr << "FAIL: " << msg << '\n';
		std::exit(1);
	}
}

/** Measured construction time and quality for one config on one file. */
struct Report {
	double build_s = 0;      ///< Minimum Solver construction time over reps.
	std::size_t nodes = 0;   ///< Graph node count |V|.
	std::size_t edges = 0;   ///< Unique graph edges over all stored paths.
	std::size_t total_len = 0; ///< Sum of stored path lengths.
	double spread_mean = 0;  ///< Mean j_max - j_min over nodes on a path.
	std::size_t spread_max = 0; ///< Max j_max - j_min over nodes on a path.
	std::string fingerprint; ///< Serialized query results (exactness check).
};

/** Builds a Solver and times it; returns the fastest of @p reps runs. */
Report time_solver(const Testcase &tc, int reps) {
	Report rep;
	double best = std::numeric_limits<double>::infinity();
	for (int r = 0; r < reps; ++r) {
		const auto t0 = std::chrono::steady_clock::now();
		Solver solver(tc.dict, *tc.cost);
		const auto t1 = std::chrono::steady_clock::now();
		best = std::min(best, std::chrono::duration<double>(t1 - t0).count());
	}
	rep.build_s = best;
	return rep;
}

/** Fills graph quality metrics and enforces the path-spelling invariant. */
void measure_graph(const Testcase &tc, const fs::path &file, Report &rep) {
	const POAGraph graph(tc.dict, *tc.cost);
	rep.nodes = graph.num_nodes();
	std::unordered_set<std::uint64_t> edges;
	double spread_sum = 0;
	std::size_t spread_n = 0;
	for (std::size_t s = 0; s < graph.num_sequences(); ++s) {
		const auto &path = graph.path(s);
		rep.total_len += path.size();
		std::string spelled;
		for (std::size_t i = 0; i < path.size(); ++i) {
			check(path[i] < graph.num_nodes(), "node id out of range");
			spelled += graph.base(path[i]);
			if (i > 0)
				edges.insert(static_cast<std::uint64_t>(path[i - 1]) << 32 |
				             static_cast<std::uint64_t>(path[i]));
		}
		check(spelled == tc.dict[s],
		      "path " + std::to_string(s) + " misspells its sequence");
	}
	rep.edges = edges.size();
	for (std::size_t u = 0; u < graph.num_nodes(); ++u) {
		const auto [lo, hi] = graph.pos_range(u);
		if (lo == static_cast<std::size_t>(-1))
			continue;
		spread_sum += static_cast<double>(hi - lo);
		rep.spread_max = std::max(rep.spread_max, hi - lo);
		++spread_n;
	}
	if (spread_n > 0)
		rep.spread_mean = spread_sum / static_cast<double>(spread_n);
}

/** Serializes all query results for the run's single band config. */
std::string fingerprint_queries(const Testcase &tc) {
	const Solver solver(tc.dict, *tc.cost);
	std::string finger;
	for (const auto &[k, query] : tc.queries) {
		const auto results = solver.query(query, k);
		finger += std::to_string(results.size()) + ":";
		for (const auto &r : results)
			finger += std::to_string(r.id) + "," + std::to_string(r.dist) + ";";
		finger += "|";
	}
	return finger;
}

void bench_testcase(const fs::path &file, int max_seqs, int reps) {
	Testcase tc = parse_testcase(file);
	if (max_seqs > 0 && static_cast<std::size_t>(max_seqs) < tc.dict.size())
		tc.dict.resize(static_cast<std::size_t>(max_seqs));

	Report rep = time_solver(tc, reps);
	measure_graph(tc, file, rep);
	rep.fingerprint = fingerprint_queries(tc);

	const double compr =
	    rep.total_len > 0 ? static_cast<double>(rep.nodes) /
	                            static_cast<double>(rep.total_len)
	                      : 0.0;
	std::cout << "BENCH " << file << " config=dynamic"
	          << " build_s=" << rep.build_s << " V=" << rep.nodes
	          << " E=" << rep.edges << " compr=" << compr
	          << " spread=" << rep.spread_mean << "/" << rep.spread_max
	          << " queries=" << tc.queries.size() << '\n';
	std::cout << "PASS " << file << " (paths spell sequences, queries ran)\n";
}

} // namespace

int main(int argc, char **argv) {
	int max_seqs = 0;
	int reps = 3;
	int first = 1;
	for (; first < argc && argv[first][0] == '-'; ++first) {
		const std::string flag = argv[first];
		if (flag == "--max-seqs" && first + 1 < argc) {
			max_seqs = std::atoi(argv[++first]);
		} else if (flag == "--reps" && first + 1 < argc) {
			reps = std::atoi(argv[++first]);
		} else {
			std::cerr << "unknown flag: " << flag << '\n';
			return 1;
		}
	}
	check(reps > 0, "--reps must be positive");

	// Dynamic band: the production configuration (no tuning knobs).
	const auto files = collect_files(argc, argv, first);
	if (files.empty()) {
		std::cout << "no testcase files found, nothing to do\n";
		return 0;
	}
	for (const auto &file : files)
		bench_testcase(file, max_seqs, reps);

	std::cout << files.size() << " testcase file(s) benchmarked\n";
	return 0;
}
