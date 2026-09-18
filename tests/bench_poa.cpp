/**
 * @file bench_poa.cpp
 * @brief Benchmarks POA graph construction time and quality across build
 *        band configurations.
 *
 * For each testcase file, a Solver and a POAGraph are built under three
 * band configs — full (exact reference: the band covers every column),
 * default (`b = 10, f = 0.01`) and narrow (`b = 2, f = 0`) — and the
 * following is reported per config: minimum construction time over
 * `--reps` runs, node/edge counts, compression (`|V| / total length`),
 * and mean/max `j_max - j_min` spread.
 *
 * Deterministic invariants are enforced (nonzero exit on violation), so
 * this doubles as a CTest: every stored path must spell its dictionary
 * sequence, and every query must return identical results under all
 * configs (build banding may change graph sharing, never search results).
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

/** One build band configuration under test. */
struct BandConfig {
	std::string name; ///< Label used in the report.
	int base;         ///< Band base half-width `b`.
	double slope;     ///< Band slope `f`.
};

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
Report time_solver(const Testcase &tc, const BandConfig &cfg, int reps) {
	Report rep;
	double best = std::numeric_limits<double>::infinity();
	for (int r = 0; r < reps; ++r) {
		const auto t0 = std::chrono::steady_clock::now();
		Solver solver(tc.dict, *tc.cost, cfg.base, cfg.slope);
		const auto t1 = std::chrono::steady_clock::now();
		best = std::min(best, std::chrono::duration<double>(t1 - t0).count());
	}
	rep.build_s = best;
	return rep;
}

/** Fills graph quality metrics and enforces the path-spelling invariant. */
void measure_graph(const Testcase &tc, const BandConfig &cfg, const fs::path &file,
                   Report &rep) {
	const POAGraph graph(tc.dict, *tc.cost, cfg.base, cfg.slope);
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

/** Serializes all query results; must agree across band configs. */
std::string fingerprint_queries(const Testcase &tc, const BandConfig &cfg) {
	const Solver solver(tc.dict, *tc.cost, cfg.base, cfg.slope);
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

void bench_testcase(const fs::path &file, const std::vector<BandConfig> &cfgs,
                    int max_seqs, int reps) {
	Testcase tc = parse_testcase(file);
	if (max_seqs > 0 && static_cast<std::size_t>(max_seqs) < tc.dict.size())
		tc.dict.resize(static_cast<std::size_t>(max_seqs));

	std::vector<Report> reps_out;
	for (const auto &cfg : cfgs) {
		Report rep = time_solver(tc, cfg, reps);
		measure_graph(tc, cfg, file, rep);
		rep.fingerprint = fingerprint_queries(tc, cfg);
		reps_out.push_back(std::move(rep));
	}
	for (std::size_t c = 1; c < cfgs.size(); ++c)
		check(reps_out[c].fingerprint == reps_out[0].fingerprint,
		      "query results differ under band " + cfgs[c].name);

	for (std::size_t c = 0; c < cfgs.size(); ++c) {
		const Report &rep = reps_out[c];
		const double compr =
		    rep.total_len > 0 ? static_cast<double>(rep.nodes) /
		                            static_cast<double>(rep.total_len)
		                      : 0.0;
		std::cout << "BENCH " << file << " config=" << cfgs[c].name
		          << " build_s=" << rep.build_s << " V=" << rep.nodes
		          << " E=" << rep.edges << " compr=" << compr
		          << " spread=" << rep.spread_mean << "/" << rep.spread_max
		          << " queries=" << tc.queries.size() << '\n';
	}
	std::cout << "PASS " << file << " (query results identical across "
	          << cfgs.size() << " band configs)\n";
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

	// Full band covers every column: the exact construction reference.
	// Default and narrow bands trade sharing for build speed.
	const std::vector<BandConfig> cfgs = {
	    {"full", std::numeric_limits<int>::max() / 2, 0.0},
	    {"default", DEFAULT_BAND_BASE, DEFAULT_BAND_SLOPE},
	    {"narrow", 2, 0.0},
	};

	const auto files = collect_files(argc, argv, first);
	if (files.empty()) {
		std::cout << "no testcase files found, nothing to do\n";
		return 0;
	}
	for (const auto &file : files)
		bench_testcase(file, cfgs, max_seqs, reps);

	std::cout << files.size() << " testcase file(s) benchmarked\n";
	return 0;
}
