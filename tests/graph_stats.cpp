/**
 * @file graph_stats.cpp
 * @brief Prints POA graph and BinaryLifter statistics for testcase files.
 *
 * For each testcase: dictionary/sequence-length summary, POA graph shape
 * (nodes, edges, heavy/light mix, compression/sharing), compressed-path
 * statistics, and per-query BinaryLifter size and construction time (the
 * per-query cost hlp_grep pays on top of the shared graph).
 *
 * Usage: graph_stats <testcase-file-or-dir>...
 */

#include <hlp_grep/hlp_grep.hpp>

#include "testcase.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace hlp_grep;

namespace {

std::string summary(const std::vector<std::size_t> &xs) {
	if (xs.empty())
		return "n/a";
	std::vector<std::size_t> s = xs;
	std::sort(s.begin(), s.end());
	return "min=" + std::to_string(s.front()) + " max=" + std::to_string(s.back())
	       + " median=" + std::to_string(s[s.size() / 2])
	       + " mean=" + std::to_string([&] {
		          std::size_t t = 0;
		          for (const std::size_t v : s)
			          t += v;
		          return t / s.size();
	          }());
}

void stats_testcase(const fs::path &file, const Testcase &tc) {
	std::cout << "=== " << file << " ===\n";

	std::vector<std::size_t> seq_lens;
	seq_lens.reserve(tc.dict.size());
	for (const auto &seq : tc.dict)
		seq_lens.push_back(seq.size());
	std::cout << "dict: " << tc.dict.size() << " sequences, "
	          << summary(seq_lens) << " (total bases "
	          << std::accumulate(seq_lens.begin(), seq_lens.end(),
		                         std::size_t{0}) << ")\n";

	const std::chrono::steady_clock::time_point start =
	    std::chrono::steady_clock::now();
	const POAGraph graph(tc.dict, *tc.cost);
	const std::chrono::steady_clock::time_point end =
	    std::chrono::steady_clock::now();

	// Unique edges and their heavy/light mix, straight from the paths.
	std::set<std::pair<POAGraph::node_id, POAGraph::node_id>> edges;
	std::size_t light_edges = 0, heavy_total = 0;
	std::vector<std::size_t> path_lens, light_per_path, heavy_per_path;
	for (std::size_t i = 0; i < graph.num_sequences(); ++i) {
		const auto &p = graph.path(i);
		path_lens.push_back(p.size());
		std::size_t heavy = 0, light = 0, dup = 0;
		for (std::size_t j = 1; j < p.size(); ++j) {
			const POAGraph::node_id u = p[j - 1], v = p[j];
			const bool fresh = edges.insert({u, v}).second;
			if (graph.edge_type(u, v) == POAGraph::EdgeType::HEAVY) {
				if (fresh)
					heavy_total++; // each heavy edge exists exactly once
				heavy++;
			} else {
				if (fresh)
					light_edges++;
				light++;
			}
			(void)dup;
		}
		heavy_per_path.push_back(heavy);
		light_per_path.push_back(light);
	}

	std::cout << "graph: " << graph.num_nodes() << " nodes, "
	          << (heavy_total + light_edges) << " unique edges ("
	          << heavy_total << " heavy, " << light_edges
	          << " light); built in "
	          << std::chrono::duration<double, std::milli>(end - start).count()
	          << " ms\n"
	          << "  path steps per sequence: " << summary(path_lens) << "\n"
	          << "  light steps per sequence: " << summary(light_per_path)
	          << "\n"
	          << "  heavy steps per sequence: " << summary(heavy_per_path)
	          << "\n"
	          << "  sharing: sum(path lengths)/nodes = "
	          << (graph.num_nodes()
			          ? static_cast<double>(std::accumulate(path_lens.begin(),
			                                                 path_lens.end(),
			                                                 std::size_t{0}))
			            / static_cast<double>(graph.num_nodes())
			          : 0.0)
	          << "\n"
	          << "  compression: total bases / nodes = "
	          << (graph.num_nodes()
			          ? static_cast<double>(std::accumulate(seq_lens.begin(),
			                                                 seq_lens.end(),
			                                                 std::size_t{0}))
			            / static_cast<double>(graph.num_nodes())
			          : 0.0)
	          << "\n";

	// Where the extra nodes live: seed path vs bubbles. The first stored
	// path seeds the graph; every later insertion event appends nodes of
	// its own, and every deletion/insertion creates branching points.
	{
		const auto &seed = graph.path(0);
		std::set<POAGraph::node_id> on_seed(seed.begin(), seed.end());
		std::size_t branching = 0, max_out = 0, chain_nodes = 0, sinks = 0;
		std::vector<std::size_t> outdeg(graph.num_nodes(), 0);
		for (const auto &[u, v] : edges) {
			(void)v;
			outdeg[u]++;
		}
		for (const std::size_t d : outdeg) {
			if (d > 1)
				branching++;
			if (d == 1)
				chain_nodes++;
			if (d == 0)
				sinks++;
			max_out = std::max(max_out, d);
		}
		const auto pct = [&graph](std::size_t c) {
			return graph.num_nodes()
			           ? 100.0 * c / graph.num_nodes()
			           : 0.0;
		};
		std::cout << "  seed path: " << seed.size()
		          << " nodes; nodes off the seed path: "
		          << graph.num_nodes() - on_seed.size() << "\n"
		          << "  branching nodes (out-degree > 1): " << branching
		          << " (" << pct(branching) << "%), single-out-edge nodes: "
		          << chain_nodes << " (" << pct(chain_nodes) << "%), sinks: "
		          << sinks << " (" << pct(sinks) << "%)"
		          << ", max out-degree: " << max_out << "\n";
	}

	// Node position spreads: the j_min/j_max drives the lifter's
	// precomputation windows [j_min - k, j_max + k + 1].
	const std::size_t START = std::numeric_limits<std::size_t>::max();
	std::vector<std::size_t> spreads;
	for (std::size_t u = 0; u < graph.num_nodes(); ++u) {
		const auto [lo, hi] = graph.pos_range(u);
		if (lo != START)
			spreads.push_back(hi - lo);
	}
	std::cout << "  node offset spread (pos_max - pos_min): "
	          << summary(spreads) << "\n";

	// Per-query BinaryLifter: rebuilt per query, its chain tables scale
	// with nodes * max_level * window width; time its construction.
	std::size_t max_level = 0;
	while ((std::size_t{1} << max_level) <= graph.num_nodes())
		++max_level;
	std::cout << "  lifter levels = " << max_level
	          << " (nodes * levels = "
	          << graph.num_nodes() * std::max<std::size_t>(max_level, 1) << ")\n";
	for (std::size_t qi = 0; qi < tc.queries.size() && qi < 3; ++qi) {
		const auto [k, query] = tc.queries[qi];
		const auto b = std::chrono::steady_clock::now();
		BinaryLifter lifter(graph, query, k);
		const auto e = std::chrono::steady_clock::now();
		(void)lifter;
		std::cout << "  query " << qi << " (len=" << query.size()
		          << ", k=" << k << "): BinaryLifter built in "
		          << std::chrono::duration<double, std::milli>(e - b).count()
		          << " ms\n";
	}
	std::cout << "\n";
}

} // namespace

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cerr << "usage: " << argv[0] << " <testcase-file-or-dir>...\n";
		return 1;
	}

	const auto files = collect_files(argc, argv);
	if (files.empty()) {
		std::cout << "no testcase files found\n";
		return 0;
	}
	for (const auto &file : files)
		stats_testcase(file, parse_testcase(file));
	return 0;
}
