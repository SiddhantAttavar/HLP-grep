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
#include <map>
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

void plot_weight_distribution(const char *label,
                             const std::vector<std::size_t> &weights,
                             std::size_t trunk) {
	std::cout << "  " << label << " (" << weights.size() << " values):\n";
	if (weights.empty()) {
		std::cout << "    n/a\n";
		return;
	}

	const std::size_t max_weight =
	    *std::max_element(weights.begin(), weights.end());
	const std::size_t bin_width =
	    std::max<std::size_t>(1, 1 + (max_weight - 1) / 20);
	const std::size_t bin_count = 1 + (max_weight - 1) / bin_width;
	std::vector<std::size_t> bins(bin_count, 0);
	for (const std::size_t weight : weights)
		++bins[(weight - 1) / bin_width];

	const std::size_t max_count =
	    *std::max_element(bins.begin(), bins.end());
	for (std::size_t i = 0; i < bins.size(); ++i) {
		const std::size_t lo = i * bin_width + 1;
		const std::size_t hi = std::min(max_weight, (i + 1) * bin_width);
		const std::size_t bar_width = bins[i] * 40 / max_count;
		std::cout << "    " << lo;
		if (hi != lo)
			std::cout << "-" << hi;
		std::cout << " | " << std::string(bar_width, '#') << " " << bins[i];
		if (trunk >= lo && trunk <= hi)
			std::cout << "  <-- trunk=" << trunk;
		std::cout << "\n";
	}
}

void stats_testcase(const fs::path &file, const Testcase &tc,
                    bool compact_nodes, bool plot_weight_dist) {
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
	const POAGraph graph(tc.dict, tc.cost, compact_nodes);
	const std::chrono::steady_clock::time_point end =
	    std::chrono::steady_clock::now();

	// Unique edges and their heavy/light mix, straight from the paths.
	std::map<std::pair<POAGraph::node_id, POAGraph::node_id>, std::size_t>
	    edge_weights;
	std::size_t light_edges = 0, heavy_total = 0;
	std::vector<std::size_t> path_lens, light_per_path, heavy_per_path;
	for (std::size_t i = 0; i < graph.num_sequences(); ++i) {
		const auto &p = graph.path(i);
		path_lens.push_back(p.size());
		std::size_t heavy = 0, light = 0;
		for (std::size_t j = 1; j < p.size(); ++j) {
			const POAGraph::node_id u = p[j - 1], v = p[j];
			auto [edge, fresh] = edge_weights.try_emplace({u, v}, 0);
			++edge->second;
			if (graph.edge_type(u, v) == POAGraph::EdgeType::HEAVY) {
				if (fresh)
					heavy_total++; // each heavy edge exists exactly once
				heavy++;
			} else {
				if (fresh)
					light_edges++;
				light++;
			}
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

	if (plot_weight_dist) {
		std::vector<std::size_t> all_edge_weights;
		all_edge_weights.reserve(edge_weights.size());
		std::vector<std::size_t> best_outgoing_weights(graph.num_nodes(), 0);
		for (const auto &[edge, weight] : edge_weights) {
			all_edge_weights.push_back(weight);
			best_outgoing_weights[edge.first] =
			    std::max(best_outgoing_weights[edge.first], weight);
		}
		best_outgoing_weights.erase(
		    std::remove(best_outgoing_weights.begin(), best_outgoing_weights.end(), 0),
		    best_outgoing_weights.end());
		const std::size_t trunk = graph.num_sequences() / 2;
		const std::size_t eligible = static_cast<std::size_t>(std::count_if(
		    best_outgoing_weights.begin(), best_outgoing_weights.end(),
		    [trunk](std::size_t weight) { return weight >= trunk; }));
		std::cout << "  weight distribution: trunk=" << trunk << ", " << eligible
		          << "/" << best_outgoing_weights.size()
		          << " node-best edges meet the cutoff\n";
		plot_weight_distribution("all edge traversal weights", all_edge_weights,
		                         trunk);
		plot_weight_distribution("best outgoing edge weight per node",
		                         best_outgoing_weights, trunk);
		std::cout << "\n";
	}

	// Where the extra nodes live: seed path vs bubbles. The first stored
	// path seeds the graph; every later insertion event appends nodes of
	// its own, and every deletion/insertion creates branching points.
	{
		const auto &seed = graph.path(0);
		std::set<POAGraph::node_id> on_seed(seed.begin(), seed.end());
		std::size_t branching = 0, max_out = 0, chain_nodes = 0, sinks = 0;
		std::vector<std::size_t> outdeg(graph.num_nodes(), 0);
		for (const auto &[edge, weight] : edge_weights) {
			(void)weight;
			outdeg[edge.first]++;
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

	// Window outliers: nodes whose pos_range sentinel/garbage values would
	// blow up the BinaryLifter window (its width would otherwise be
	// clamp-wise proportional to the spread).
	std::size_t wide_nodes = 0, start_nodes = 0;
	std::size_t widest = 0;
	POAGraph::node_id wide_id = 0;
	std::pair<std::size_t, std::size_t> wide_range{0, 0};
	for (std::size_t u = 0; u < graph.num_nodes(); ++u) {
		const auto [lo, hi] = graph.pos_range(u);
		if (lo == START) {
			start_nodes++;
			continue;
		}
		const std::size_t width = hi - lo;
		if (width > widest) {
			widest = width;
			wide_id = u;
			wide_range = {lo, hi};
		}
		if (width > 1000)
			wide_nodes++;
	}
	std::cout << "  window outliers: START nodes=" << start_nodes
	          << " nodes wider than 1000: " << wide_nodes
	          << " widest node " << wide_id << " := ["
	          << wide_range.first << ", " << wide_range.second
	          << ") label len "
	          << graph.node(wide_id).seq.size() << "\n\n";

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
		BinaryLifter lifter(graph);
		const auto e = std::chrono::steady_clock::now();
		(void)lifter;
		std::cout << "  query " << qi << " (len=" << query.size()
		          << ", k=" << k << "): BinaryLifter up-table built in "
		          << std::chrono::duration<double, std::milli>(e - b).count()
		          << " ms\n";
	}
	std::cout << "\n";
}

} // namespace

int main(int argc, char **argv) {
	bool compact_nodes = true;
	bool plot_weight_dist = false;
	std::vector<const char *> args;
	for (int i = 1; i < argc; ++i) {
		if (std::string(argv[i]) == "--no-compact") {
			compact_nodes = false;
		} else if (std::string(argv[i]) == "--plot-weight-dist") {
			plot_weight_dist = true;
		} else {
			args.push_back(argv[i]);
		}
	}
	if (args.empty()) {
		std::cerr << "usage: " << argv[0]
		          << " [--no-compact] [--plot-weight-dist]"
		             " <testcase-file-or-dir>...\n";
		return 1;
	}

	// collect_files() indexes argv terms from 1 (argv[0] is the program
	// name), so keep a dummy first element.
	auto arg_ptrs = [&] {
		std::vector<char *> p{const_cast<char *>(argv[0])};
		for (const char *a : args)
			p.push_back(const_cast<char *>(a));
		return p;
	}();
	const auto files = collect_files(static_cast<int>(arg_ptrs.size()),
	                                 arg_ptrs.data());	if (files.empty()) {
		std::cout << "no testcase files found\n";
		return 0;
	}
	std::cout << "run compaction: " << (compact_nodes ? "on" : "off") << "\n";
	for (const auto &file : files) {
		stats_testcase(file, parse_testcase(file), compact_nodes,
		               plot_weight_dist);
	}
	return 0;
}
