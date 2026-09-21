/**
 * @file lifter_stats.cpp
 * @brief Statistics on how much of the eagerly precomputed BinaryLifter
 *        chain tables is actually used, aggregated over the queries of
 *        a testcase.
 *
 * Answers the lazy-build question: every (u, t) table entry the
 * constructor builds eagerly would need, when computed lazily, its two
 * half-chain children (u, t-1) and (u ^ 2^(t-1), t-1) first. Stats:
 * allocated entries and cells, used entries during the queries, used
 * union across queries, per-level histograms, and the dependency
 * closure of the used set (tables a lazy implementation would still
 * have to compute even if it only realizes the used tables).
 *
 * The scoring walk mirrors Solver::query (initial row + jump/step +
 * early-abort) so usage matches real query behavior exactly; it must
 * stay in sync with hlp_grep.hpp.
 *
 * Usage: lifter_stats <testcase-file-or-dir>...
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

using UT = std::pair<std::size_t, std::size_t>; // (node, level)

/// Initial DP row for a path's start node; mirrors Solver::initial_row.
std::vector<int> initial_row(const POAGraph &graph, const CostModel &cost,
                             POAGraph::node_id start,
                             const std::string &query, int k) {
	const auto [jmin, jmax] = graph.pos_range(start);
	const long m = static_cast<long>(query.size());
	const long mk = static_cast<long>(k);
	const long lo = std::max(0L, static_cast<long>(jmin) - mk);
	const long hi = std::min(m + 1, static_cast<long>(jmax) + mk + 2);

	std::vector<int> row;
	row.reserve(static_cast<std::size_t>(hi - lo));
	long best_alt = DistMatrix::INF;
	long sum_ins = 0;
	for (long i = 0; i < hi; ++i) {
		if (i > 0) {
			sum_ins += cost.ins(query[i - 1]);
			best_alt = std::min(
			    best_alt,
			    static_cast<long>(cost.consume(graph.base(start),
			                                   query[i - 1])) -
			        static_cast<long>(cost.ins(query[i - 1])));
		}
		long cur = static_cast<long>(cost.del(graph.base(start))) + sum_ins;
		if (i > 0)
			cur = std::min(cur, best_alt + sum_ins);
		if (cur > DistMatrix::INF)
			cur = DistMatrix::INF;
		if (i >= lo)
			row.push_back(static_cast<int>(cur));
	}
	return row;
}

void stats_testcase(const fs::path &file, const Testcase &tc) {
	std::cout << "=== " << file << " ===\n";

	std::vector<std::size_t> seq_lens;
	for (const auto &seq : tc.dict)
		seq_lens.push_back(seq.size());
	std::cout << "dict: " << tc.dict.size() << " sequences (total bases "
	          << std::accumulate(seq_lens.begin(), seq_lens.end(),
		                         std::size_t{0}) << ", "
	          << tc.queries.size() << " queries)\n";
	if (tc.queries.empty()) {
		std::cout << "no queries\n\n";
		return;
	}

	const POAGraph graph(tc.dict, *tc.cost);
	const std::size_t n = graph.num_nodes();
	std::size_t max_level = 0;
	while ((std::size_t{1} << max_level) <= n)
		++max_level;
	if (n == 0) {
		std::cout << "empty graph\n\n";
		return;
	}

	// Heavy-chain layout: all tables (u, t) with a sub-chain are exactly
	// the (u, off(u) + 2^t <= chain_len) entries along each heavy chain
	// (mirrors BinaryLifter's up table, built bottom-up identically
	// without a topological pass).
	const auto &[k0, q0] = tc.queries[0];
	const long m0 = static_cast<long>(q0.size());
	std::vector<std::size_t> chain_id(n, SIZE_MAX), chain_off(n, SIZE_MAX);
	std::vector<std::size_t> chain_len; // number of nodes per chain
	{
		std::vector<std::size_t> next(n, SIZE_MAX);
		for (std::size_t u = 0; u < n; ++u)
			if (graph.node(u).heavy_neighbour.has_value())
				next[u] = *graph.node(u).heavy_neighbour;
		std::vector<char> is_start(n, true);
		std::set<std::size_t> targets;
		for (std::size_t u = 0; u < n; ++u)
			if (next[u] != SIZE_MAX)
				targets.insert(next[u]);
		for (const std::size_t t : targets)
			is_start[t] = false;
		for (std::size_t u = 0; u < n; ++u) {
			if (!is_start[u])
				continue;
			const std::size_t cid = chain_len.size();
			chain_len.push_back(0);
			std::size_t v = u;
			while (true) {
				chain_id[v] = cid;
				chain_off[v] = chain_len[cid]++;
				if (next[v] == SIZE_MAX)
					break;
				v = next[v];
			}
		}
	}

	// Heavy-chain steps available from u: (chain length - 1 - offset) is
	// the max number of heavy edges walkable from u.
	auto steps_from = [&](std::size_t u) -> std::size_t {
		if (chain_id[u] == SIZE_MAX || chain_off[u] == SIZE_MAX)
			return 0;
		return chain_len[chain_id[u]] - 1 - chain_off[u];
	};

	// Dependency walk: heavy_neighbour chain, used to find the node
	// 2^s steps ahead of u (mirrors up[u][t]).
	auto ahead = [&graph](std::size_t u, std::size_t steps) {
		POAGraph::node_id v = u;
		for (std::size_t i = 0; i < steps; ++i) {
			const auto &node = graph.node(v);
			if (!node.heavy_neighbour.has_value())
				return POAGraph::node_id(
				    std::numeric_limits<std::size_t>::max());
			v = *node.heavy_neighbour;
		}
		return v;
	};

	// Window width of a node under the first (probe) query; mirrors
	// BinaryLifter::window with (k0, m0).
	auto window_width = [&](std::size_t u) -> std::size_t {
		const auto [jmin, jmax] = graph.pos_range(u);
		const long lo = std::max(0L, static_cast<long>(jmin) -
		                              static_cast<long>(k0));
		const long hi = std::min(m0 + 1, static_cast<long>(jmax) +
		                            static_cast<long>(k0) + 2);
		return static_cast<std::size_t>(std::max(hi, lo + 1) - lo);
	};

	// Allocated (u, t) entries are the ones BinaryLifter builds eagerly:
	// node u must have 2^t heavy steps available.
	std::size_t allocated = 0, allocated_cells = 0;
	std::vector<std::size_t> alloc_by_level(max_level, 0);
	for (std::size_t u = 0; u < n; ++u)
		for (std::size_t t = 0; t < max_level; ++t) {
			const std::size_t steps = std::size_t{1} << t;
			if (steps_from(u) < steps)
				continue;
			alloc_by_level[t]++;
			allocated++;
			allocated_cells += window_width(u) * window_width(ahead(u, steps));
		}

	// Replay Solver::query for every query on a per-query lifter,
	// tracking which (u, bit) tables each jump applies (jump applies
	// the bits of st.length in order, walking the chain).
	std::size_t jumps = 0, applied = 0;
	std::vector<std::size_t> used_by_level(max_level, 0);
	std::set<UT> used_union;
	for (const auto &[k, query] : tc.queries) {
		BinaryLifter lifter(graph, query, k);
		const long m = static_cast<long>(query.size());
		for (std::size_t i = 0; i < tc.dict.size(); ++i) {
			if (std::abs(static_cast<long>(tc.dict[i].size()) - m) > k)
				continue;
			const auto &cp = graph.compressed_path(i);
			std::vector<int> row = initial_row(graph, *tc.cost, cp.start,
			                                   query, k);
			POAGraph::node_id cur = cp.start;
			for (const auto &st : cp.steps) {
				if (st.type == POAGraph::EdgeType::HEAVY) {
					for (std::size_t len = st.length, bit = 0;
					     len > 0; len >>= 1, ++bit) {
						if ((len & 1) == 0)
							continue;
						applied++;
						if (used_union.emplace(cur, bit).second)
							used_by_level[bit]++;
					}
					cur = lifter.jump(cur, st.length, row);
					jumps++;
				} else {
					row = lifter.step(cur, st.next, row);
					cur = st.next;
				}
				if (*std::min_element(row.begin(), row.end()) > k)
					break;
			}
		}
	}
	std::cout << "allocated tables: " << allocated << " entries, "
	          << allocated_cells << " cells of possible "
	          << (n * max_level) << "\n";
	std::cout << "used tables (union over " << tc.queries.size()
	          << " queries): " << used_union.size() << " ("
	          << (100.0 * used_union.size() / std::max(allocated, 1UL))
	          << "% of allocated); " << jumps << " heavy jumps, "
	          << applied << " table applications\n";

	// Lazy dependency closure: realizing (u, t) recursively needs
	// (u, t - 1) and (u ^ 2^(t-1), t - 1). Count the distinct (u', t')
	// tables even a lazy implementation would have to compute.
	std::set<UT> closure;
	std::vector<UT> work(used_union.begin(), used_union.end());
	while (!work.empty()) {
		const auto [u, t] = work.back();
		work.pop_back();
		if (!closure.emplace(u, t).second)
			continue;
		if (t == 0)
			continue;
		const std::size_t half = 1UL << (t - 1);
		work.emplace_back(u, t - 1);
		const POAGraph::node_id v = ahead(u, half);
		if (v < n)
			work.emplace_back(v, t - 1);
	}
	std::cout << "lazy closure of used set: " << closure.size()
	          << " entries ("
	          << (100.0 * closure.size() / std::max(allocated, 1UL))
	          << "% of allocated)\n";
	for (std::size_t t = 0; t < max_level; ++t) {
		if (alloc_by_level[t] == 0 && used_by_level[t] == 0)
			continue;
		std::cout << "  level " << t << " (2^" << t << " steps): alloc "
		          << alloc_by_level[t] << ", used " << used_by_level[t]
		          << " ("
		          << (100.0 * used_by_level[t] /
		              std::max(alloc_by_level[t], 1UL))
		          << "%)\n";
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
