/**
 * @file binary_lifter.hpp
 * @brief Binary lifting over the heavy chains of a POAGraph.
 *
 * Precomputes, for every graph node and every power of two, the node reached
 * by walking that many heavy edges. Alongside the node table, per-level
 * chain distances ed(u, 2^t, a, b) are precomputed as min-plus products of
 * level (t - 1) blocks. The chain tables depend on the
 * query string and threshold k, so both are supplied at construction time
 * and fixed for the lifetime of the object. Jump queries advance a
 * caller-provided DP row across the jumped chain with min-plus products.
 */
#pragma once
#include <hlp_grep/dist_matrix.hpp>
#include <hlp_grep/poa_graph.hpp>

#include <omp.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace hlp_grep {

/**
 * @brief Binary lifting table over the heavy chains of a POAGraph.
 *
 * Entry (u, t) of the node table is the node 2^t heavy steps from u, or
 * empty once fewer than 2^t heavy steps remain from u. Nodes with no
 * outgoing heavy edge (sinks) start the empty entries. No sentinel node id
 * is used, so the table cannot be confused with POAGraph's internal
 * markers.
 *
 * The lifter is built in two stages. The constructor builds only the
 * query-independent node table; decompose() expands the dictionary's
 * compressed heavy runs into the exact blocks used by jump() and records
 * which (node, level) chain blocks queries can need. build(query, k)
 * materializes only those blocks for one (query, k) pair.
 */
class BinaryLifter {
public:
	/// Id of a graph node (mirrors POAGraph::node_id).
	using node_id = POAGraph::node_id;

	/**
	 * @brief Constructs the lifter's query-independent top half.
	 *
	 * The node table has one row per graph node and one column per level
	 * 0 .. max_level - 1, where max_level = floor(log2(num_nodes)) + 1 is
	 * an upper bound on the useful power-of-two step counts. Built level
	 * by level from the per-node heavy outgoing edge pointers:
	 * up[u][t] = up[up[u][t - 1]][t - 1] whenever up[u][t - 1] is set, so
	 * no topological order is needed.
	 *
	 * The chain tables are NOT built here: decompose() records which blocks
	 * queries touch and build(query, k) materializes exactly those. jump
	 * and step must not be called before build().
	 *
	 * @param graph POAGraph whose heavy chains are indexed; must outlive
	 *              this object and must not be modified while it is alive.
	 * @param num_threads Worker threads for the chain table build; 0 means
	 *              use up to half of the hardware concurrency.
	 */
	explicit BinaryLifter(const POAGraph &graph,
	                      std::size_t num_threads = 0)
	    : graph(graph), cost_model(graph.cost_model()),
	      num_threads(num_threads) {
		const std::size_t n = graph.num_nodes();
		max_level = 0;
		while ((std::size_t{1} << max_level) <= n)
			++max_level; // floor(log2(n)) + 1; 0 for n == 0
		up.assign(n, std::vector<std::optional<node_id>>(max_level));
		for (std::size_t u = 0; u < n; ++u)
			up[u][0] = graph.node(u).heavy_neighbour;
		for (std::size_t t = 1; t < max_level; ++t)
			for (std::size_t u = 0; u < n; ++u)
				if (up[u][t - 1])
					up[u][t] = up[*up[u][t - 1]][t - 1];
	}

	/**
	 * @brief Decomposes heavy runs in the paths and records which (node,
	 *        level) chain blocks queries touch.
	 *
	 * Each HEAVY step is replaced in-place by the greedy blocks jump()
	 * applies: at node v with rem steps left, the level is the smaller of
	 * the lowest set bit of heavy_length(v) and the highest set bit of rem.
	 * This depends only on the graph/path structure, so the decomposition
	 * and touched set are query-independent. Scoring can abort rows whose
	 * minimum exceeds k; decompose() cannot evaluate that, so it records
	 * the abort-free superset for the dictionary. The
	 * dependency closure is included: block (u, t) needs (u, t - 1) and
	 * (*up[u][t - 1], t - 1), recursively down to level 0, so the
	 * per-level build lists contain exactly the blocks the products
	 * read. decompose() must run before build(); it can be called again
	 * because already decomposed power-of-two steps are fixed points. It
	 * also pre-sizes up_mat: node
	 * u's matrix vector holds levels 0 .. its largest touched level, so
	 * build() fills slots in place instead of re-assigning the whole
	 * structure per query.
	 *
	 * @param paths Compressed paths, one per dictionary sequence; modified
	 *              in place to split heavy runs into jump blocks.
	 */
	void decompose(std::vector<POAGraph::CompressedPath> &paths) {
		const std::size_t n = graph.num_nodes();
		touched.assign(n, std::vector<char>(max_level, 0));
		for (auto &cp : paths) {
			node_id cur = cp.start;
			std::vector<POAGraph::CompressedEdge> decomposed;
			decomposed.reserve(cp.steps.size());
			for (const auto &st : cp.steps) {
				if (st.type == POAGraph::EdgeType::HEAVY) {
					std::size_t rem = st.length;
					while (rem > 0) {
						const std::size_t level =
						    std::min(
						        static_cast<unsigned long long>(
						            __builtin_ctzll(graph.node(
						                cur).heavy_length)),
						        static_cast<unsigned long long>(
						            63 - __builtin_clzll(rem)));
						const std::size_t block = std::size_t{1} << level;
						POAGraph::CompressedEdge block_step;
						block_step.type = POAGraph::EdgeType::HEAVY;
						block_step.length = static_cast<int>(block);
						decomposed.push_back(block_step);
						touched[cur][level] = true;
						cur = *up[cur][level];
						rem -= block;
					}
				} else {
					decomposed.push_back(st);
					cur = st.next;
				}
			}
			cp.steps = std::move(decomposed);
		}
		// Dependency closure: (u, t) needs (u, t - 1) and its chain
		// midpoint's (t - 1) block. One sweep per level from the top
		// catches every transitive requirement.
		for (std::size_t t = max_level; t-- > 1;) {
			for (std::size_t u = 0; u < n; ++u)
				if (touched[u][t]) {
					touched[u][t - 1] = true;
					touched[*up[u][t - 1]][t - 1] = true;
				}
		}
		// Per-level build lists (ascending node order).
		build_nodes.assign(max_level, {});
		for (std::size_t t = 0; t < max_level; ++t) {
			for (std::size_t u = 0; u < n; ++u)
				if (touched[u][t])
					build_nodes[t].push_back(u);
		}
		// Pre-size up_mat: since (u, t) touched drags (u, t - 1) into
		// the closure, each node's touched levels form a prefix of
		// 0..its largest touched level; one vector per node sized to
		// the top of that prefix. build() then fills the slots in
		// place, and build() may not run without decompose().
		up_mat.assign(n, {});
		for (std::size_t u = 0; u < n; ++u) {
			std::size_t last = 0;
			for (std::size_t t = max_level; t-- > 0;)
				if (touched[u][t]) {
					last = t + 1;
					break;
				}
			up_mat[u].resize(last, DistMatrix(0, 0));
		}
	}

	/**
	 * @brief Materializes the chain tables for one (query, k) pair.
	 *
	 * Only the blocks marked by decompose() are built — level 0 holds the
	 * per-edge blocks of the touched heavy outgoing edges (the edge
	 * transition plus the destination node's label), and each higher
	 * level is the min-plus product of two adjacent level (t - 1) blocks
	 * sharing a chain midpoint. Untouched (u, t) slots stay empty.
	 * Must be called after decompose() (which pre-sizes up_mat) and before
	 * jump()/step()/window() reads; rebuilds every marked slot
	 * unconditionally on every call.
	 *
	 * @param query Query string the chain distance tables index into.
	 * @param k     Edit distance threshold widening the windows.
	 */
	void build(const std::string &query, int k) {
		// The per-node blocks are one task each with uniform cost; past
		// ~8 threads barrier and spawn overhead dominates the level loop.
		const int threads = static_cast<int>(
		    num_threads == 0 ? std::max(omp_get_max_threads() / 2, 1)
		                     : num_threads);
		this->query = query;
		this->k = k;
#pragma omp parallel num_threads(threads)
		{
			for (std::size_t t = 0; t < max_level; ++t) {
#pragma omp for schedule(static)
				for (std::size_t idx = 0;
				     idx < build_nodes[t].size(); ++idx) {
					const std::size_t u = build_nodes[t][idx];
					if (t == 0)
						up_mat[u][0] =
						    edge_matrix(u, *up[u][0]);
					else {
						const std::size_t v =
						    *up[u][t - 1];
						up_mat[u][t] =
						    DistMatrix::min_plus_product(
						        up_mat[u][t - 1],
						        up_mat[v][t - 1]);
					}
				}
				// The for's implicit barrier separates
				// consecutive levels (a level's blocks may be
				// read by the next one).
			}
		}
	}

	/**
	 * @brief Jumps l heavy steps from u (u^l; l == 0 is
	 *        the identity) and advances the caller's DP row across the
	 *        jumped chain.
	 *
	 * @p cost and @p pos_range hold the start node's DP row, possibly
	 * clipped to a contiguous subrange of its after-label window. They are
	 * advanced with one DistMatrix::min_plus_apply over the precomputed power-of-two
	 * block. decompose() has already split path runs so each call's length
	 * is a power of two and the corresponding level exists at @p u.
	 * Afterwards @p cost holds one entry per allowed position of the reached
	 * node, then clipped to the contiguous span between the first and last
	 * positions whose costs are at most k. In effect, after the call
	 *   cost(b) = min_a cost_in(a) + ed(u, l, a, b)
	 * for the chain [u, u^1, ..., u^l] (each heavy edge crossing its
	 * destination node's full label). Cells no transition realizes carry
	 * DistMatrix::INF, which dominates any real cost while keeping every
	 * entry within the saturating product's overflow-free range.
	 *
	 * @param u    Node to start from.
	 * @param l    Number of heavy steps to walk.
	 * @param cost DP row of the start node; replaced by the clipped row
	 *             of the reached node.
	 * @param pos_range Position range represented by @p cost; replaced by
	 *                  the clipped range of the reached node.
	 * @return The node l heavy steps from u.
	 * @pre u is a valid node id, the row matches a subrange of u's
	 *      window, l is a power-of-two block recorded by decompose(),
	 *      and the block was materialized by build().
	 */
	node_id jump(node_id u, int l, std::vector<int> &cost,
	             std::pair<std::size_t, std::size_t> &pos_range) const {
		assert(u < up.size());
		assert_row(u, cost, pos_range);
		assert(l > 0);
		const std::size_t block = static_cast<std::size_t>(l);
		assert((block & (block - 1)) == 0);
		const std::size_t level = 63 - __builtin_clzll(block);
		assert(level < max_level && up[u][level].has_value());
		assert(level < up_mat[u].size() &&
		       up_mat[u][level].num_cols() > 0);
		const DistMatrix &mat = up_mat[u][level];
		std::vector<int> row =
		    mat.min_plus_apply(cost, pos_range);
		const node_id v = *up[u][level];
		pos_range = mat.col_range();
		assert(pos_range == window(v));
		clip_row(row, pos_range);
		cost = std::move(row);
		return v;
	}

	/**
	 * @brief Steps a DP row across the single edge (u, v): applies the
	 *        edge transition plus v's whole label directly to the input
	 *        row and returns the advanced row.
	 *
	 * @p cost and @p pos_range hold the allowed query positions at u,
	 * possibly clipped to a contiguous subrange; on return they hold v's row
	 * clipped to a subrange as well. The row is computed by a semiglobal ED DP that
	 * superposes all source
	 * positions: row 0 seeds the source costs within the entry band
	 * window(v, 0) — sources below it die, since every continuation
	 * leaves all bands — chained left to right with insertions, and
	 * sweep_label() then advances each of v's label characters over its
	 * own band window(v, j) with the recurrence
	 *   f(i, b) = min( f(i - 1, b) + del(label[i - 1]),
	 *                  f(i, b - 1) + ins(query[b - 1]),
	 *                  f(i - 1, b - 1) + consume(label[i - 1],
	 *                                             query[b - 1]) )
	 * After |label| sweeps the row holds, at each b,
	 *   out(b) = min_a cost(a) + ed(u, v, a, b)
	 * — the exact min over all entry positions of the same transition
	 * edge_matrix() materializes, without building the matrix.
	 * Only the reachable span is swept: retained sources lie in
	 * [s_lo, s_hi) with costs within k, and reaching (j, b) from a
	 * costs at least b - a - j (loose ins = 1 form), so bands are
	 * clamped to [s_lo, s_hi + j + (k - cmin)) and the output to
	 * [s_lo, s_hi + |label| + (k - cmin)). Cells outside can only
	 * exceed k and are left INF for clip_row to drop. Single-character
	 * labels use a one-pass fast path over the same span instead of
	 * the generic sweep.
	 *
	 * @param u    Source node of the edge.
	 * @param v    Destination node of the edge.
	 * @param cost DP row of u; replaced by the clipped DP row at v.
	 * @param pos_range Position range represented by @p cost; replaced by
	 *                  the clipped range at v.
	 * @pre The row is clipped: non-empty (callers stop walking a path as
	 *      soon as it is empty) and holding an entry within k. It matches
	 *      a subrange of u's window, and there is an edge (u, v) in the
	 *      graph (not verified here).
	 */
	void step(node_id u, node_id v, std::vector<int> &cost,
	          std::pair<std::size_t, std::size_t> &pos_range) const {
		const std::string &label = graph.seq(v);
		const std::size_t label_len = label.size();
		const std::size_t m = query.size();
		if (label_len == 1) {
			auto [s_lo, s_hi] = pos_range;
			if (s_hi <= m) {
				s_hi++;
				cost.push_back(DistMatrix::INF);
			}
			const int del_c = cost_model.del();
			const char c = label[0];
			int cur = cost[0] + del_c;
			int diag = cost[0], left = cur;
			cost[0] = cur;
			const bool first = cur > k;
			for (long b = s_lo + 1; b < s_hi; ++b) {
				const int up = cost[b - s_lo];
				cur = std::min(up + del_c, left + cost_model.ins());
				cur = std::min(cur,
				   diag + cost_model.consume(c, query[b - 1]));
				cost[b - s_lo - first] = cur;
				left = cur;
				diag = up;
			}
			if (first) {
				s_lo++;
				cost.pop_back();
			}
			pos_range = {s_lo, s_hi};
			clip_row(cost, pos_range);
			return;
		}

		const auto [s_lo, s_hi] = pos_range;

		const long hL = std::min(m + 1, s_hi + label_len);
		cost.resize(static_cast<std::size_t>(hL - s_lo),
			  DistMatrix::INF);
		long li = s_lo;
		for (long i = 1; i <= label_len; ++i) {
			const long hi_i = std::min(hL, static_cast<long>(s_hi + i));
			const char c = label[i - 1];
			const int del_c = cost_model.del();
			int cur = cost[li - s_lo] + del_c;
			int diag = cost[li - s_lo], left = cur;
			cost[li - s_lo] = cur;
			for (long b = li + 1; b < hi_i; ++b) {
				const int up = cost[b - s_lo]; // f[i - 1][b]
				cur = std::min(up + del_c, left + cost_model.ins());
				cur = std::min(cur, diag + cost_model.consume(c, query[b - 1]));
				cost[b - s_lo] = cur;
				left = cur;
				diag = up;
			}
			while (li < hL && cost[li - s_lo] > k) {
				li++;
			}
			if (li >= hL) {
				cost.clear();
				pos_range = {hL, hL};
				return;
			}
		}
		std::move(cost.begin() + (li - s_lo), cost.end(), cost.begin());
		cost.resize(hL - li);
		pos_range = {li, hL};
		clip_row(cost, pos_range);
	}

	/** Half-open query-position window represented by a node's DP row. */
	std::pair<std::size_t, std::size_t> window(node_id u) const {
		const auto [lo, hi] = window(u, graph.seq(u).size());
		return {lo, std::max(hi, lo + 1)};
	}

	/** Removes leading and trailing row cells above the current threshold.
	 * Remaining query positions form a contiguous range; later transitions
	 * only advance query positions, so discarded endpoints cannot contribute
	 * to a future result within k.
	 */
	void clip_row(std::vector<int> &row,
	              std::pair<std::size_t, std::size_t> &pos_range) const {
		assert(row.size() == pos_range.second - pos_range.first);
		const std::size_t old_lo = pos_range.first;
		std::size_t first = 0;
		std::size_t last = row.size();
		while (first < last && row[first] > k)
			++first;
		while (last > first && row[last - 1] > k)
			--last;
		if (first == last) {
			row.clear();
			pos_range.first = pos_range.second;
			return;
		}
		if (first > 0)
			std::move(row.begin() + first, row.begin() + last, row.begin());
		row.resize(last - first);
		pos_range = {old_lo + first, old_lo + last};
	}

private:
	/** Debug-only check that a row matches a subrange of a node window. */
	void assert_row(
	    node_id u, const std::vector<int> &row,
	    const std::pair<std::size_t, std::size_t> &pos_range) const {
		const auto full_range = window(u);
		assert(pos_range.first >= full_range.first &&
		       pos_range.second <= full_range.second &&
		       pos_range.first <= pos_range.second &&
		       row.size() == pos_range.second - pos_range.first);
	}

	/**
	 * @brief The DistMatrix of the heavy edge (u, v): the transition
	 *        across the edge followed by v's whole label, i.e.
	 *        ed(u, v, a, b) between the edge landing on v and the query
	 *        interval [a, b].
	 *
	 * Used only at construction, to build the level-0 chain tables;
	 * step() applies the same transition directly without the matrix.
	 *
	 * Rows are the allowed query positions at u, columns those at v, i.e.
	 * the after-label windows of the two nodes (see window()). Entry
	 * (a, b) is the min cost to go from DP state (u, a) to (v, b) — the
	 * source position a seeds a semiglobal DP that aligns v's label
	 * against query[a..b) with insertions, deletions, and consumes,
	 * sweeping each label character over its own band window(v, j) so
	 * the construction costs O(rows * |label| * (spread + 2k)) instead
	 * of quadratic-in-|label| work. For a one-character label this
	 * reduces to the closed single-edge form; in general the matrix
	 * equals the min-plus product of the label's one-character
	 * transitions restricted to the band rectangle, so it satisfies the
	 * ED-DAG layer precondition of DistMatrix::min_plus_product, with
	 * INF cells (b < a, seeds below the entry band, band-infeasible
	 * states) being unreachable constants of that DAG.
	 */
	DistMatrix edge_matrix(node_id u, node_id v) const {
		const auto [lo_u, hi_u] = window(u);
		const auto [lo_v, hi_v] = window(v);
		const std::string &label = graph.seq(v);
		DistMatrix mat(static_cast<std::size_t>(hi_u - lo_u),
		               static_cast<std::size_t>(hi_v - lo_v), DistMatrix::INF,
		               {lo_u, hi_u}, {lo_v, hi_v});
		const auto [lL, hL] = window(v, label.size());
		if (lL >= hL) {
			mat.cache_col_finite();
			return mat; // final band empty: every crossing exceeds k
		}
		const auto [l0, h0] = window(v, 0);
		const std::size_t label_len = label.size();
		const std::size_t m = query.size();
		const std::size_t end = std::min(hi_u, hL);
		if (label_len == 1) {
			const char c = label[0];
			for (std::size_t a = std::max(lo_u, l0); a < end; ++a) {
				if (a >= lo_v)
					mat(a - lo_u, a - lo_v) = cost_model.del();
				int cur = 0, sub = cost_model.mismatch;
				std::size_t b_max = std::min(hi_v, a + k + 2);
				for (std::size_t b = a + 1; b < b_max; ++b) {
					if (query[b - 1] == c)
						sub = cost_model.match;
					if (b >= lo_v)
						mat(a - lo_u, b - lo_v) = cur + sub;
					cur += cost_model.ins();
				}
			}
			mat.cache_col_finite();
			return mat;
		}
		// DP row over [l0, hL); band j of the label occupies
		// window(v, j), and row L's band is exactly [lo_v, hi_v).
		std::vector<long> f(hL - l0);
		for (long a = std::max(lo_u, l0); a < end; ++a) {
			std::fill(f.begin(), f.end(), DistMatrix::INF);
			// Row 0: seed at a, chain insertions within the entry band.
			f[a - l0] = 0;
			for (std::size_t b = a + 1; b < h0; ++b)
				f[b - l0] = f[b - l0 - 1] + cost_model.ins();
			for (long i = 1; i <= label_len; ++i) {
				const long li = std::max(a, static_cast<long>(l0) + i - k);
				const long hi_i = std::min(hL, h0 + i);
				if (li >= hi_i) {
					std::fill(f.begin(), f.end(), DistMatrix::INF);
					break;
				}
				const char c = label[i - 1];
				const int del_c = cost_model.del();
				long cur = f[li - l0] + del_c;
				if (li > static_cast<long>(l0)) {
					cur = std::min(cur,
						f[li - l0 - 1] + cost_model.consume(c, query[li - 1]));
				}
				long diag = f[li - l0], left = cur;
				f[li - l0] = cur;
				for (long b = li + 1; b < hi_i; ++b) {
					const long up = f[b - l0]; // f[i - 1][b]
					cur = std::min(up + del_c, left + cost_model.ins());
					cur = std::min(cur, diag + cost_model.consume(c, query[b - 1]));
					f[b - l0] = cur;
					left = cur;
					diag = up;
				}
			}
			for (long b = std::max(a, static_cast<long>(lo_v)); b < hL; ++b)
				mat(a - lo_u, b - lo_v) = static_cast<int>(f[b - l0]);
		}
		mat.cache_col_finite();
		return mat;
	}

	/**
	 * @brief Half-open [lo, hi) band of query positions feasible at the
	 *        state @p offset label characters into u's label (offset 0 =
	 *        before the label, |label| = after it): the node occurs at
	 *        path offsets j_min + offset .. j_max + offset (bases
	 *        consumed up to that state), and any alignment of cost <= k
	 *        satisfies |b - consumed| <= k, so the band is
	 *        [j_min + offset - k, j_max + offset + k + 1) clamped to
	 *        [0, |query|]. Width spread + 2k + 1, independent of the
	 *        label length — sweeping each label character over its own
	 *        band keeps the DP linear in |label| times the window width.
	 *        The band may be empty (label running past |query|); the
	 *        bands of consecutive offsets overlap so that transitions
	 *        between them stay covered.
	 */
	std::pair<std::size_t, std::size_t> window(node_id u, std::size_t offset) const {
		const auto [jmin, jmax] = graph.pos_range(u);
		const std::size_t mk = static_cast<std::size_t>(k);
		const std::size_t m = query.size();
		const std::size_t lo = jmin + offset >= mk ? jmin + offset - mk : 0;
		const std::size_t hi = std::min(m + 1, jmax + offset + mk + 1);
		return {lo, hi};
	}

	/// Graph the lifter was built from; its heavy edges must not be
	/// modified while this object is alive.
	const POAGraph &graph;
	/// Query string the chain tables were built for.
	std::string query;
	/// Cost model referenced from the graph; must outlive the lifter.
	const CostModel &cost_model;
	/// Edit distance threshold the chain tables were built with; cells
	/// no transition realizes take DistMatrix::INF as their bound.
	int k = -1;
	/// Worker threads used for the chain table build (0 = auto).
	std::size_t num_threads = 0;
	/// Number of table levels: floor(log2(num_nodes)) + 1 (0 for an empty
	/// graph).
	std::size_t max_level = 0;
	/// up[u][t] is the node 2^t heavy steps from u, or empty if the chain
	/// ends before then.
	std::vector<std::vector<std::optional<node_id>>> up;
	/// up_mat[u][t] is the DistMatrix ed(u, 2^t, a, b) for every (u, t)
	/// touched by queries and built by the last build() call; untouched
	/// slots stay empty. decompose() sizes each node's vector to the top of
	/// its touched prefix, so build() fills slots in place.
	std::vector<std::vector<DistMatrix>> up_mat;
	/// touched[u][t] (mark-time) flags chain blocks the dry run marked;
	/// consumed into the per-level build lists by decompose().
	std::vector<std::vector<char>> touched;
	/// Per-level (ascending node order) list of blocks build() builds.
	std::vector<std::vector<std::size_t>> build_nodes;
};

} // namespace hlp_grep
