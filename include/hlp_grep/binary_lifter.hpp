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

#include <algorithm>
#include <cstddef>
#include <optional>
#include <stdexcept>
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
 */
class BinaryLifter {
public:
	/// Id of a graph node (mirrors POAGraph::node_id).
	using node_id = POAGraph::node_id;

	/**
	 * @brief Constructs the lifter and builds both tables.
	 *
	 * The node table has one row per graph node and one column per level
	 * 0 .. max_level - 1, where max_level = floor(log2(num_nodes)) + 1 is
	 * an upper bound on the useful power-of-two step counts. Built level
	 * by level from the per-node heavy outgoing edge pointers:
	 * up[u][t] = up[up[u][t - 1]][t - 1] whenever up[u][t - 1] is set, so
	 * no topological order is needed.
	 *
	 * The chain tables are built right after: level 0 holds the
	 * per-edge blocks of the actual heavy outgoing edges (the edge
	 * transition plus the destination node's label), and each higher
	 * level is the min-plus product of two adjacent level (t - 1)
	 * blocks sharing a chain midpoint.
	 *
	 * @param graph POAGraph whose heavy chains are indexed; must outlive
	 *              this object and must not be modified while it is alive.
	 * @param query Query string the chain distance tables index into.
	 * @param k     Edit distance threshold widening the windows.
	 */
	explicit BinaryLifter(const POAGraph &graph, const std::string &query,
	                      int k)
	    : graph(graph), query(query), k(k),
	      cost_model(graph.cost_model()) {
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

		// Per-level chain DistMatrix tables: level 0 holds the per-edge
		// blocks of the actual heavy outbound edges — the transition
		// across the edge (u, v) followed by v's whole label. Each higher
		// level is the min-plus product of two adjacent level (t - 1)
		// blocks sharing a chain midpoint.
		// Level t of node u is computed only when u's heavy_length is a
		// multiple of 2^t: the property is inherited by the two level
		// (t - 1) halves (their chain lengths are heavy_length(u) and
		// heavy_length(u) - 2^(t-1), both multiples of 2^(t-1)), so the
		// products below only ever read blocks that were built. A chain
		// of length L stores one block per power of two dividing one of
		// its nodes' heavy lengths (~2L blocks in total), which bounds
		// the whole table to O(V) matrices.
		up_mat.assign(
		    n, std::vector<DistMatrix>(max_level, DistMatrix(0, 0)));
		for (std::size_t u = 0; u < n; ++u)
			if (up[u][0].has_value())
				up_mat[u][0] = edge_matrix(u, *up[u][0]);
		for (std::size_t t = 1; t < max_level; ++t) {
			const std::size_t p = std::size_t{1} << t;
			for (std::size_t u = 0; u < n; ++u)
				if (up[u][t].has_value() &&
				    graph.node(u).heavy_length % p == 0)
					up_mat[u][t] = DistMatrix::min_plus_product(
					    up_mat[u][t - 1],
					    up_mat[*up[u][t - 1]][t - 1],
					    cost_model.is_monge());
		}
	}

	/**
	 * @brief Jumps l heavy steps from u (u^l; l == 0 is
	 *        the identity) and advances the caller's DP row across the
	 *        jumped chain.
	 *
	 * @p cost must hold one entry per allowed query position of the
	 * start node — its after-label window (see window()). It is advanced
	 * with DistMatrix::min_plus_apply
	 * over precomputed power-of-two chain blocks chosen greedily: at each
	 * node v with `rem` steps left, the block size is
	 * 2^min(ctz(heavy_length(v)), floor(log2(rem))) — the smaller of the
	 * lowest set bit of the node's heavy length and the highest set bit
	 * of the remaining step count. The first term guarantees the block
	 * was precomputed (level t exists exactly when 2^t divides the
	 * heavy length), the second keeps the walked total at l; the blocks
	 * apply in path order. Afterwards @p cost holds one entry per
	 * allowed position of the reached node. In effect, after the call
	 *   cost(b) = min_a cost_in(a) + ed(u, l, a, b)
	 * for the chain [u, u^1, ..., u^l] (each heavy edge crossing its
	 * destination node's full label). Cells no transition realizes carry
	 * DistMatrix::INF, which dominates any real cost while keeping every
	 * entry within the saturating product's overflow-free range.
	 *
	 * @param u    Node to start from.
	 * @param l    Number of heavy steps to walk.
	 * @param cost DP row of the start node; replaced by the advanced row
	 *             of the reached node.
	 * @return The node l heavy steps from u.
	 * @throws std::out_of_range if u is not a valid node id or the heavy
	 *         chain ends before l steps.
	 * @throws std::invalid_argument if cost.size() does not match the
	 *         start node's window.
	 */
	node_id jump(node_id u, int l, std::vector<int> &cost) const {
		if (u >= up.size())
			throw std::out_of_range("BinaryLifter::jump: invalid node id");
		node_id v = u;
		if (l == 0)
			return u;
		std::vector<int> row(cost);
		std::size_t rem = static_cast<std::size_t>(l);
		while (rem > 0) {
			if (!graph.node(v).heavy_neighbour.has_value())
				throw std::out_of_range(
				    "BinaryLifter::jump: heavy chain too short");
			// Block size 2^level: the smaller of the lowest set bit of
			// the node's heavy length (guarantees the block exists) and
			// the highest set bit of the remaining step count (keeps the
			// walked total at l).
			const std::size_t level =
			    std::min(static_cast<unsigned long long>(
			                 __builtin_ctzll(
			                     graph.node(v).heavy_length)),
			             static_cast<unsigned long long>(
			                 63 - __builtin_clzll(rem)));
			if (level >= max_level || !up[v][level].has_value())
				throw std::out_of_range(
				    "BinaryLifter::jump: heavy chain too short");
			row = up_mat[v][level].min_plus_apply(row,
			                                      cost_model.is_monge());
			v = *up[v][level];
			rem -= std::size_t{1} << level;
		}
		cost = std::move(row);
		return v;
	}

	/**
	 * @brief Steps a DP row across the single edge (u, v): applies the
	 *        edge transition plus v's whole label directly to the input
	 *        row and returns the advanced row.
	 *
	 * @p cost holds the allowed query positions at u (its after-label
	 * window, see window()); the returned row holds those at v. The row
	 * is computed by a semiglobal ED DP that superposes all source
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
	 * edge_matrix() materializes, without building the matrix, in
	 * O(|label| * (spread + 2k) + |cost|) time. Row |label|'s band is
	 * exactly v's window, so the output reads the swept cells directly.
	 *
	 * @param u    Source node of the edge.
	 * @param v    Destination node of the edge.
	 * @param cost DP row of u; must hold one entry per position of u's
	 *             window.
	 * @return The DP row at v, one entry per position of v's window.
	 * @throws std::invalid_argument if cost.size() does not match u's
	 *         window.
	 * @pre There is an edge (u, v) in the graph; not verified here.
	 */
	std::vector<int> step(node_id u, node_id v,
	                      const std::vector<int> &cost) const {
		const auto [lo_u, hi_u] = window(u);
		const auto [lo_v, hi_v] = window(v);
		if (cost.size() != static_cast<std::size_t>(hi_u - lo_u))
			throw std::invalid_argument(
			    "BinaryLifter::step: cost size does not match u's window");
		const std::string &label = graph.seq(v);
		std::vector<int> out(static_cast<std::size_t>(hi_v - lo_v),
		                     DistMatrix::INF);
		const auto [lL, hL] = window(v, static_cast<long>(label.size()));
		if (lL >= hL)
			return out; // final band empty: every crossing exceeds k
		const auto [l0, h0] = window(v, 0);
		std::vector<long> f(static_cast<std::size_t>(hL - l0),
		                    DistMatrix::INF);
		// Row 0: the source row clipped to the entry band, then chained
		// with insertions at v. Sources outside [l0, h0) die: below the
		// band every continuation leaves all bands, above it they start
		// past the query (possible only for clamped all-INF rows).
		for (long b = std::max(l0, lo_u); b < std::min(h0, hi_u); ++b)
			f[b - l0] = cost[b - lo_u];
		for (long b = l0 + 1; b < h0; ++b) {
			const long prev = f[b - 1 - l0];
			if (prev >= DistMatrix::INF)
				continue;
			long cur = prev + cost_model.ins(query[b - 1]);
			if (cur > DistMatrix::INF)
				cur = DistMatrix::INF;
			if (cur < f[b - l0])
				f[b - l0] = cur;
		}
		sweep_label(v, f, l0);
		for (long b = lo_v; b < hi_v; ++b)
			out[static_cast<std::size_t>(b - lo_v)] =
			    static_cast<int>(f[b - l0]);
		return out;
	}

private:
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
	 * transitions restricted to the band rectangle, so it is Monge
	 * whenever the cost model is (see CostModel::is_monge), with INF
	 * cells (b < a, seeds below the entry band, band-infeasible states)
	 * holding the inequalities trivially.
	 */
	DistMatrix edge_matrix(node_id u, node_id v) const {
		const auto [lo_u, hi_u] = window(u);
		const auto [lo_v, hi_v] = window(v);
		const std::string &label = graph.seq(v);
		DistMatrix mat(static_cast<std::size_t>(hi_u - lo_u),
		               static_cast<std::size_t>(hi_v - lo_v),
		               DistMatrix::INF);
		const auto [lL, hL] = window(v, static_cast<long>(label.size()));
		if (lL >= hL)
			return mat; // final band empty: every crossing exceeds k
		const auto [l0, h0] = window(v, 0);
		// DP row over [l0, hL); band j of the label occupies
		// window(v, j), and row L's band is exactly [lo_v, hi_v).
		std::vector<long> f(static_cast<std::size_t>(hL - l0),
		                    DistMatrix::INF);
		for (long a = lo_u; a < hi_u; ++a) {
			if (a < l0 || a >= h0)
				continue; // entry outside the seed band: no
				          // transition survives the bands
			std::fill(f.begin(), f.end(), DistMatrix::INF);
			// Row 0: seed at a, chain insertions within the entry band.
			f[a - l0] = 0;
			for (long b = a + 1; b < h0; ++b)
				f[b - l0] = std::min(
				    f[b - 1 - l0] + cost_model.ins(query[b - 1]),
				    static_cast<long>(DistMatrix::INF));
			sweep_label(v, f, l0);
			for (long b = std::max(a, lo_v); b < hi_v; ++b)
				mat(static_cast<std::size_t>(a - lo_u),
				    static_cast<std::size_t>(b - lo_v)) =
				    static_cast<int>(f[b - l0]);
		}
		return mat;
	}

	/**
	 * @brief Advances the DP row @p f across v's whole label, sweeping
	 *        each label character over its own band window(v, j).
	 *
	 * @p f holds one entry per position of [l0, hL) — the span from the
	 * entry band's bottom window(v, 0).lo to the after-label band's top
	 * window(v, |label|).hi — and its row 0 must already hold the entry
	 * states (seeds or source costs, insertion-chained within the entry
	 * band). Afterwards f holds the after-label states: row |label|'s
	 * band is exactly the node window [window(v).lo, window(v).hi).
	 * Cells a band never covers stay INF; they are below-band states
	 * whose continuations all exceed the threshold k.
	 */
	void sweep_label(node_id v, std::vector<long> &f, long l0) const {
		const std::string &label = graph.seq(v);
		for (long j = 1; j <= static_cast<long>(label.size()); ++j) {
			const auto [lj, hj] = window(v, j);
			const char c = label[j - 1];
			const int del_c = cost_model.del(c);
			// f(j - 1, lj - 1): the band bottom's diagonal source; INF
			// when the bands are clamped apart (no query char there).
			long diag = lj - 1 >= l0 ? f[lj - 1 - l0] : DistMatrix::INF;
			long left = DistMatrix::INF; // f(j, b - 1): below-band at lj
			for (long b = lj; b < hj; ++b) {
				const char qb = query[b - 1];
				const long up = f[b - l0]; // f(j - 1, b)
				long cur = up + del_c;
				if (left < DistMatrix::INF)
					cur = std::min(cur, left + cost_model.ins(qb));
				if (diag < DistMatrix::INF)
					cur = std::min(cur, diag + cost_model.consume(c, qb));
				if (cur > DistMatrix::INF)
					cur = DistMatrix::INF;
				f[b - l0] = cur;
				left = cur;
				diag = up;
			}
		}
	}

	/**
	 * @brief Half-open [lo, hi) window of query positions allowed at a
	 *        node under the stored query context: the after-label band
	 *        window(u, |label|) (see the offset overload), with the
	 *        lower end clamped to a non-empty range. The upper end
	 *        includes position |query| (all query characters consumed),
	 *        which Solver::query reads as the final row's last entry for
	 *        any path matched within k.
	 *
	 * The non-empty clamp keeps matrix dimensions positive for nodes
	 * whose feasible band lies entirely outside [0, |query|]: their rows
	 * carry valid (large) transition values only.
	 */
	std::pair<long, long> window(node_id u) const {
		const auto [lo, hi] =
		    window(u, static_cast<long>(graph.seq(u).size()));
		return {lo, std::max(hi, lo + 1)};
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
	std::pair<long, long> window(node_id u, long offset) const {
		const auto [jmin, jmax] = graph.pos_range(u);
		const long mk = static_cast<long>(k);
		const long m = static_cast<long>(query.size());
		const long lo =
		    std::max(0L, static_cast<long>(jmin) + offset - mk);
		const long hi =
		    std::min(m + 1, static_cast<long>(jmax) + offset + mk + 1);
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
	/// Number of table levels: floor(log2(num_nodes)) + 1 (0 for an empty
	/// graph).
	std::size_t max_level = 0;
	/// up[u][t] is the node 2^t heavy steps from u, or empty if the chain
	/// ends before then.
	std::vector<std::vector<std::optional<node_id>>> up;
	/// up_mat[u][t] is the DistMatrix ed(u, 2^t, a, b) for every (u, t)
	/// with a non-empty sub-chain.
	std::vector<std::vector<DistMatrix>> up_mat;
};

} // namespace hlp_grep
