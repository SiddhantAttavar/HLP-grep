/**
 * @file binary_lifter.hpp
 * @brief Binary lifting over the heavy chains of a POAGraph.
 *
 * Precomputes, for every graph node and every power of two, the node reached
 * by walking that many heavy edges. Alongside the node table, per-level
 * chain distances ed(u, 2^t, a, b) are precomputed as min-plus products of
 * level (t - 1) blocks (notes/idea.md:82). The chain tables depend on the
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
	 * single-edge blocks of the actual heavy outgoing edges
	 * (ed(u, 1, a, b)), and each higher level is the min-plus product of
	 * two adjacent level (t - 1) blocks sharing a chain midpoint.
	 *
	 * @param graph POAGraph whose heavy chains are indexed; must outlive
	 *              this object and must not be modified while it is alive.
	 * @param query Query string the chain distance tables index into.
	 * @param k     Edit distance threshold widening the windows.
	 */
	explicit BinaryLifter(const POAGraph &graph, const std::string &query,
	                      int k)
	    : graph(graph), query(query), k(k) {
		cost_model = graph.cost_model();
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

		// Per-level chain DistMatrix tables: level 0 holds the
		// single-edge blocks of the actual heavy outbound edges
		// (ed(u, 1, a, b)); each higher level is the min-plus product of
		// two adjacent level (t - 1) blocks sharing a chain midpoint.
		up_mat.assign(
		    n, std::vector<DistMatrix>(max_level, DistMatrix(0, 0)));
		for (std::size_t u = 0; u < n; ++u)
			if (up[u][0].has_value())
				up_mat[u][0] = edge_matrix(u, *up[u][0]);
		for (std::size_t t = 1; t < max_level; ++t)
			for (std::size_t u = 0; u < n; ++u)
				if (up[u][t].has_value())
					up_mat[u][t] = DistMatrix::min_plus_product(
					    up_mat[u][t - 1], up_mat[*up[u][t - 1]][t - 1]);
	}

	/**
	 * @brief Jumps l heavy steps from u (u^l in notes/idea.md; l == 0 is
	 *        the identity) and advances the caller's DP row across the
	 *        jumped chain.
	 *
	 * @p cost must hold one entry per allowed query position of the
	 * start node — its clamped pos_range window [j_min - k, j_max + k]
	 * over [0, |query|]. It is advanced with DistMatrix::min_plus_apply
	 * over the precomputed power-of-two chain blocks from the most
	 * significant bit of l down, so the blocks apply in path order.
	 * Afterwards @p cost holds one entry per allowed position of the
	 * reached node. In effect, after the call
	 *   cost(b) = min_a cost_in(a) + ed(u, l, a, b)
	 * for the chain [u, u^1, ..., u^l]. Entries stay within the block
	 * cost bound k + 1, so repeated applications cannot overflow.
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
	node_id jump(node_id u, std::size_t l, std::vector<int> &cost) const {
		if (u >= up.size())
			throw std::out_of_range("BinaryLifter::jump: invalid node id");
		node_id v = u;
		if (l == 0)
			return u;
		std::vector<int> row(cost);
		for (std::size_t t = max_level; t-- > 0;) {
			if (((l >> t) & 1) == 0)
				continue;
			if (!up[v][t].has_value())
				throw std::out_of_range(
				    "BinaryLifter::jump: heavy chain too short");
			row = up_mat[v][t].min_plus_apply(row);
			v = *up[v][t];
		}
		cost = std::move(row);
		return v;
	}

	/**
	 * @brief Steps a DP row across the single edge (u, v): applies the
	 *        one-edge transition (the same recurrence as the chain
	 *        tables' level-0 blocks) to @p cost on demand and returns the
	 *        advanced row.
	 *
	 * Rows are the allowed query positions at u, columns those at v, i.e.
	 * the pos_range windows [j_min - k, j_max + k] clamped to
	 * [0, |query|] (half-open). In effect, after the call
	 *   out(b) = min_a cost(a) + ed(u, 1, a, b)
	 * for the one-edge chain landing on v.
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
		return edge_matrix(u, v).min_plus_apply(cost);
	}

private:
	/**
	 * @brief The DistMatrix of the single edge (u, v): ed(u, 1, a, b)
	 *        between the one-edge chain landing on v and the query
	 *        interval [a, b] (notes/idea.md).
	 *
	 * Rows are the allowed query positions at u, columns those at v, i.e.
	 * the pos_range windows [j_min - k, j_max + k] clamped to
	 * [0, |query|] (half-open). Entry (a, b) is the min cost to go from
	 * DP state (u, a) to (v, b) across the edge, built left to right
	 * along b with the recurrence
	 *   mat(a, b) = min( mat(a, b - 1) + ins,
	 *                    del + ins * (b - a),
	 *                    match + ins * (b - a - 1) )
	 * where match is the running min of match(base(v), query[j]) over
	 * j in [a, b - 1]. Cells with b < a cannot be realized by any
	 * transition; they take the bound k + 1, an over-estimate of every
	 * cost a query matched within k can reach, so they never undercut
	 * real transitions in a min-plus composition.
	 */
	DistMatrix edge_matrix(node_id u, node_id v) const {
		const auto [lo_u, hi_u] = window(u);
		const auto [lo_v, hi_v] = window(v);
		const char base_v = graph.node(v).base;
		DistMatrix mat(static_cast<std::size_t>(hi_u - lo_u),
		               static_cast<std::size_t>(hi_v - lo_v), k + 1);
		for (long a = lo_u; a < hi_u; ++a) {
			int best_match = k + 1;
			int prev = k + 1;
			for (long b = std::max(lo_v, a); b < hi_v; ++b) {
				const int span = static_cast<int>(b - a);
				if (span > 0)
					best_match = std::min(
					    best_match,
					    cost_model.match(base_v, query[b - 1]));
				int cur;
				if (span > 0) {
					cur = std::min(
					    cost_model.del + cost_model.ins * span,
					    std::min(best_match + cost_model.ins * (span - 1),
					             prev + cost_model.ins));
				} else {
					cur = cost_model.del;
				}
				mat(static_cast<std::size_t>(a - lo_u),
				    static_cast<std::size_t>(b - lo_v)) = cur;
				prev = cur;
			}
		}
		return mat;
	}

	/**
	 * @brief Half-open [lo, hi) window of query positions allowed at a
	 *        node under the stored query context: [j_min - k, j_max + k + 1]
	 *        clamped to [0, |query|] (pos_range shifted by k and clamped).
	 *        The upper end includes position |query| (all query characters
	 *        consumed), which Solver::query reads as the final row's last
	 *        entry for any path matched within k.
	 */
	std::pair<long, long> window(node_id u) const {
		const auto [jmin, jmax] = graph.pos_range(u);
		const long mk = static_cast<long>(k);
		const long m = static_cast<long>(query.size());
		const long lo = std::max(0L, static_cast<long>(jmin) - mk);
		const long hi = std::min(m + 1, static_cast<long>(jmax) + mk + 2);
		// Nodes outside every threshold-feasible path (their pos_range
		// window lies entirely above |query|) would otherwise get an empty
		// window, which corrupts the chain-table products at construction.
		// Their rows carry valid (large) transition values only.
		return {lo, std::max(hi, lo + 1)};
	}

	/// Graph the lifter was built from; its heavy edges must not be
	/// modified while this object is alive.
	const POAGraph &graph;
	/// Query string the chain tables were built for.
	std::string query;
	/// Cost model snapshotted from the graph at construction.
	CostModel cost_model;
	/// Edit distance threshold the chain tables were built with; cells
	/// no transition realizes take k + 1 as their bound.
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
