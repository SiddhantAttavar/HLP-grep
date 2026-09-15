/**
 * @file poa_graph.hpp
 * @brief Partial order alignment (POA) graph built from a dictionary of DNA
 *        sequences.
 *
 * Constructs a pangenome-style graph incrementally: each dictionary sequence
 * is aligned to the current graph (min-cost alignment under CostModel) and its
 * alignment is merged into the graph. Every dictionary sequence is stored as
 * an explicit path of node ids, as required by heavy-light path compression.
 */
#pragma once
#include <hlp_grep/cost_model.hpp>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace hlp_grep {

/**
 * @brief A POA (partial order alignment) graph of a sequence dictionary.
 *
 * Nodes represent aligned bases, edges carry visiting frequencies, and each
 * dictionary sequence is represented as a path through the graph. The DAG is
 * kept in a topological order maintained by the graph.
 */
class POAGraph {
public:
	/// Id of a graph node.
	using node_id = std::size_t;

	/**
	 * @brief Heavy/light classification used by the heavy-light
	 *        decomposition of sequence paths (notes/idea.md).
	 */
	enum class EdgeType : unsigned char { HEAVY, LIGHT };

	/**
	 * @brief A step along a compressed sequence path (notes/idea.md).
	 *
	 * Either a HEAVY chain advancing `length` heavy edges (all from the same
	 * chain), or a single LIGHT edge landing on node `next`. The two variants
	 * share storage in the union; `type` discriminates them.
	 */
	struct CompressedEdge {
		EdgeType type = EdgeType::LIGHT;

		/// HEAVY: number of heavy edges advanced along the chain.
		/// LIGHT: destination node of the single light edge.
		union {
			node_id length;
			node_id next;
		};
	};

	/**
	 * @brief Per-edge metadata stored in the graph: neighbor node, its
	 *        heavy/light classification (finalized by mark_heavy_edges()),
	 *        and its unique edge id (shared by the out-edge and in-edge
	 *        records of the same edge). Visiting frequencies are temporary
	 *        values and are not stored in the struct.
	 */
	struct Edge {
		EdgeType type = EdgeType::LIGHT; ///< Heavy or light edge.
		node_id neighbor = 0;            ///< Other end of the edge.
		std::size_t edge_id = 0;         ///< Unique id assigned at creation.
	};

	/**
	 * @brief A node of the graph: its base character and the range of
	 *        offsets at which it occurs among the sequence paths passing
	 *        through it (j_min / j_max in notes/idea.md). Offsets use
	 *        kStart as a "never on a path" sentinel.
	 */
	struct Node {
		char base = 0;            ///< Base character represented by this node.
		std::size_t pos_min = 0;  ///< Minimum path offset through this node.
		std::size_t pos_max = 0;  ///< Maximum path offset through this node.
		/// Destination node of the heaviest outgoing edge, or empty if the
		/// node has no outgoing edges. Set by mark_heavy_edges() at the end
		/// of build; valid only while the node's out-edge list is not
		/// modified.
		std::optional<node_id> heavy_neighbour;
	};

	/**
	 * @brief Outgoing edges of a node.
	 *
	 * @param u Node id.
	 * @return The list of edges leaving @p u (empty for a sink node).
	 */
	const std::vector<Edge> &outgoing_edges(node_id u) const {
		return out_edges[u];
	}

	/**
	 * @brief Builds the POA graph from the given dictionary.
	 * Sequences are added in input order. The first sequence seeds the graph;
	 * each subsequent sequence is aligned to the current graph and merged.
	 * After the last sequence the heavy/light edge classification is finalized.
	 *
	 * @param dict Dictionary of DNA sequences to index.
	 * @param cost Cost model used for the sequence-to-graph alignments.
	 *             Defaults to the unit-cost model (`CostModel{}`).
	 */
	explicit POAGraph(const std::vector<std::string> &dict,
	                  const CostModel &cost = {})
	    : cost(cost) {
		build(dict);
	}

	/**
	 * @brief Returns the number of nodes in the graph.
	 */
	std::size_t num_nodes() const {
		return nodes.size();
	}

	/**
	 * @brief Returns the base character represented by a node.
	 *
	 * @param u Node id.
	 */
	char base(node_id u) const {
		return nodes[u].base;
	}

	/**
	 * @brief The cost model this graph was built with (used internally for
	 *        the sequence-to-graph alignments).
	 */
	const CostModel &cost_model() const {
		return cost;
	}

	/**
	 * @brief Read-only access to a node's full data: base character,
	 *        position range, and its heavy outgoing Edge pointer (null if
	 *        the node has no outgoing edges).
	 *
	 * @param u Node id.
	 */
	const Node &node(node_id u) const {
		return nodes[u];
	}

	/**
	 * @brief Number of dictionary sequences inserted into the graph.
	 */
	std::size_t num_sequences() const {
		return paths.size();
	}

	/**
	 * @brief Returns the graph path of a dictionary sequence.
	 *
	 * Consecutive nodes are connected by graph edges. Node ids
	 * need not be monotonic along a path.
	 *
	 * @param seq Index of the sequence in the dictionary as passed to the
	 *            constructor (0-based by convention, matching Result::id).
	 */
	const std::vector<node_id> &path(std::size_t seq) const {
		return paths[seq];
	}

	/**
	 * @brief Temporary visiting frequency of edge (u, v): the number of
	 *        stored sequence paths traversing it, or 0 if the edge is
	 *        never traversed. Weights are not stored in the edges; this
	 *        counts on demand and costs O(total path length) per call.
	 */
	std::size_t edge_weight(node_id u, node_id v) const {
		std::size_t weight = 0;
		for (const auto &path : paths)
			for (std::size_t i = 1; i < path.size(); ++i)
				if (path[i - 1] == u && path[i] == v)
					++weight;
		return weight;
	}

	/**
	 * @brief Heavy/light classification of edge (u, v); LIGHT for absent
	 *        edges.
	 */
	EdgeType edge_type(node_id u, node_id v) const {
		if (const Edge *e = find_edge(out_edges[u], v))
			return e->type;
		return EdgeType::LIGHT;
	}

	/**
	 * @brief The single heavy outgoing edge of a node, or kStart if the node
	 *        has no outgoing edges.
	 */
	node_id heavy_edge(node_id u) const {
		for (const Edge &e : out_edges[u])
			if (e.type == EdgeType::HEAVY)
				return e.neighbor;
		return kStart;
	}

	/**
	 * @brief Position range of a node: (min, max) 0-based offsets of the node
	 *        among the sequence paths passing through it, i.e. the
	 *        j_min/j_max pair used by notes/idea.md to bound query
	 *        precomputation windows to [j_min - k, j_max + k].
	 *
	 * The offset of a node in a path is its 0-based index along that path.
	 * Both values equal kNone (SIZE_MAX) for nodes on no stored path (an
	 * empty dictionary sequence).
	 */
	std::pair<std::size_t, std::size_t> pos_range(node_id u) const {
		return {nodes[u].pos_min, nodes[u].pos_max};
	}

	/**
	 * @brief Compressed representation of a sequence path (notes/idea.md):
	 *        consecutive heavy edges of the path merge into single HEAVY
	 *        steps (`length` = number of heavy edges advanced); every light
	 *        edge becomes a LIGHT step (`next` = the node it lands on, so
	 *        the walk does not depend on unique heavy-chain successors).
	 *
	/**
	 * @brief Compressed representation of a sequence path (notes/idea.md):
	 *        consecutive heavy edges of the path merge into single HEAVY
	 *        steps (`length` = number of heavy edges advanced); every light
	 *        edge becomes a LIGHT step (`next` = the node it lands on, so
	 *        the walk does not depend on unique heavy-chain successors).
	 */
	struct CompressedPath {
		node_id start = 0;                 ///< Node the walk starts from (`path(seq).front()`).
		std::vector<CompressedEdge> steps; ///< Steps, in path order.
	};

	CompressedPath compressed_path(std::size_t seq) const {
		CompressedPath out;
		const auto &p = paths[seq];
		out.start = p.front();
		out.steps.reserve(p.size());
		std::size_t i = 1;
		while (i < p.size()) {
			CompressedEdge step;
			if (edge_type(p[i - 1], p[i]) == EdgeType::HEAVY) {
				step.type = EdgeType::HEAVY;
				node_id len = 1;
				++i;
				while (i < p.size() &&
				       edge_type(p[i - 1], p[i]) == EdgeType::HEAVY) {
					++len;
					++i;
				}
				step.length = len;
			} else {
				step.next = p[i];
				++i;
			}
			out.steps.push_back(step);
		}
		return out;
	}

private:
	/**
	 * @brief One step of a sequence-to-graph alignment.
	 */
	struct AlignmentOp {
		enum class Type : unsigned char { INSERT, DELETE, CONSUME };

		Type type;
		node_id node; ///< Graph node referenced by DELETE/CONSUME.
		char ch;      ///< Inserted or consumed character.
	};

	/**
	 * @brief Result of a sequence-to-graph alignment.
	 */
	struct Alignment {
		int cost = 0;                 ///< Minimum alignment cost.
		std::vector<AlignmentOp> ops; ///< Alignment steps, in sequence order.
	};

	/// Virtual start node; edges from it are implicit and cost nothing.
	static constexpr node_id kStart = static_cast<node_id>(-1);

	/// Per node: base character and position range.
	std::vector<Node> nodes;
	/// Outgoing edges per node.
	std::vector<std::vector<Edge>> out_edges;
	/// Incoming edges per node.
	std::vector<std::vector<Edge>> in_edges;
	/// Id assigned to the next created edge.
	std::size_t next_edge_id = 0;
	/// Per dictionary sequence: node ids visited by that sequence.
	std::vector<std::vector<node_id>> paths;
	/// Node ids in topological order of the DAG.
	std::vector<node_id> topo;

	CostModel cost;

	// -- construction ---------------------------------------------------------

	/** Builds the graph from the dictionary, in input order. */
	void build(const std::vector<std::string> &dict) {
		if (dict.empty())
			return;

		node_id prev = kStart;
		std::vector<node_id> seed;
		seed.reserve(dict[0].size());
		for (const char c : dict[0]) {
			const node_id u = append_node(c);
			add_edge(prev, u);
			seed.push_back(u);
			prev = u;
		}
		paths.push_back(std::move(seed));
		rebuild_topo();

		for (std::size_t id = 1; id < dict.size(); ++id) {
			const Alignment al = align_to_graph(dict[id]);
			add_alignment(al.ops);
			rebuild_topo();
		}

		mark_heavy_edges();
		compute_pos_ranges();
	}

	/** Creates a new node with base @p c at the end of the graph. */
	node_id append_node(char c) {
		const node_id u = nodes.size();
		nodes.push_back({c, kStart, kStart});
		out_edges.emplace_back();
		in_edges.emplace_back();
		return u;
	}

	/** Finds an existing edge in an adjacency list, or returns null. */
	static Edge *find_edge(std::vector<Edge> &edges, node_id v) {
		for (Edge &e : edges)
			if (e.neighbor == v)
				return &e;
		return nullptr;
	}

	/** Const overload of find_edge(). */
	static const Edge *find_edge(const std::vector<Edge> &edges, node_id v) {
		for (const Edge &e : edges)
			if (e.neighbor == v)
				return &e;
		return nullptr;
	}

	/**
	 * @brief Creates edge (u, v) with a fresh unique edge id if missing.
	 *        Both the out-edge and in-edge records share the id. Edges
	 *        from the virtual start node are implicit and therefore not
	 *        stored. Visiting frequencies are temporary values and are
	 *        computed on demand, not incremented here.
	 */
	void add_edge(node_id u, node_id v) {
		if (u == kStart)
			return;
		if (!find_edge(out_edges[u], v)) {
			const std::size_t id = next_edge_id++;
			out_edges[u].push_back({EdgeType::LIGHT, v, id});
			in_edges[v].push_back({EdgeType::LIGHT, u, id});
		}
	}

	/**
	 * @brief Recomputes topo (Kahn's algorithm), deterministically.
	 *
	 * Alignment-merged graphs are not topologically ordered by node id alone
	 * (mid-path insertions create edges like u_low -> u_new -> u_mid), so a
	 * real order is maintained for DP-based alignments and for clients.
	 */
	void rebuild_topo() {
		const std::size_t V = nodes.size();
		std::vector<std::size_t> indeg(V, 0);
		for (node_id u = 0; u < V; ++u)
			for (const Edge &e : out_edges[u])
				indeg[e.neighbor]++;
		topo.clear();
		topo.reserve(V);
		for (node_id u = 0; u < V; ++u) // stable source order
			if (indeg[u] == 0)
				topo.push_back(u);
		for (std::size_t h = 0; h < topo.size(); ++h) {
			const node_id u = topo[h];
			for (const Edge &e : out_edges[u]) {
				std::size_t &d = indeg[e.neighbor];
				if (--d == 0)
					topo.push_back(e.neighbor);
			}
		}
	}

	/**
	 * @brief Aligns a sequence to the current graph.
	 *
	 * Minimum-cost global alignment of @p s against the DAG: DP over
	 * (topological position, prefix length) with a virtual start row, using
	 * insertion/deletion/substitution costs from cost. Cost ties prefer
	 * exact consumes, then deletions, then insertions.
	 *
	 * @param s Sequence to align.
	 * @return Alignment record: cost and alignment steps, to be merged with
	 *         add_alignment().
	 */
	Alignment align_to_graph(const std::string &s) {
		const std::size_t V = nodes.size();
		const std::size_t m = s.size();
		constexpr int kInf = std::numeric_limits<int>::max() / 2;
		const std::size_t R = V + 1; // DP rows; last row = virtual start.
		const std::size_t C = m + 1; // DP columns.
		const node_id vr = static_cast<node_id>(V);
		const node_id kInsCode =
		    static_cast<unsigned char>(AlignmentOp::Type::INSERT) + 1;
		const node_id kDelCode =
		    static_cast<unsigned char>(AlignmentOp::Type::DELETE) + 1;
		const node_id kConCode =
		    static_cast<unsigned char>(AlignmentOp::Type::CONSUME) + 1;

		std::vector<int> dp(R * C, kInf);
		std::vector<unsigned char> op(R * C, 0);
		std::vector<node_id> from(R * C, kStart);
		std::vector<node_id> row_of(V, 0); // node id -> topological row.
		for (std::size_t p = 0; p < V; ++p)
			row_of[topo[p]] = p;

		// Start row: dp[vr][j] = j * ins (insert-only prefix).
		dp[vr * C] = 0;
		for (std::size_t j = 1; j <= m; ++j) {
			const std::size_t k = vr * C + j;
			dp[k] = dp[k - 1] + cost.ins;
			op[k] = kInsCode;
			from[k] = vr;
		}
		// dp[u][0]: delete graph bases between the start and node u.
		for (std::size_t p = 0; p < V; ++p) {
			int d = kInf;
			if (in_edges[topo[p]].empty()) {
				d = cost.del;
			} else {
				for (const Edge &e : in_edges[topo[p]]) {
					const int c = dp[row_of[e.neighbor] * C] + cost.del;
					if (c < d)
						d = c;
				}
			}
			dp[p * C] = d;
		}

		// Main DP over (topological row, prefix length).
		for (std::size_t p = 0; p < V; ++p) {
			const node_id u = topo[p];
			for (std::size_t i = 1; i <= m; ++i) {
				const std::size_t k = p * C + i;
				int best = dp[k - 1] + cost.ins; // insert s[i-1] at u
				unsigned char bop = kInsCode;
				node_id bfrom = p;
				int bpri = 0;

				auto consider_del = [&](node_id prow) {
					const int c = dp[prow * C + i] + cost.del;
					if (c < best || (c == best && bpri < 1)) {
						best = c;
						bop = kDelCode;
						bfrom = prow;
						bpri = 1;
					}
				};
				auto consider_con = [&](node_id prow) {
					const int c =
					    dp[prow * C + i - 1] + cost.match(nodes[u].base, s[i - 1]);
					int pri = 2;
					if (nodes[u].base == s[i - 1] &&
					    cost.match(nodes[u].base, s[i - 1]) == 0)
						pri = 3;
					if (c < best || (c == best && bpri < pri)) {
						best = c;
						bop = kConCode;
						bfrom = prow;
						bpri = pri;
					}
				};

				if (in_edges[u].empty()) {
					consider_del(vr);
					consider_con(vr);
				} else {
					for (const Edge &e : in_edges[u]) {
						consider_del(row_of[e.neighbor]);
						consider_con(row_of[e.neighbor]);
					}
				}
				dp[k] = best;
				op[k] = bop;
				from[k] = bfrom;
			}
		}

		// Answer: min over all end nodes (graph may end anywhere). Cost ties
		// prefer the same ordering as inside a cell, so equal-cost paths reuse
		// existing graph structure instead of spawning new nodes.
		std::size_t br = vr;
		int bc = dp[vr * C + m];
		int bpri = -1;
		for (std::size_t p = 0; p < V; ++p) {
			const int c = dp[p * C + m];
			const node_id u = topo[p];
			int pri = 0;
			switch (op[p * C + m]) {
			case 2:
				pri = 1;
				break;
			case 3:
				pri = nodes[u].base == s[m - 1] ? 3 : 2;
				break;
			default:
				break;
			}
			if (c < bc || (c == bc && pri > bpri)) {
				bc = c;
				br = p;
				bpri = pri;
			}
		}

		// Traceback parents to the empty-prefix column (rows, not node ids).
		Alignment al;
		al.cost = bc;
		std::vector<AlignmentOp> ops;
		std::size_t r = br;
		for (std::size_t i = m; i > 0;) {
			const std::size_t k = r * C + i;
			const node_id row2node = (r == vr) ? kStart : topo[r];
			switch (op[k]) {
			case 1:
				ops.push_back({AlignmentOp::Type::INSERT, kStart, s[i - 1]});
				--i;
				break;
			case 2:
				ops.push_back({AlignmentOp::Type::DELETE, row2node, '\0'});
				r = from[k];
				break;
			default:
				ops.push_back({AlignmentOp::Type::CONSUME, row2node, s[i - 1]});
				r = from[k];
				--i;
				break;
			}
		}
		std::reverse(ops.begin(), ops.end());
		al.ops = std::move(ops);
		return al;
	}

	/**
	 * @brief Merges a sequence's alignment into the graph.
	 *
	 * Exact-match consumption reuses nodes and reuses existing edges;
	 * substitutions and insertions create new nodes. Deletions are skipped
	 * and do not touch edges. Appends the resulting path to paths and
	 * refreshes the topological order.
	 */
	void add_alignment(const std::vector<AlignmentOp> &ops) {
		std::vector<node_id> path;
		path.reserve(ops.size());
		node_id prev = kStart;
		for (const AlignmentOp &a : ops) {
			switch (a.type) {
			case AlignmentOp::Type::DELETE:
				break; // the new path skips this graph node entirely
			case AlignmentOp::Type::INSERT: {
				const node_id u = append_node(a.ch);
				add_edge(prev, u);
				path.push_back(u);
				prev = u;
				break;
			}
			case AlignmentOp::Type::CONSUME: {
				const node_id u = (nodes[a.node].base == a.ch) ? a.node
				                                           : append_node(a.ch);
				add_edge(prev, u);
				path.push_back(u);
				prev = u;
				break;
			}
			}
		}
		paths.push_back(std::move(path));
		rebuild_topo();
	}

	/**
	 * @brief Classifies each node's outgoing edges from temporary visiting
	 *        frequencies.
	 *
	 * Per-edge traversal counts are recomputed from the stored sequence
	 * paths (temporary values, not stored in the edges). For every node
	 * with outgoing edges, the edge with the highest count becomes HEAVY
	 * (exactly one per node); ties are broken by the largest destination
	 * node id. All other outgoing edges stay LIGHT.
	 */
	void mark_heavy_edges() {
		std::vector<std::size_t> weight(next_edge_id, 0);
		for (const auto &path : paths)
			for (std::size_t i = 1; i < path.size(); ++i)
				weight[find_edge(out_edges[path[i - 1]], path[i])
				           ->edge_id]++;

		for (std::size_t u = 0; u < out_edges.size(); ++u) {
			auto &edges = out_edges[u];
			Edge *heavy = nullptr;
			for (Edge &e : edges) {
				e.type = EdgeType::LIGHT;
				if (!heavy || weight[e.edge_id] > weight[heavy->edge_id] ||
				    (weight[e.edge_id] == weight[heavy->edge_id] &&
				     e.neighbor > heavy->neighbor))
					heavy = &e;
			}
			if (heavy)
				heavy->type = EdgeType::HEAVY;
			nodes[u].heavy_neighbour =
			    heavy ? std::optional<node_id>(heavy->neighbor)
			          : std::nullopt;
		}
	}

	/**
	 * @brief Computes Node::pos_min / Node::pos_max: for every node, the
	 *        minimum and maximum offset at which the node occurs among the
	 *        stored sequence paths (j_min / j_max in notes/idea.md).
	 *        Offsets are 0-based positions along each path.
	 */
	void compute_pos_ranges() {
		for (Node &n : nodes) {
			n.pos_min = kStart;
			n.pos_max = kStart;
		}
		for (const auto &path : paths)
			for (std::size_t j = 0; j < path.size(); ++j) {
				Node &n = nodes[path[j]];
				if (n.pos_min == kStart || j < n.pos_min)
					n.pos_min = j;
				if (n.pos_max == kStart || j > n.pos_max)
					n.pos_max = j;
			}
	}
};

} // namespace hlp_grep
