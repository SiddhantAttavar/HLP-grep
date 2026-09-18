/**
 * @file poa_graph.hpp
 * @brief Partial order alignment (POA) graph built from a dictionary of DNA
 *        sequences.
 *
 * Constructs a POA graph incrementally: each dictionary sequence
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

/// Default base half-width for the adaptive build band `w = b + f * L`.
inline constexpr int DEFAULT_BAND_BASE = 10;
/// Default slope for the adaptive build band `w = b + f * L`.
inline constexpr double DEFAULT_BAND_SLOPE = 0.01;

/**
 * @brief A POA (partial order alignment) graph of a sequence dictionary.
 *
 * Nodes represent aligned bases, edges are classified heavy/light by
 * visiting frequencies, and each dictionary sequence is represented
 * as a path through the graph. The DAG is kept in a topological order
 * maintained by the graph.
 */
class POAGraph {
public:
	/// Id of a graph node.
	using node_id = std::size_t;

	/**
	 * @brief Heavy/light classification used by the heavy-light
	 *        decomposition of sequence paths.
	 */
	enum class EdgeType : unsigned char { HEAVY, LIGHT };

	/**
	 * @brief A step along a compressed sequence path.
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
			int length;
			node_id next;
		};
	};

	/**
	 * @brief A node of the graph: its base character and the range of
	 *        offsets at which it occurs among the sequence paths passing
	 *        through it (j_min / j_max). Offsets use
	 *        START as a "never on a path" sentinel.
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
	 * @brief Builds the POA graph from the given dictionary.
	 * Sequences are added in input order. The first sequence seeds the graph;
	 * each subsequent sequence is aligned to the current graph and merged.
	 * After the last sequence the heavy/light edge classification is finalized.
	 *
	 * The sequence-to-graph alignments use an adaptive abPOA-style band:
	 * per graph node only query columns in
	 * `[pos_min - w, pos_max + w]` (clamped to `[0, |s|]`) are computed,
	 * where `w = band_base + band_slope * |s|` and `pos_min/pos_max` are
	 * the incremental offsets of previously inserted paths. Narrow bands
	 * are heuristic: divergent sequences reuse less graph structure but
	 * their stored paths still spell the dictionary exactly.
	 *
	 * @param dict Dictionary of DNA sequences to index.
	 * @param cost Cost model used for the sequence-to-graph alignments.
	 *             Defaults to the unit-cost model. The referenced model
	 *             must outlive the graph.
	 * @param band_base Base half-width `b` of the adaptive build band.
	 * @param band_slope Slope `f` of the adaptive build band.
	 */
	explicit POAGraph(const std::vector<std::string> &dict,
	                  const CostModel &cost = DEFAULT_COST_MODEL,
	                  int band_base = DEFAULT_BAND_BASE,
	                  double band_slope = DEFAULT_BAND_SLOPE)
	    : cost(cost), band_base(band_base), band_slope(band_slope) {
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
	 * @brief Returns the cost model this graph was built with (used internally for
	 *        the sequence-to-graph alignments).
	 */
	const CostModel &cost_model() const {
		return cost;
	}

	/**
	 * @brief Read-only access to a node's full data: base character,
	 *        position range, and its heavy outgoing neighbour (empty if
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
	 * @brief The single heavy outgoing edge of a node, or START if the node
	 *        has no outgoing edges.
	 */
	node_id heavy_edge(node_id u) const {
		for (const Edge &e : out_edges[u])
			if (e.type == EdgeType::HEAVY)
				return e.neighbor;
		return START;
	}

	/**
	 * @brief Position range of a node: (min, max) 0-based offsets of the node
	 *        among the sequence paths passing through it, i.e. the
	 *        j_min/j_max pair used to bound query
	 *        precomputation windows to [j_min - k, j_max + k + 1] clamped
	 *        to [0, |query|].
	 *
	 * The offset of a node in a path is its 0-based index along that path.
	 * Both values equal SIZE_MAX for nodes on no stored path (an
	 * empty dictionary sequence).
	 */
	std::pair<std::size_t, std::size_t> pos_range(node_id u) const {
		return {nodes[u].pos_min, nodes[u].pos_max};
	}

	/**
	 * @brief Compressed representation of a sequence path: consecutive
	 *        heavy edges of the path merge into single HEAVY
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
				int len = 1;
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

	/// Virtual start node; edges from it are implicit and cost nothing.
	static constexpr node_id START = static_cast<node_id>(-1);

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

	const CostModel &cost; ///< Cost model of the graph alignments; must outlive the graph.
	/// Base half-width `b` of the adaptive build band `w = b + f * L`.
	int band_base = DEFAULT_BAND_BASE;
	/// Slope `f` of the adaptive build band `w = b + f * L`.
	double band_slope = DEFAULT_BAND_SLOPE;

	// -- construction ---------------------------------------------------------

	/** Builds the graph from the dictionary, in input order. */
	void build(const std::vector<std::string> &dict) {
		if (dict.empty())
			return;

		node_id prev = START;
		std::vector<node_id> seed;
		seed.reserve(dict[0].size());
		for (const char c : dict[0]) {
			const node_id u = append_node(c);
			add_edge(prev, u);
			seed.push_back(u);
			prev = u;
		}
		paths.push_back(std::move(seed));
		update_pos_ranges(paths.back());
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
		nodes.push_back({c, START, START});
		out_edges.emplace_back();
		in_edges.emplace_back();
		return u;
	}

	/**
	 * @brief Extends Node::pos_min / Node::pos_max with one stored path.
	 *
	 * For every offset `j` of @p path, node `path[j]` gets
	 * `pos_min = min(pos_min, j)` and `pos_max = max(pos_max, j)`.
	 * Nodes skipped by the new path keep their previous ranges.
	 */
	void update_pos_ranges(const std::vector<node_id> &path) {
		for (std::size_t j = 0; j < path.size(); ++j) {
			Node &n = nodes[path[j]];
			if (n.pos_min == START || j < n.pos_min)
				n.pos_min = j;
			if (n.pos_max == START || j > n.pos_max)
				n.pos_max = j;
		}
	}

	/**
	 * @brief Adaptive build band half-width for a sequence of length @p len:
	 *        `w = band_base + band_slope * len`, clamped to `>= 0`.
	 */
	long band_width(std::size_t len) const {
		const long w = static_cast<long>(band_base) +
		               static_cast<long>(band_slope * static_cast<double>(len));
		return w < 0 ? 0 : w;
	}

	/**
	 * @brief Inclusive band `[lo, hi]` of query columns computed for node @p u.
	 *
	 * `[pos_min - w, pos_max + w]` clamped to `[0, m]`. Nodes on no stored
	 * path yet (sentinel ranges) cover the full row. Empty when the node
	 * range lies entirely outside the clamped interval; callers skip those
	 * rows (their cells stay INF).
	 */
	std::pair<long, long> band_bounds(node_id u, long m, long w) const {
		const Node &n = nodes[u];
		if (n.pos_min == START || n.pos_max == START)
			return {0, m};
		long lo = static_cast<long>(n.pos_min) - w;
		long hi = static_cast<long>(n.pos_max) + w;
		if (lo < 0)
			lo = 0;
		if (hi > m)
			hi = m;
		return {lo, hi};
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
		if (u == START)
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
	 * @brief Aligns a sequence to the current graph under the adaptive band.
	 *
	 * Minimum-cost global alignment of @p s against the DAG restricted to
	 * the per-node bands `[pos_min - w, pos_max + w]` clamped to `[0, |s|]`,
	 * where `w = band_base + band_slope * |s|` and the ranges are the
	 * incremental offsets of previously inserted paths. DP over
	 * (topological position, prefix length) with a virtual start row, using
	 * insertion/deletion/substitution costs from cost. Cost ties prefer
	 * exact consumes, then deletions, then insertions. Out-of-band cells
	 * stay INF; the band is heuristic, so divergent sequences simply reuse
	 * less graph structure.
	 *
	 * @param s Sequence to align.
	 * @return Alignment record: cost and alignment steps, to be merged with
	 *         add_alignment().
	 */
	Alignment align_to_graph(const std::string &s) {
		const std::size_t V = nodes.size();
		const std::size_t m = s.size();
		const long ml = static_cast<long>(m);
		constexpr int INF = std::numeric_limits<int>::max() / 2;
		const std::size_t R = V + 1; // DP rows; last row = virtual start.
		const node_id vr = static_cast<node_id>(V);
		const node_id INS_CODE =
		    static_cast<unsigned char>(AlignmentOp::Type::INSERT) + 1;
		const node_id DEL_CODE =
		    static_cast<unsigned char>(AlignmentOp::Type::DELETE) + 1;
		const node_id CON_CODE =
		    static_cast<unsigned char>(AlignmentOp::Type::CONSUME) + 1;

		std::vector<node_id> row_of(V, 0); // node id -> topological row.
		for (std::size_t p = 0; p < V; ++p)
			row_of[topo[p]] = p;

		// Inclusive per-row bands over columns 0..m; the virtual start row
		// is unbanded. Rows whose band is empty store nothing.
		const long w = band_width(m);
		std::vector<long> lo(R, 1), hi(R, 0);
		lo[vr] = 0;
		hi[vr] = ml;
		for (std::size_t p = 0; p < V; ++p) {
			const auto [blo, bhi] = band_bounds(topo[p], ml, w);
			lo[p] = blo;
			hi[p] = bhi;
		}

		// Banded row-local storage: row p holds columns [lo[p], hi[p]] back
		// to back at row_off[p]; cell (p, i) is at row_off[p] + (i - lo[p]).
		// Total cells are O(V * w) instead of O(V * m). Every access below
		// is guarded by the same band-membership checks as the compute.
		std::vector<std::size_t> row_off(R + 1, 0);
		for (std::size_t p = 0; p < R; ++p)
			row_off[p + 1] =
			    row_off[p] +
			    (hi[p] >= lo[p]
			         ? static_cast<std::size_t>(hi[p] - lo[p] + 1)
			         : 0);
		std::vector<int> dp(row_off[R], INF);
		std::vector<unsigned char> op(row_off[R], 0);
		std::vector<node_id> from(row_off[R], START);

		// Start row: dp[vr][j] = sum of insertion costs over s[0..j).
		dp[row_off[vr]] = 0;
		for (std::size_t j = 1; j <= m; ++j) {
			const std::size_t k = row_off[vr] + j;
			dp[k] = dp[k - 1] + cost.ins(s[j - 1]);
			op[k] = INS_CODE;
			from[k] = vr;
		}
		// dp[u][0]: delete graph bases between the start and node u.
		for (std::size_t p = 0; p < V; ++p) {
			if (lo[p] > 0 || hi[p] < 0)
				continue; // column 0 out of band
			int d = INF;
			if (in_edges[topo[p]].empty()) {
				d = cost.del(nodes[topo[p]].base);
			} else {
				for (const Edge &e : in_edges[topo[p]]) {
					const node_id prow = row_of[e.neighbor];
					if (lo[prow] > 0 || hi[prow] < 0)
						continue; // predecessor has no column 0
					const int c =
					    dp[row_off[prow]] +
					    cost.del(nodes[topo[p]].base);
					if (c < d)
						d = c;
				}
			}
			if (d > INF)
				d = INF;
			dp[row_off[p]] = d;
		}

		// Main banded DP over (topological row, prefix length): only columns
		// i in [lo[p], hi[p]] are computed; transitions from out-of-band
		// cells are skipped (they hold INF).
		for (std::size_t p = 0; p < V; ++p) {
			const node_id u = topo[p];
			const long i_lo = std::max(lo[p], 1L);
			for (long il = i_lo; il <= hi[p]; ++il) {
				const std::size_t i = static_cast<std::size_t>(il);
				const std::size_t k =
				    row_off[p] + static_cast<std::size_t>(il - lo[p]);
				int best = INF;
				unsigned char bop = 0;
				node_id bfrom = START;
				int bpri = -1;

				if (il - 1 >= lo[p]) {
					best = dp[k - 1] + cost.ins(s[i - 1]); // insert s[i-1] at u
					bop = INS_CODE;
					bfrom = p;
					bpri = 0;
				}

				auto consider_del = [&](node_id prow) {
					if (prow != vr && (il < lo[prow] || il > hi[prow]))
						return;
					const int c =
					    dp[row_off[prow] +
					       static_cast<std::size_t>(il - lo[prow])] +
					    cost.del(nodes[u].base);
					if (c < best || (c == best && bpri < 1)) {
						best = c;
						bop = DEL_CODE;
						bfrom = prow;
						bpri = 1;
					}
				};
				auto consider_con = [&](node_id prow) {
					if (prow != vr && (il - 1 < lo[prow] || il - 1 > hi[prow]))
						return;
					const int c =
					    dp[row_off[prow] +
					       static_cast<std::size_t>(il - 1 - lo[prow])] +
					    cost.consume(nodes[u].base, s[i - 1]);
					int pri = 2;
					if (nodes[u].base == s[i - 1] &&
					    cost.consume(nodes[u].base, s[i - 1]) == 0)
						pri = 3;
					if (c < best || (c == best && bpri < pri)) {
						best = c;
						bop = CON_CODE;
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
				if (best > INF)
					best = INF;
				dp[k] = best;
				op[k] = bop;
				from[k] = bfrom;
			}
		}

		// Answer: min over all end nodes holding column m (graph may end
		// anywhere). Cost ties prefer the same ordering as inside a cell,
		// so equal-cost paths reuse existing graph structure instead of
		// spawning new nodes. Rows without column m are skipped; when no
		// row reaches m (band drift on divergent sequences) the sequence
		// is inserted as fresh nodes.
		std::size_t br = vr;
		int bc = dp[row_off[vr] + m];
		int bpri = -1;
		for (std::size_t p = 0; p < V; ++p) {
			if (ml < lo[p] || ml > hi[p])
				continue;
			const std::size_t km =
			    row_off[p] + static_cast<std::size_t>(ml - lo[p]);
			const int c = dp[km];
			const node_id u = topo[p];
			int pri = 0;
			switch (op[km]) {
			case 2:
				pri = 1;
				break;
			case 3:
				pri = (m > 0 && nodes[u].base == s[m - 1]) ? 3 : 2;
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

		// Degenerate band miss: no end node reaches column m.
		if (bc >= INF) {
			Alignment al;
			al.cost = bc;
			al.ops.reserve(m);
			for (std::size_t j = 0; j < m; ++j)
				al.ops.push_back({AlignmentOp::Type::INSERT, START, s[j]});
			return al;
		}

		// Traceback parents to the empty-prefix column (rows, not node ids).
		// op == 0 marks an uncomputed cell on a heuristic path: insert the
		// query character instead of following a bogus parent.
		Alignment al;
		al.cost = bc;
		std::vector<AlignmentOp> ops;
		std::size_t r = br;
		for (std::size_t i = m; i > 0;) {
			const std::size_t k =
			    row_off[r] + (i - static_cast<std::size_t>(lo[r]));
			const node_id row2node = (r == vr) ? START : topo[r];
			switch (op[k]) {
			case 1:
				ops.push_back({AlignmentOp::Type::INSERT, START, s[i - 1]});
				--i;
				break;
			case 2:
				ops.push_back({AlignmentOp::Type::DELETE, row2node, '\0'});
				r = from[k];
				break;
			case 3:
				ops.push_back({AlignmentOp::Type::CONSUME, row2node, s[i - 1]});
				r = from[k];
				--i;
				break;
			default:
				ops.push_back({AlignmentOp::Type::INSERT, START, s[i - 1]});
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
		node_id prev = START;
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
		update_pos_ranges(paths.back());
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
	 *        stored sequence paths (j_min / j_max).
	 *        Offsets are 0-based positions along each path.
	 */
	void compute_pos_ranges() {
		for (Node &n : nodes) {
			n.pos_min = START;
			n.pos_max = START;
		}
		for (const auto &path : paths)
			for (std::size_t j = 0; j < path.size(); ++j) {
				Node &n = nodes[path[j]];
				if (n.pos_min == START || j < n.pos_min)
					n.pos_min = j;
				if (n.pos_max == START || j > n.pos_max)
					n.pos_max = j;
			}
	}
};

} // namespace hlp_grep
