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
#include <string>
#include <utility>
#include <vector>

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
	 * @brief Builds the POA graph from the given dictionary.
	 *
	 * Sequences are added in input order. The first sequence seeds the graph;
	 * each subsequent sequence is aligned to the current graph and merged.
	 *
	 * @param dict Dictionary of DNA sequences to index.
	 * @param cost Cost model used for the sequence-to-graph alignments.
	 *             Defaults to the unit-cost model (`CostModel{}`).
	 */
	explicit POAGraph(const std::vector<std::string> &dict,
	                  const CostModel &cost = {})
	    : cost_(cost) {
		build(dict);
	}

	/**
	 * @brief Returns the number of nodes in the graph.
	 */
	std::size_t num_nodes() const {
		return bases_.size();
	}

	/**
	 * @brief Returns the base character represented by a node.
	 *
	 * @param u Node id.
	 */
	char base(node_id u) const {
		return bases_[u];
	}

	/**
	 * @brief Number of dictionary sequences inserted into the graph.
	 */
	std::size_t num_sequences() const {
		return paths_.size();
	}

	/**
	 * @brief Returns the graph path of a dictionary sequence.
	 *
	 * Consecutive nodes are connected by edges of positive weight. Node ids
	 * need not be monotonic along a path.
	 *
	 * @param seq Index of the sequence in the dictionary as passed to the
	 *            constructor (0-based by convention, matching Result::id).
	 */
	const std::vector<node_id> &path(std::size_t seq) const {
		return paths_[seq];
	}

	/**
	 * @brief Visiting frequency of edge (u, v), or 0 if the edge is absent.
	 */
	std::size_t edge_weight(node_id u, node_id v) const {
		if (const Edge *e = find_edge(out_edges_[u], v))
			return e->weight;
		return 0;
	}

private:
	/**
	 * @brief An outgoing/incoming edge: head/tail node and its visiting
	 *        frequency (number of alignment traversals).
	 */
	struct Edge {
		node_id neighbor;   ///< Other end of the edge.
		std::size_t weight; ///< Number of path traversals of this edge.
	};

	/**
	 * @brief One step of a sequence-to-graph alignment.
	 */
	struct AlignmentOp {
		enum class Type : unsigned char { kInsert, kDelete, kConsume };

		Type type;
		node_id node; ///< Graph node referenced by kDelete/kConsume.
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

	/// Base character of every node.
	std::vector<char> bases_;
	/// Outgoing edges per node.
	std::vector<std::vector<Edge>> out_edges_;
	/// Incoming edges per node.
	std::vector<std::vector<Edge>> in_edges_;
	/// Per dictionary sequence: node ids visited by that sequence.
	std::vector<std::vector<node_id>> paths_;
	/// Node ids in topological order of the DAG.
	std::vector<node_id> topo_;

	CostModel cost_;

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
		paths_.push_back(std::move(seed));
		rebuild_topo();

		for (std::size_t id = 1; id < dict.size(); ++id) {
			const Alignment al = align_to_graph(dict[id]);
			add_alignment(al.ops);
			rebuild_topo();
		}
	}

	/** Creates a new node with base @p c at the end of the graph. */
	node_id append_node(char c) {
		const node_id u = bases_.size();
		bases_.push_back(c);
		out_edges_.emplace_back();
		in_edges_.emplace_back();
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
	 * @brief Increments the visiting frequency of edge (u, v), creating it
	 *        if missing. Edges from the virtual start node are implicit and
	 *        therefore not stored.
	 */
	void add_edge(node_id u, node_id v) {
		if (u == kStart)
			return;
		if (Edge *e = find_edge(out_edges_[u], v)) {
			e->weight++;
			find_edge(in_edges_[v], u)->weight++;
		} else {
			out_edges_[u].push_back({v, 1});
			in_edges_[v].push_back({u, 1});
		}
	}

	/**
	 * @brief Recomputes topo_ (Kahn's algorithm), deterministically.
	 *
	 * Alignment-merged graphs are not topologically ordered by node id alone
	 * (mid-path insertions create edges like u_low -> u_new -> u_mid), so a
	 * real order is maintained for DP-based alignments and for clients.
	 */
	void rebuild_topo() {
		const std::size_t V = bases_.size();
		std::vector<std::size_t> indeg(V, 0);
		for (node_id u = 0; u < V; ++u)
			for (const Edge &e : out_edges_[u])
				indeg[e.neighbor]++;
		topo_.clear();
		topo_.reserve(V);
		for (node_id u = 0; u < V; ++u) // stable source order
			if (indeg[u] == 0)
				topo_.push_back(u);
		for (std::size_t h = 0; h < topo_.size(); ++h) {
			const node_id u = topo_[h];
			for (const Edge &e : out_edges_[u]) {
				std::size_t &d = indeg[e.neighbor];
				if (--d == 0)
					topo_.push_back(e.neighbor);
			}
		}
	}

	/**
	 * @brief Aligns a sequence to the current graph.
	 *
	 * Minimum-cost global alignment of @p s against the DAG: DP over
	 * (topological position, prefix length) with a virtual start row, using
	 * insertion/deletion/substitution costs from cost_. Cost ties prefer
	 * exact consumes, then deletions, then insertions.
	 *
	 * @param s Sequence to align.
	 * @return Alignment record: cost and alignment steps, to be merged with
	 *         add_alignment().
	 */
	Alignment align_to_graph(const std::string &s) {
		const std::size_t V = bases_.size();
		const std::size_t m = s.size();
		constexpr int kInf = std::numeric_limits<int>::max() / 2;
		const std::size_t R = V + 1; // DP rows; last row = virtual start.
		const std::size_t C = m + 1; // DP columns.
		const node_id vr = static_cast<node_id>(V);
		const node_id kInsCode =
		    static_cast<unsigned char>(AlignmentOp::Type::kInsert) + 1;
		const node_id kDelCode =
		    static_cast<unsigned char>(AlignmentOp::Type::kDelete) + 1;
		const node_id kConCode =
		    static_cast<unsigned char>(AlignmentOp::Type::kConsume) + 1;

		std::vector<int> dp(R * C, kInf);
		std::vector<unsigned char> op(R * C, 0);
		std::vector<node_id> from(R * C, kStart);
		std::vector<node_id> row_of(V, 0); // node id -> topological row.
		for (std::size_t p = 0; p < V; ++p)
			row_of[topo_[p]] = p;

		// Start row: dp[vr][j] = j * ins (insert-only prefix).
		dp[vr * C] = 0;
		for (std::size_t j = 1; j <= m; ++j) {
			const std::size_t k = vr * C + j;
			dp[k] = dp[k - 1] + cost_.ins;
			op[k] = kInsCode;
			from[k] = vr;
		}
		// dp[u][0]: delete graph bases between the start and node u.
		for (std::size_t p = 0; p < V; ++p) {
			int d = kInf;
			if (in_edges_[topo_[p]].empty()) {
				d = cost_.del;
			} else {
				for (const Edge &e : in_edges_[topo_[p]]) {
					const int c = dp[row_of[e.neighbor] * C] + cost_.del;
					if (c < d)
						d = c;
				}
			}
			dp[p * C] = d;
		}

		// Main DP over (topological row, prefix length).
		for (std::size_t p = 0; p < V; ++p) {
			const node_id u = topo_[p];
			for (std::size_t i = 1; i <= m; ++i) {
				const std::size_t k = p * C + i;
				int best = dp[k - 1] + cost_.ins; // insert s[i-1] at u
				unsigned char bop = kInsCode;
				node_id bfrom = p;
				int bpri = 0;

				auto consider_del = [&](node_id prow) {
					const int c = dp[prow * C + i] + cost_.del;
					if (c < best || (c == best && bpri < 1)) {
						best = c;
						bop = kDelCode;
						bfrom = prow;
						bpri = 1;
					}
				};
				auto consider_con = [&](node_id prow) {
					const int c =
					    dp[prow * C + i - 1] + cost_.match(bases_[u], s[i - 1]);
					int pri = 2;
					if (bases_[u] == s[i - 1] &&
					    cost_.match(bases_[u], s[i - 1]) == 0)
						pri = 3;
					if (c < best || (c == best && bpri < pri)) {
						best = c;
						bop = kConCode;
						bfrom = prow;
						bpri = pri;
					}
				};

				if (in_edges_[u].empty()) {
					consider_del(vr);
					consider_con(vr);
				} else {
					for (const Edge &e : in_edges_[u]) {
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
			const node_id u = topo_[p];
			int pri = 0;
			switch (op[p * C + m]) {
			case 2:
				pri = 1;
				break;
			case 3:
				pri = bases_[u] == s[m - 1] ? 3 : 2;
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
			const node_id row2node = (r == vr) ? kStart : topo_[r];
			switch (op[k]) {
			case 1:
				ops.push_back({AlignmentOp::Type::kInsert, kStart, s[i - 1]});
				--i;
				break;
			case 2:
				ops.push_back({AlignmentOp::Type::kDelete, row2node, '\0'});
				r = from[k];
				break;
			default:
				ops.push_back({AlignmentOp::Type::kConsume, row2node, s[i - 1]});
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
	 * Exact-match consumption reuses nodes and increments edge weights;
	 * substitutions and insertions create new nodes. Deletions are skipped
	 * and do not touch weights. Appends the resulting path to paths_ and
	 * refreshes the topological order.
	 */
	void add_alignment(const std::vector<AlignmentOp> &ops) {
		std::vector<node_id> path;
		path.reserve(ops.size());
		node_id prev = kStart;
		for (const AlignmentOp &a : ops) {
			switch (a.type) {
			case AlignmentOp::Type::kDelete:
				break; // the new path skips this graph node entirely
			case AlignmentOp::Type::kInsert: {
				const node_id u = append_node(a.ch);
				add_edge(prev, u);
				path.push_back(u);
				prev = u;
				break;
			}
			case AlignmentOp::Type::kConsume: {
				const node_id u = (bases_[a.node] == a.ch) ? a.node
				                                           : append_node(a.ch);
				add_edge(prev, u);
				path.push_back(u);
				prev = u;
				break;
			}
			}
		}
		paths_.push_back(std::move(path));
		rebuild_topo();
	}
};
