/**
 * @file poa_test.cpp
 * @brief Validates the POAGraph sequence-to-graph alignment.
 *
 * Checks the structural invariants of the POA graph built from a dictionary:
 * every sequence path spells back its dictionary sequence over connected
 * graph edges, and known small dictionaries produce the exact expected graph
 * shape (node count, bases, and shared path prefixes).
 */

#include <hlp_grep/poa_graph.hpp>

#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>


using namespace hlp_grep;

namespace {

/** Exits with an error message on failure, like the testcase helpers. */
void check(bool ok, const std::string &msg) {
	if (!ok) {
		std::cerr << "FAIL: " << msg << '\n';
		std::exit(1);
	}
}

/**
 * @brief Builds the (u, v) -> traversal count map over all stored paths.
 *
 * Replaces the removed POAGraph::edge_weight(): the graph no longer
 * exposes on-demand weight queries, but the test still verifies
 * heavy/light classification against the raw path data.
 */
std::unordered_map<unsigned long long, std::size_t>
path_weights(const POAGraph &g) {
	std::unordered_map<unsigned long long, std::size_t> weights;
	for (std::size_t s = 0; s < g.num_sequences(); ++s) {
		const auto &p = g.path(s);
		for (std::size_t i = 1; i < p.size(); ++i)
			++weights[static_cast<unsigned long long>(p[i - 1]) << 32 |
			          static_cast<unsigned long long>(p[i])];
	}
	return weights;
}

/** Weight of edge (u, v) from a path_weights() map; 0 if never traversed. */
std::size_t weight_of(
    const std::unordered_map<unsigned long long, std::size_t> &weights,
    POAGraph::node_id u, POAGraph::node_id v) {
	const auto it = weights.find(static_cast<unsigned long long>(u) << 32 |
	                             static_cast<unsigned long long>(v));
	return it == weights.end() ? 0 : it->second;
}

/** Asserts every path spells its dictionary string, connected by graph edges. */
void check_paths(const POAGraph &g, const std::vector<std::string> &dict,
                 const std::string &name) {
	check(g.num_sequences() == dict.size(), name + ": sequence count mismatch");
	const auto weights = path_weights(g);
	for (std::size_t i = 0; i < dict.size(); ++i) {
		std::string spelled;
		bool first = true;
		POAGraph::node_id prev = POAGraph::node_id(-1);
		for (const auto u : g.path(i)) {
			check(u < g.num_nodes(), name + ": node id out of range");
			if (!first)
				check(weight_of(weights, prev, u) > 0,
				      name + ": path " + std::to_string(i) + " has no edge " +
				          std::to_string(prev) + "->" + std::to_string(u));
			spelled += g.seq(u);
			prev = u;
			first = false;
		}
		check(spelled == dict[i], name + ": path " + std::to_string(i) +
		                              " spells '" + spelled + "', expected '" +
		                              dict[i] + "'");
	}
}

/**
 * @brief Validates the compressed path of one sequence: replaying the steps
 *        from the stored start node must spell the dictionary sequence, and
 *        every step must be well-formed and match the graph edge types.
 */
void check_compressed(const POAGraph &g, std::size_t seq,
                      const std::vector<std::string> &dict,
                      const std::string &name) {
	const auto p = g.path(seq);
	const auto cp = g.compressed_path(seq);
	check(!p.empty() && cp.start == p.front(),
	      name + ": start node must be path front");
	std::string spelled = g.seq(cp.start);
	POAGraph::node_id cur = cp.start;
	bool prev_was_heavy_end = false;
	for (const auto &st : cp.steps) {
		switch (st.type) {
		case POAGraph::EdgeType::HEAVY: {
			check(st.length > 0,
			      name + ": heavy step must advance > 0");
			// A HEAVY step must not occur right after another step whose
			// heavy chain could have continued: no two steps may share a
			// heavy chain (otherwise compression missed a merge).
			const auto chain_next = g.node(cur).heavy_neighbour;
			check(!prev_was_heavy_end || !chain_next ||
			          g.edge_type(cur, *chain_next) ==
			              POAGraph::EdgeType::LIGHT,
			      name + ": adjacent heavy steps on the same chain missed");
			for (int h = 0; h < st.length; ++h) {
				const auto hn = g.node(cur).heavy_neighbour;
				check(hn.has_value() &&
				              g.edge_type(cur, *hn) ==
				                  POAGraph::EdgeType::HEAVY,
				      name + ": heavy step traverses a non-heavy edge");
				cur = *hn;
				spelled += g.seq(cur);
			}
			prev_was_heavy_end = true;
			break;
		}
		case POAGraph::EdgeType::LIGHT: {
			check(st.next < g.num_nodes(),
			      name + ": light step lands out of range");
			check(g.edge_type(cur, st.next) == POAGraph::EdgeType::LIGHT,
			      name + ": light step traverses a heavy edge");
			cur = st.next;
			spelled += g.seq(cur);
			prev_was_heavy_end = false;
			break;
		}
		}
	}
	check(spelled == dict[seq], name + ": replay of path " +
	                                std::to_string(seq) + " spells '" +
	                                spelled + "', expected '" + dict[seq] +
	                                "'");
}

const std::vector<POAGraph::node_id> PATH_AC_T = {0, 1};   // ACG + T
const std::vector<POAGraph::node_id> PATH_AC_A = {0, 2};   // ACG + A
const std::vector<POAGraph::node_id> PATH_AC_G_T = {0, 1, 2}; // AC + G + T
const std::vector<POAGraph::node_id> PATH_AC_T2 = {0, 2};  // AC + T
const std::vector<POAGraph::node_id> PATH_0123 = {0, 1, 2, 3};

} // namespace

int main() {
	// Substitution: shared prefix, one alternate path for the last base.
	// The shared prefix ACG compacts into one node; the two last bases
	// branch off it.
	{
		const std::vector<std::string> dict = {"ACGT", "ACGA"};
		const POAGraph g(dict);
		check(g.num_nodes() == 3, "substitution: expected 3 nodes");
		check(g.seq(1) == "T", "substitution: node 1 must be 'T'");
		check(g.seq(2) == "A", "substitution: node 2 must be 'A'");
		check(g.path(0) == PATH_AC_T, "substitution: path0");
		check(g.path(1) == PATH_AC_A, "substitution: path1");
		check_paths(g, dict, "substitution");
		std::cout << "PASS substitution\n";
	}

	// Duplicates reuse the entire path; no new nodes. The trunk ACG
	// compacts; the shared last base T stays alone (it ends both paths).
	{
		const std::vector<std::string> dict = {"ACGT", "ACGT"};
		const POAGraph g(dict);
		check(g.num_nodes() == 2, "duplicate: expected 2 nodes");
		check(g.path(0) == g.path(1), "duplicate: paths must be identical");
		check_paths(g, dict, "duplicate");
		std::cout << "PASS duplicate\n";
	}

	// Deletion: second sequence skips the G node. A and C merge into
	// one node (C's two out-edges keep the skip edge intact); G and T
	// cannot merge (T has two in-edges).
	{
		const std::vector<std::string> dict = {"ACGT", "ACT"};
		const POAGraph g(dict);
		check(g.num_nodes() == 3, "deletion: expected 3 nodes");
		check(g.path(1) == PATH_AC_T2, "deletion: expected path [0,2]");
		check(g.path(0) == PATH_AC_G_T, "deletion: expected path [0,1,2]");
		check_paths(g, dict, "deletion");
		std::cout << "PASS deletion\n";
	}

	// Insertion: new suffix nodes, shared first two nodes. Nothing
	// merges: a path ends on C (blocking A->C and the C-G-T run).
	{
		const std::vector<std::string> dict = {"AC", "ACGT"};
		const POAGraph g(dict);
		check(g.num_nodes() == 4, "insertion: expected 4 nodes");
		check(g.path(1) == PATH_0123, "insertion: expected path [0,1,2,3]");
		check(g.path(0) == std::vector<POAGraph::node_id>{0, 1},
		      "insertion: path0 expected [0,1]");
		check(g.seq(1) == "C", "insertion: node 1 must be 'C'");
		check_paths(g, dict, "insertion");
		std::cout << "PASS insertion\n";
	}

	// Weighted insertion/deletion costs still produce the same graph shape.
	{
		const std::vector<std::string> dict = {"ACGT", "ACGA"};
		const CostModel weighted(2, 3); // reference must outlive graph
		const POAGraph g(dict, weighted);
		check(g.num_nodes() == 3, "weighted: expected 3 nodes");
		check(g.path(1) == PATH_AC_A, "weighted: path1 expected [0,2]");
		check_paths(g, dict, "weighted");
		std::cout << "PASS weighted\n";
	}

	// Empty string dictionary entry: empty path, no extra nodes.
	{
		const std::vector<std::string> dict = {"", "ACGT"};
		const POAGraph g(dict);
		check(g.num_nodes() == 2, "empty-string: expected 2 nodes");
		check(g.path(0).empty(), "empty-string: path 0 must be empty");
		check_paths(g, dict, "empty-string");
		std::cout << "PASS empty-string\n";
	}

	// Empty dictionary: no nodes, no crash.
	{
		const POAGraph g(std::vector<std::string>{});
		check(g.num_nodes() == 0, "empty-dict: expected 0 nodes");
		check(g.num_sequences() == 0, "empty-dict: expected 0 sequences");
		std::cout << "PASS empty-dict\n";
	}

	// A mutation family shares everything except the final base: the
	// shared trunk compacts into one node per in-tact run (the two
	// trailing Cs gain extra in-edges from a mid-trunk insertion, so
	// they stay separate), plus one node per distinct last base.
	{
		std::vector<std::string> dict;
		for (int i = 0; i < 8; ++i)
			dict.push_back("ACGCC" + std::string(1, char('A' + i)));
		const POAGraph g(dict);
		check(g.num_nodes() == 11, "mutation-family: expected 11 nodes");
		check(g.seq(0) == "ACG", "mutation-family: trunk label");
		check_paths(g, dict, "mutation-family");
		std::cout << "PASS mutation-family\n";
	}

	// Heavy/light classification: exactly one HEAVY outgoing edge per node
	// with outgoing edges, chosen as the max-weight edge, ties broken by the
	// largest destination node id.
	{
		// Dictionary designed so each node has multiple, tie-heavy outgoing
		// edges: AACT and AAAT share "AA", then both branch to C/T, and
		// C-sequence continues with CT while the other ends; the last split
		// (C -> T) has equal weight to C -> ... from its own path.
		const std::vector<std::string> dict = {"AACT", "AACT", "AAAT", "AACTA"};
		const POAGraph g(dict);
		check_paths(g, dict, "heavy-ties");
		const auto weights = path_weights(g);

		for (std::size_t u = 0; u < g.num_nodes(); ++u) {
			const auto heavy = g.node(u).heavy_neighbour;
			if (!heavy.has_value()) {
				// No outgoing edges: every candidate edge must be LIGHT.
				for (std::size_t v = 0; v < g.num_nodes(); ++v)
					check(g.edge_type(u, v) == POAGraph::EdgeType::LIGHT,
					      "heavy-ties: node " + std::to_string(u) +
					          " has no outgoing edges but an edge is HEAVY");
				continue;
			}
			check(weight_of(weights, u, *heavy) > 0,
			      "heavy-ties: heavy edge must exist");
			check(g.edge_type(u, *heavy) == POAGraph::EdgeType::HEAVY,
			      "heavy-ties: heavy_neighbour disagrees with edge_type()");
			// HEAVY edge must be (a) max weight and (b) tie-broken by
			// largest destination id among equal-weight edges.
			for (std::size_t v = 0; v < g.num_nodes(); ++v) {
				if (v == *heavy || weight_of(weights, u, v) == 0)
					continue;
				const std::size_t wh = weight_of(weights, u, *heavy);
				const std::size_t wv = weight_of(weights, u, v);
				check(wv < wh || (wv == wh && v < *heavy),
				      "heavy-ties: edge " + std::to_string(u) + "->" +
				          std::to_string(v) + " (w=" + std::to_string(wv) +
				          ") outranks heavy edge " + std::to_string(u) + "->" +
				          std::to_string(*heavy) + " (w=" + std::to_string(wh) +
				          ")");
				check(g.edge_type(u, v) == POAGraph::EdgeType::LIGHT,
				      "heavy-ties: more than one HEAVY edge from node " +
				          std::to_string(u));
			}
		}

		// The trunk AA chain merges into one node whose edge to C is
		// taken by 3 of the 4 paths and is heavy.
		check(weight_of(weights, 0, 1) == 3, "heavy-ties: edge 0->1 weight");
		check(g.edge_type(0, 1) == POAGraph::EdgeType::HEAVY,
		      "heavy-ties: edge 0->1 must be heavy");
		check(g.edge_type(3, 0) == POAGraph::EdgeType::LIGHT,
		      "heavy-ties: absent edge must be LIGHT");
		std::cout << "PASS heavy-ties\n";
	}

	// Compressed path: every sequence's compressed representation must
	// reproduce the original node path when replayed (HEAVY steps follow the
	// unique heavy chain; LIGHT steps land explicitly on `next`).
	{
		const std::vector<std::string> dict = {"AACT", "AACT", "AAAT",
		                                       "AACTA", "ACGCCA", "ACGCCC"};
		const POAGraph g(dict);
		check_paths(g, dict, "compressed-path");

		// Steps count: every path edge with n edges compresses into at most
		// n steps, and every step covers at least one edge.
		for (std::size_t s = 0; s < dict.size(); ++s) {
			const auto &p = g.path(s);
			check(g.compressed_path(s).steps.size() <=
			              p.size() - (p.empty() ? 0 : 1),
			      "compressed-path: too many steps");
		}
		for (std::size_t s = 0; s < dict.size(); ++s)
			check_compressed(g, s, dict, "compressed-path");

		// Extended coverage: sequences with backtracking substitutions and
		// varied shared prefixes to shake out heavy/light misclassification
		// in the compressed representation.
		{
			const std::vector<std::string> dict2 = {
			    "ACGTACGT", "ACGTTCGT", "ACGTACGTA", "ACG",   "ACGA",
			    "AAAA",     "AACGT",   "TTTT",     "TTTTA", "ACGTACGT"};
			const POAGraph g2(dict2);
			check_paths(g2, dict2, "compressed-path-extra");
			for (std::size_t s = 0; s < dict2.size(); ++s)
				check_compressed(g2, s, dict2, "compressed-path-extra");
		}

		std::cout << "PASS compressed-path\n";
	}

	// Position ranges: per node, min/max 0-based offsets across every path
	// that passes through it (j_min/j_max).
	{
		const std::vector<std::string> dict = {"ACGT", "ACGA", "AC"};
		const POAGraph g(dict);
		check_paths(g, dict, "pos-range");

		// Recompute expected ranges directly from the paths: an offset
		// counts the bases consumed before a node's label starts.
		const std::size_t NONE = static_cast<std::size_t>(-1);
		std::vector<std::pair<std::size_t, std::size_t>> expected(
		    g.num_nodes(), {NONE, NONE});
		for (std::size_t s = 0; s < dict.size(); ++s) {
			std::size_t offset = 0;
			for (const auto u : g.path(s)) {
				auto &er = expected[u];
				if (er.first == NONE || offset < er.first)
					er.first = offset;
				if (er.second == NONE || offset > er.second)
					er.second = offset;
				offset += g.seq(u).size();
			}
		}
		for (std::size_t u = 0; u < g.num_nodes(); ++u) {
			const auto [lo, hi] = g.pos_range(u);
			check(lo == expected[u].first && hi == expected[u].second,
			      "pos-range: node " + std::to_string(u) + " expected (" +
			          std::to_string(expected[u].first) + "," +
			          std::to_string(expected[u].second) + ") got (" +
			          std::to_string(lo) + "," + std::to_string(hi) + ")");
			check(lo <= hi, "pos-range: min > max");
			check(lo < g.num_nodes() ? hi >= lo : true,
			      "pos-range: sentinel inconsistency");
		}

		// The mid-graph nodes are in every path at every offset 0..3:
		// node 0 (A) occurs at offset 0 in all paths.
		{
			const auto [lo, hi] = g.pos_range(0);
			check(lo == 0 && hi == 0, "pos-range: node 0 expected (0,0)");
		}
		std::cout << "PASS pos-range\n";
	}

	// Compaction maximality: no contractible chain edge may survive
	// compaction. An edge (u, v) is contractible exactly when u's
	// out-degree is 1, v's in-degree is 1, and no stored path forbids
	// it: u must not end a path, and v must neither start nor end one
	// (paths may begin or end mid-run via deletions, and contraction
	// would force those across label characters they never had).
	{
		const std::vector<std::string> dict = {"ACGT", "ACT",  "ACGA",
		                                       "CGT",  "ACGTAC", "AC"};
		const POAGraph g(dict);
		check_paths(g, dict, "compaction-maximal");
		for (std::size_t s = 0; s < dict.size(); ++s)
			check_compressed(g, s, dict, "compaction-maximal");

		const auto weights = path_weights(g);
		std::unordered_map<POAGraph::node_id, std::size_t> outdeg, indeg;
		for (const auto &[edge, w] : weights) {
			(void)w;
			outdeg[static_cast<POAGraph::node_id>(edge >> 32)]++;
			indeg[static_cast<POAGraph::node_id>(edge & 0xffffffffUL)]++;
		}
		std::vector<char> is_start(g.num_nodes(), 0), is_end(g.num_nodes(), 0);
		for (std::size_t s = 0; s < dict.size(); ++s) {
			const auto &p = g.path(s);
			if (p.empty())
				continue;
			is_start[p.front()] = 1;
			is_end[p.back()] = 1;
		}
		for (const auto &[edge, w] : weights) {
			(void)w;
			const auto u = static_cast<POAGraph::node_id>(edge >> 32);
			const auto v =
			    static_cast<POAGraph::node_id>(edge & 0xffffffffUL);
			const bool contractible = outdeg[u] == 1 && indeg[v] == 1 &&
			                          !is_start[v] && !is_end[v] &&
			                          !is_end[u];
			check(!contractible,
			      "compaction-maximal: contractible edge " +
			          std::to_string(u) + "->" + std::to_string(v) +
			          " survived");
		}
		std::cout << "PASS compaction-maximal\n";
	}

	std::cout << "all POA alignment checks passed\n";
	return 0;
}
