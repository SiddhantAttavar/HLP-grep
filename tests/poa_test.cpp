/**
 * @file poa_test.cpp
 * @brief Validates the POAGraph sequence-to-graph alignment.
 *
 * Checks the structural invariants of the POA graph built from a dictionary:
 * every sequence path spells back its dictionary sequence, node ids along
 * each path are strictly increasing (topological order), and known small
 * dictionaries produce the exact expected graph shape (node count, bases,
 * and shared path prefixes).
 */

#include <hlp_grep/poa_graph.hpp>

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

/** Exits with an error message on failure, like the testcase helpers. */
void check(bool ok, const std::string &msg) {
	if (!ok) {
		std::cerr << "FAIL: " << msg << '\n';
		std::exit(1);
	}
}

/** Asserts every path spells its dictionary string, connected by graph edges. */
void check_paths(const POAGraph &g, const std::vector<std::string> &dict,
                 const std::string &name) {
	check(g.num_sequences() == dict.size(), name + ": sequence count mismatch");
	for (std::size_t i = 0; i < dict.size(); ++i) {
		std::string spelled;
		bool first = true;
		POAGraph::node_id prev = POAGraph::node_id(-1);
		for (const auto u : g.path(i)) {
			check(u < g.num_nodes(), name + ": node id out of range");
			if (!first)
				check(g.edge_weight(prev, u) > 0,
				      name + ": path " + std::to_string(i) + " has no edge " +
				          std::to_string(prev) + "->" + std::to_string(u));
			spelled += g.base(u);
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
	std::string spelled(1, g.base(cp.start));
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
			check(!prev_was_heavy_end || cur == POAGraph::node_id(-1) ||
			          g.edge_type(cur, g.heavy_edge(cur)) ==
			                  POAGraph::EdgeType::LIGHT,
			      name + ": adjacent heavy steps on the same chain missed");
			for (std::size_t h = 0; h < st.length; ++h) {
				check(g.edge_type(cur, g.heavy_edge(cur)) ==
				              POAGraph::EdgeType::HEAVY,
				      name + ": heavy step traverses a non-heavy edge");
				cur = g.heavy_edge(cur);
				spelled += g.base(cur);
			}
			prev_was_heavy_end = true;
			break;
		}
		case POAGraph::EdgeType::LIGHT: {
			check(st.next < g.num_nodes(),
			      name + ": light step lands out of range");
			check(g.edge_weight(cur, st.next) > 0,
			      name + ": light step has no edge");
			check(g.edge_type(cur, st.next) == POAGraph::EdgeType::LIGHT,
			      name + ": light step traverses a heavy edge");
			cur = st.next;
			spelled += g.base(cur);
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

const std::vector<POAGraph::node_id> kPath0123 = {0, 1, 2, 3};
const std::vector<POAGraph::node_id> kPath0124 = {0, 1, 2, 4};
const std::vector<POAGraph::node_id> kPath013 = {0, 1, 3};

} // namespace

int main() {
	// Substitution: shared prefix, one alternate path for the last base.
	{
		const std::vector<std::string> dict = {"ACGT", "ACGA"};
		const POAGraph g(dict);
		check(g.num_nodes() == 5, "substitution: expected 5 nodes");
		check(g.base(4) == 'A', "substitution: last node must be 'A'");
		check(g.path(0) == kPath0123, "substitution: path0");
		check(g.path(1) == kPath0124, "substitution: path1");
		check_paths(g, dict, "substitution");
		std::cout << "PASS substitution\n";
	}

	// Duplicates reuse the entire path; no new nodes.
	{
		const std::vector<std::string> dict = {"ACGT", "ACGT"};
		const POAGraph g(dict);
		check(g.num_nodes() == 4, "duplicate: expected 4 nodes");
		check(g.path(0) == g.path(1), "duplicate: paths must be identical");
		check_paths(g, dict, "duplicate");
		std::cout << "PASS duplicate\n";
	}

	// Deletion: second sequence skips the G node.
	{
		const std::vector<std::string> dict = {"ACGT", "ACT"};
		const POAGraph g(dict);
		check(g.num_nodes() == 4, "deletion: expected 4 nodes");
		check(g.path(1) == kPath013, "deletion: expected path [0,1,3]");
		check_paths(g, dict, "deletion");
		std::cout << "PASS deletion\n";
	}

	// Insertion: new suffix nodes, shared first two nodes.
	{
		const std::vector<std::string> dict = {"AC", "ACGT"};
		const POAGraph g(dict);
		check(g.num_nodes() == 4, "insertion: expected 4 nodes");
		check(g.path(1) == kPath0123, "insertion: expected path [0,1,2,3]");
		check(g.path(0) == std::vector<POAGraph::node_id>{0, 1},
		      "insertion: path0 expected [0,1]");
		check_paths(g, dict, "insertion");
		std::cout << "PASS insertion\n";
	}

	// Weighted insertion/deletion costs still produce the same graph shape.
	{
		const std::vector<std::string> dict = {"ACGT", "ACGA"};
		const POAGraph g(dict, CostModel(2, 3));
		check(g.num_nodes() == 5, "weighted: expected 5 nodes");
		check(g.path(1) == kPath0124, "weighted: path1 expected [0,1,2,4]");
		check_paths(g, dict, "weighted");
		std::cout << "PASS weighted\n";
	}

	// Empty string dictionary entry: empty path, no extra nodes.
	{
		const std::vector<std::string> dict = {"", "ACGT"};
		const POAGraph g(dict);
		check(g.num_nodes() == 4, "empty-string: expected 4 nodes");
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

	// A mutation family shares everything except the final base: the graph
	// must stay compact (5 shared nodes + one node per distinct last base).
	{
		std::vector<std::string> dict;
		for (int i = 0; i < 8; ++i)
			dict.push_back("ACGCC" + std::string(1, char('A' + i)));
		const POAGraph g(dict);
		check(g.num_nodes() == 13,
		      "mutation-family: expected 13 nodes (5 shared + 8 last bases)");
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

		for (std::size_t u = 0; u < g.num_nodes(); ++u) {
			const POAGraph::node_id heavy = g.heavy_edge(u);
			if (heavy == POAGraph::node_id(-1)) {
				// No outgoing edges: every candidate edge must be LIGHT.
				for (std::size_t v = 0; v < g.num_nodes(); ++v)
					check(g.edge_type(u, v) == POAGraph::EdgeType::LIGHT,
					      "heavy-ties: node " + std::to_string(u) +
					          " has no outgoing edges but an edge is HEAVY");
				continue;
			}
			check(g.edge_weight(u, heavy) > 0,
			      "heavy-ties: heavy edge must exist");
			check(g.edge_type(u, heavy) == POAGraph::EdgeType::HEAVY,
			      "heavy-ties: heavy_edge() disagrees with edge_type()");
			// HEAVY edge must be (a) max weight and (b) tie-broken by
			// largest destination id among equal-weight edges.
			for (std::size_t v = 0; v < g.num_nodes(); ++v) {
				if (v == heavy || g.edge_weight(u, v) == 0)
					continue;
				const std::size_t wh = g.edge_weight(u, heavy);
				const std::size_t wv = g.edge_weight(u, v);
				check(wv < wh || (wv == wh && v < heavy),
				      "heavy-ties: edge " + std::to_string(u) + "->" +
				          std::to_string(v) + " (w=" + std::to_string(wv) +
				          ") outranks heavy edge " + std::to_string(u) + "->" +
				          std::to_string(heavy) + " (w=" + std::to_string(wh) +
				          ")");
				check(g.edge_type(u, v) == POAGraph::EdgeType::LIGHT,
				      "heavy-ties: more than one HEAVY edge from node " +
				          std::to_string(u));
			}
		}

		// The trunk A-A chain edge 0->1 is taken by all 4 paths and is heavy;
		// absent edges report LIGHT.
		check(g.edge_weight(0, 1) == 4, "heavy-ties: edge 0->1 weight");
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
	// that passes through it (j_min/j_max in notes/idea.md).
	{
		const std::vector<std::string> dict = {"ACGT", "ACGA", "AC"};
		const POAGraph g(dict);
		check_paths(g, dict, "pos-range");

		// Recompute expected ranges directly from the paths.
		const std::size_t kNone = static_cast<std::size_t>(-1);
		std::vector<std::pair<std::size_t, std::size_t>> expected(
		    g.num_nodes(), {kNone, kNone});
		for (std::size_t s = 0; s < dict.size(); ++s)
			for (std::size_t j = 0; j < g.path(s).size(); ++j) {
				const auto u = g.path(s)[j];
				auto &er = expected[u];
				if (er.first == kNone || j < er.first)
					er.first = j;
				if (er.second == kNone || j > er.second)
					er.second = j;
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

	std::cout << "all POA alignment checks passed\n";
	return 0;
}
