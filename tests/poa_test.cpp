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

	std::cout << "all POA alignment checks passed\n";
	return 0;
}
