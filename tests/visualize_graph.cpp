/**
 * @file visualize_graph.cpp
 * @brief Renders POA graphs as Graphviz DOT for testcase files.
 *
 * For each testcase, builds the POA graph from its dictionary and writes one
 * `digraph` to stdout (multiple testcases yield multiple graphs, which `dot`
 * accepts in a single input file).
 *
 * Usage: visualize_graph [--no-compact] [--weights] [--heavy-light]
 *                        [--pos-ranges] [--path N]... <testcase-file-or-dir>...
 *
 *   --no-compact   skip merging of single-in/single-out node runs
 *   --weights      label each edge with its traversal count
 *   --heavy-light  heavy edges red+bold, light edges gray
 *   --pos-ranges   append each node's [j_min, j_max] to its label
 *   --path N       highlight dictionary sequence N's path (repeatable)
 *
 * Render with e.g.:
 *   visualize_graph --weights --heavy-light --pos-ranges --path 1 \
 *       tests/testcases/manual/01_basic.txt | dot -Tsvg -o graph.svg
 */

#include <hlp_grep/hlp_grep.hpp>

#include "testcase.hpp"

#include <cctype>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace hlp_grep;

namespace {

constexpr const char *EDGE_COLORS[] = {"blue",      "darkgreen", "purple",
                                       "darkorange", "brown",    "deeppink"};
constexpr const char *NODE_FILLS[] = {"lightblue", "palegreen", "thistle",
                                      "moccasin",  "mistyrose", "lightcyan"};
constexpr std::size_t NUM_COLORS =
    sizeof(EDGE_COLORS) / sizeof(EDGE_COLORS[0]);
/// pos_min/pos_max value of nodes on no stored path.
constexpr std::size_t NO_PATH = static_cast<std::size_t>(-1);

struct Options {
	bool compact = true;
	bool weights = false;
	bool heavy_light = false;
	bool pos_ranges = false;
	std::vector<std::size_t> paths;
};

/** Escapes a label for a double-quoted DOT string. */
std::string dot_escape(const std::string &s) {
	std::string out;
	for (const char c : s) {
		if (c == '"' || c == '\\')
			out += '\\';
		out += c;
	}
	return out;
}

void emit_graph(const fs::path &file, const Testcase &tc,
                const Options &opt) {
	const POAGraph graph(tc.dict, tc.cost, opt.compact);
	for (const std::size_t seq : opt.paths)
		if (seq >= graph.num_sequences()) {
			std::cerr << "--path " << seq << " out of range for " << file
			          << " (" << graph.num_sequences() << " sequences)\n";
			std::exit(1);
		}

	// Edge traversal counts (edges store no weights; count the paths).
	std::map<std::pair<POAGraph::node_id, POAGraph::node_id>, std::size_t>
	    weight;
	for (std::size_t s = 0; s < graph.num_sequences(); ++s) {
		const auto &path = graph.path(s);
		for (std::size_t i = 1; i < path.size(); ++i)
			weight[{path[i - 1], path[i]}]++;
	}

	// Highlight membership: first --path flag wins on overlaps.
	std::map<std::pair<POAGraph::node_id, POAGraph::node_id>, std::size_t>
	    edge_path;
	std::map<POAGraph::node_id, std::size_t> node_path;
	for (std::size_t k = 0; k < opt.paths.size(); ++k) {
		const auto &path = graph.path(opt.paths[k]);
		for (const POAGraph::node_id u : path)
			node_path.emplace(u, k);
		for (std::size_t i = 1; i < path.size(); ++i)
			edge_path.emplace(std::make_pair(path[i - 1], path[i]), k);
	}

	std::string name = file.stem().string();
	for (char &c : name)
		if (!std::isalnum(static_cast<unsigned char>(c)))
			c = '_';

	std::cout << "// " << file << "\ndigraph \"" << name << "\" {\n"
	          << "  rankdir=LR;\n";

	for (POAGraph::node_id u = 0; u < graph.num_nodes(); ++u) {
		std::string label =
		    std::to_string(u) + ":" + dot_escape(graph.seq(u));
		if (opt.pos_ranges) {
			const auto [lo, hi] = graph.pos_range(u);
			label += "\\n";
			label += (lo == NO_PATH) ? "[-,-]"
			                         : "[" + std::to_string(lo) + "," +
			                               std::to_string(hi) + "]";
		}
		std::cout << "  " << u << " [label=\"" << label << "\"";
		if (const auto it = node_path.find(u); it != node_path.end())
			std::cout << " style=filled fillcolor=\""
			          << NODE_FILLS[it->second % NUM_COLORS] << "\"";
		std::cout << "];\n";
	}

	for (const auto &[ends, w] : weight) {
		const auto [u, v] = ends;
		const bool heavy =
		    graph.edge_type(u, v) == POAGraph::EdgeType::HEAVY;
		std::cout << "  " << u << " -> " << v;
		std::string attrs;
		if (opt.weights)
			attrs += " label=\"" + std::to_string(w) + "\"";
		// Path highlights replace the heavy/light coloring; the
		// heavy/light distinction survives as solid vs dashed.
		if (const auto it = edge_path.find(ends); it != edge_path.end())
			attrs += std::string(" color=") +
			         EDGE_COLORS[it->second % NUM_COLORS] +
			         " penwidth=3.0" +
			         (heavy ? " style=solid" : " style=dashed");
		else if (opt.heavy_light)
			attrs += heavy ? " color=red penwidth=2.0"
			               : " color=gray50";
		if (!attrs.empty())
			std::cout << " [" << attrs.substr(1) << "]";
		std::cout << ";\n";
	}
	std::cout << "}\n";
}

} // namespace

int main(int argc, char **argv) {
	Options opt;
	std::vector<const char *> args;
	for (int i = 1; i < argc; ++i) {
		const std::string flag = argv[i];
		if (flag == "--no-compact") {
			opt.compact = false;
		} else if (flag == "--weights") {
			opt.weights = true;
		} else if (flag == "--heavy-light") {
			opt.heavy_light = true;
		} else if (flag == "--pos-ranges") {
			opt.pos_ranges = true;
		} else if (flag == "--path" && i + 1 < argc) {
			try {
				opt.paths.push_back(
				    std::stoul(argv[++i]));
			} catch (...) {
				std::cerr << "--path needs an integer sequence id\n";
				return 1;
			}
		} else if (!flag.empty() && flag[0] == '-') {
			std::cerr << "unknown flag: " << flag << '\n';
			return 1;
		} else {
			args.push_back(argv[i]);
		}
	}
	if (args.empty()) {
		std::cerr << "usage: " << argv[0]
		          << " [--no-compact] [--weights] [--heavy-light]"
		             " [--pos-ranges] [--path N]... <testcase-file-or-dir>...\n";
		return 1;
	}

	// collect_files() indexes argv terms from 1 (argv[0] is the program
	// name), so keep a dummy first element.
	auto arg_ptrs = [&] {
		std::vector<char *> p{argv[0]};
		for (const char *a : args)
			p.push_back(const_cast<char *>(a));
		return p;
	}();
	const auto files = collect_files(static_cast<int>(arg_ptrs.size()),
	                                 arg_ptrs.data());
	if (files.empty()) {
		std::cout << "no testcase files found\n";
		return 0;
	}
	for (const auto &file : files)
		emit_graph(file, parse_testcase(file), opt);
	return 0;
}
