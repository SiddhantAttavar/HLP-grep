/**
 * @file bench.cpp
 * @brief Benchmarks a single testcase with one or more methods: build time,
 *        peak memory and per-query runtime for each method, written as one
 *        JSON file.
 *
 * Usage: bench --method hlp_grep [--method naive ...] <testcase-file>
 *              [--out out.json] [--index graph.hlpg]
 *
 * The hlp_grep method accepts a prebuilt index file (written by `save`)
 * via --index: the solver is reconstructed from it without building the
 * dictionary, and the reported build time is zero.
 *
 * Each measured method runs in its own forked child process; the parent
 * recovers the child's rusage (ru_maxrss = peak RSS of the whole run, parse
 * included) via wait4 and assembles the JSON. Per query, the child reports
 * the index, threshold k, runtime, match count and an FNV-1a hash of the
 * result vector, so runs of different methods can be cross-checked for
 * identical results.
 */

#include <hlp_grep/hlp_grep.hpp>

#include "naive_solver.hpp"
#include "testcase.hpp"

#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

// WFA2-lib headers are C-only (no extern "C" inside): wrap the inclusion.
#ifdef HLP_GREP_WITH_WFA2
extern "C" {
#include <wavefront/wavefront_align.h>
#include <wavefront/wavefront_attributes.h>
}
#endif

#ifdef HLP_GREP_WITH_DT_PATRICIA
#include <dt_patricia/dt_patricia.hpp>
#endif

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <exception>
#include <iomanip>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <cctype>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace hlp_grep;

namespace {

/** FNV-1a 64-bit hash of a results vector: bytes of each (id, dist) pair. */
std::uint64_t results_hash(const std::vector<Result> &results) {
	std::uint64_t h = 0xcbf29ce484222325ULL;
	for (const auto &r : results) {
		const auto feed = [&h](auto value) {
			unsigned char bytes[sizeof(value)];
			std::memcpy(bytes, &value, sizeof(value));
			for (const unsigned char b : bytes) {
				h ^= b;
				h *= 0x100000001b3ULL;
			}
		};
		feed(r.id);
		feed(r.dist);
	}
	return h;
}

double elapsed_ms(const std::chrono::steady_clock::time_point begin,
                  const std::chrono::steady_clock::time_point end) {
	return std::chrono::duration<double, std::milli>(end - begin).count();
}

/** Serialised per-query record sent from the child to the parent. */
struct QueryStat {
	std::size_t index;
	int k;
	double runtime_ms;
	std::uint64_t num_matches;
	std::uint64_t hash;
};

/** Collected measurements of one method run. */
struct MethodRun {
	std::string name;
	double build_ms;
	long peak_bytes;
	std::vector<QueryStat> queries;
};

void write_all(int fd, const std::string &s) {
	std::size_t done = 0;
	while (done < s.size()) {
		const long n = write(fd, s.data() + done, s.size() - done);
		if (n <= 0)
			::_exit(1);
		done += static_cast<std::size_t>(n);
	}
}

void report_query(int fd, std::size_t index, int k,
                  const std::vector<Result> &results,
                  const std::chrono::steady_clock::time_point begin,
                  const std::chrono::steady_clock::time_point end) {
	const double ms = elapsed_ms(begin, end);
	write_all(fd, "QUERY " + std::to_string(index) + " " +
	                  std::to_string(k) + " " + std::to_string(ms) + " " +
	                  std::to_string(results_hash(results)) + " " +
	                  std::to_string(results.size()) + "\n");
}

#ifdef HLP_GREP_WITH_DT_PATRICIA
/** Converts DT-Patricia results into library Result records. */
std::vector<Result> from_dt_patricia(
    const std::vector<dt_patricia::AlignmentResult> &rs) {
	std::vector<Result> out;
	out.reserve(rs.size());
	for (const auto &r : rs)
		out.push_back({static_cast<std::size_t>(r.string_id),
		               static_cast<int>(r.score)});
	return out;
}
#endif

/** Answers every testcase query with an hlp_grep solver and reports it. */
void report_hlp_queries(int fd, const Solver &solver, const Testcase &tc) {
	for (std::size_t i = 0; i < tc.queries.size(); ++i) {
		const auto [k, query] = tc.queries[i];
		const auto b = std::chrono::steady_clock::now();
		const auto results = solver.query(query, k);
		const auto e = std::chrono::steady_clock::now();
		report_query(fd, i, k, results, b, e);
	}
}

/** Loads a prebuilt index file or exits the child on failure. */
Solver load_hlp_index(const std::string &index_file, const Testcase &tc) {
	try {
		return Solver::load_graph(index_file, tc.cost);
	} catch (const std::exception &e) {
		std::cerr << "cannot load index file " << index_file << ": "
		          << e.what() << '\n';
		::_exit(1);
	}
}

/**
 * @brief In-memory B^ed-tree gram-count index (SIGMOD 2010, Zhang et al.).
 *
 * bed-tree (https://github.com/ZhangZhenjie/bed-tree) is a disk-based B+
 * tree executable without a reusable library API, so bench reimplements
 * its gram-count order range query natively: a q-gram inverted index over
 * the dictionary with exact count-filter pruning and exact verification
 * via edit_distance. A single edit destroys at most BED_Q q-grams, so a
 * dictionary sequence within distance k of the query shares at least
 * |grams(query)| - k * BED_Q distinct q-grams with it; anything below
 * that count (or failing the |len(query) - len(s)| <= k length filter)
 * cannot match. Unit-cost model only (validated in main).
 */
class BedTreeIndex {
public:
	explicit BedTreeIndex(const std::vector<std::string> &dict) {
		for (std::size_t i = 0; i < dict.size(); ++i) {
			if (dict[i].size() < BED_Q)
				continue; // too short to carry a q-gram: verified directly
			std::unordered_set<std::string> seen;
			for (std::size_t p = 0; p + BED_Q <= dict[i].size(); ++p)
				seen.insert(dict[i].substr(p, BED_Q));
			for (const auto &g : seen)
				postings[g].push_back(i); // ids ascending: already sorted
		}
	}

	std::vector<Result> query(const std::vector<std::string> &dict,
	                          const CostModel &cost, const std::string &alphabet,
	                          const std::string &query_str, int k) const {
		const auto [c_min, g_min] = model_mins(cost, alphabet);
		// Distinct query q-grams; an empty set (short query) disables the
		// count filter and every length-passing sequence is verified.
		std::unordered_set<std::string> qgrams;
		if (query_str.size() >= BED_Q)
			for (std::size_t p = 0; p + BED_Q <= query_str.size(); ++p)
				qgrams.insert(query_str.substr(p, BED_Q));
		const long need =
		    static_cast<long>(qgrams.size()) - static_cast<long>(k) * BED_Q;

		std::vector<int> counts(dict.size(), 0);
		std::vector<std::size_t> touched;
		touched.reserve(dict.size());
		for (const auto &g : qgrams) {
			const auto it = postings.find(g);
			if (it == postings.end())
				continue;
			for (const std::size_t id : it->second) {
				if (counts[id] == 0)
					touched.push_back(id);
				++counts[id];
			}
		}

		std::vector<Result> results;
		for (std::size_t i = 0; i < dict.size(); ++i) {
			if (std::labs(static_cast<long>(dict[i].size()) -
			              static_cast<long>(query_str.size())) > k)
				continue;
			if (counts[i] < need)
				continue; // count filter: shares too few q-grams
			const int dist =
			    edit_distance(query_str, dict[i], cost, k, c_min, g_min);
			if (dist <= k)
				results.push_back({i, dist});
		}
		return results; // ascending ids: already sorted
	}

private:
	static constexpr std::size_t BED_Q = 3; ///< q-gram length.
	/// q-gram -> sorted ids of dictionary sequences containing it.
	std::unordered_map<std::string, std::vector<std::size_t>> postings;
};

/**
 * @brief Vantage-point tree for exact threshold edit distance search.
 *
 * vpsearch (https://github.com/enthought/vpsearch) ships as Python/Cython
 * on top of parasail alignment scores, which cannot plug into this C++
 * benchmark; bench reimplements its core algorithm natively over the
 * Levenshtein edit distance instead: the tree is built by picking a
 * vantage point per node (fixed-seed RNG, as in vpsearch's default
 * RandomState) and splitting the remaining points at the median distance
 * mu, and range queries prune with the triangle inequality (inside iff
 * d <= mu + k, outside iff d >= mu - k). Internal vantage distances use
 * the exact (uncapped) edit distance, since pruning needs true values;
 * matches are reported with d <= k. Exact for any metric cost model;
 * bench restricts it to the unit-cost model (validated in main).
 */
class VpSearchTree {
public:
	VpSearchTree(const std::vector<std::string> &dict, const CostModel &cost,
	             const std::string &alphabet)
	    : dict(dict), cost(cost), alphabet(alphabet), rng(VP_SEED) {
		std::vector<std::size_t> ids(dict.size());
		std::iota(ids.begin(), ids.end(), 0);
		root = build(std::move(ids));
	}

	std::vector<Result> query(const std::string &query_str, int k) const {
		std::vector<Result> results;
		std::vector<const Node *> stack;
		if (root)
			stack.push_back(root.get());
		while (!stack.empty()) {
			const Node *node = stack.back();
			stack.pop_back();
			const int d = exact(query_str, dict[node->id]);
			if (d <= k)
				results.push_back({node->id, d});
			if (node->inside && d <= node->mu + k)
				stack.push_back(node->inside.get());
			if (node->outside && d >= node->mu - k)
				stack.push_back(node->outside.get());
		}
		std::sort(results.begin(), results.end(),
		          [](const Result &a, const Result &b) { return a.id < b.id; });
		return results;
	}

private:
	struct Node {
		std::size_t id; ///< Dictionary index of the vantage point.
		int mu = 0;     ///< Median distance from the vantage point.
		std::unique_ptr<Node> inside, outside;
	};

	static constexpr std::uint32_t VP_SEED = 0x8d278e75; ///< Fixed RNG seed.

	/** Exact (uncapped) edit distance: k = |a| + |b| exceeds any alignment. */
	int exact(const std::string &a, const std::string &b) const {
		const auto [c_min, g_min] = model_mins(cost, alphabet);
		const int cap =
		    static_cast<int>(a.size() + b.size()); // exact upper bound
		return edit_distance(a, b, cost, cap, c_min, g_min);
	}

	std::unique_ptr<Node> build(std::vector<std::size_t> ids) {
		if (ids.empty())
			return nullptr;
		std::uniform_int_distribution<std::size_t> pick(0, ids.size() - 1);
		const std::size_t vantage = ids[pick(rng)];
		ids.erase(std::remove(ids.begin(), ids.end(), vantage), ids.end());
		auto node = std::make_unique<Node>();
		node->id = vantage;
		if (ids.empty())
			return node;
		std::vector<int> dists;
		dists.reserve(ids.size());
		for (const std::size_t id : ids)
			dists.push_back(exact(dict[vantage], dict[id]));
		std::vector<int> order = dists; // median without disturbing ids order
		const std::size_t mid = order.size() / 2;
		std::nth_element(order.begin(), order.begin() + mid, order.end());
		node->mu = order[mid];
		std::vector<std::size_t> inner, outer;
		inner.reserve(mid);
		outer.reserve(ids.size() - mid);
		for (std::size_t i = 0; i < ids.size(); ++i) {
			if (dists[i] < node->mu)
				inner.push_back(ids[i]);
			else
				outer.push_back(ids[i]);
		}
		node->inside = build(std::move(inner));
		node->outside = build(std::move(outer));
		return node;
	}

	const std::vector<std::string> &dict;
	const CostModel &cost;
	const std::string &alphabet;
	std::mt19937 rng;
	std::unique_ptr<Node> root;
};

/**
 * @brief Runs the selected method and streams per-query records to @p fd.
 *
 * First line: `BUILD <build_ms>`; then one `QUERY <index> <k> <ms> <hash>
 * <matches>` line per query, in testcase order. When @p index_file is set,
 * the hlp_grep method reconstructs its solver from that prebuilt index
 * (written by `save`) instead of building the dictionary, and reports a
 * zero build time like the index-free methods.
 */
void run_method(int fd, const std::string &method, const Testcase &tc,
                bool compact_nodes, const std::string &index_file) {
	if (method == "hlp_grep") {
		if (index_file.empty()) {
			const auto build_b = std::chrono::steady_clock::now();
			const Solver solver(tc.dict, tc.cost, compact_nodes);
			const auto build_e = std::chrono::steady_clock::now();
			write_all(fd, "BUILD " + std::to_string(elapsed_ms(build_b, build_e)) +
			                  "\n");
			report_hlp_queries(fd, solver, tc);
		} else {
			const Solver solver = load_hlp_index(index_file, tc);
			write_all(fd, "BUILD 0\n");
			report_hlp_queries(fd, solver, tc);
		}
#ifdef HLP_GREP_WITH_WFA2
	} else if (method == "wfa") {
		// Threshold edit distance search with WFA2-lib: align the query
		// against every dictionary sequence (Levenshtein distance via the
		// built-in edit metric) and keep scores <= k. No index to build;
		// the aligner is reused across pairs. Requires the unit-cost
		// model (validated in the parent process).
		//
		// Early stopping: in the edit metric the wavefront count equals
		// the score, so capping the aligner at k + 1 alignment steps
		// (wavefronts) aborts every pair whose edit distance exceeds k:
		// pairs with dist <= k always complete a full alignment within
		// the cap, while the rest return WF_STATUS_MAX_STEPS_REACHED
		// without the work of walking to their final score.
		write_all(fd, "BUILD 0\n");

		wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
		attributes.distance_metric = edit;
		attributes.alignment_scope = compute_score;
		attributes.memory_mode = wavefront_memory_high;
		wavefront_aligner_t *wf = wavefront_aligner_new(&attributes);
		if (wf == nullptr)
			::_exit(1);

		const auto report_candidates = [&](const std::string &query, int k) {
			wavefront_aligner_set_max_alignment_steps(wf, k + 1);
			std::vector<Result> results;
			for (std::size_t i = 0; i < tc.dict.size(); ++i) {
				const std::string &seq = tc.dict[i];
				if (std::abs(static_cast<long>(seq.size()) -
				             static_cast<long>(query.size())) > k)
					continue;
				const int status =
				    wavefront_align(wf, query.data(),
				                    static_cast<int>(query.size()),
				                    seq.data(), static_cast<int>(seq.size()));
				if (status == WF_STATUS_MAX_STEPS_REACHED) {
					continue; // edit distance exceeds k: threshold-pruned
				}
				if (status != WF_STATUS_ALG_COMPLETED &&
				    status != WF_STATUS_ALG_PARTIAL)
					::_exit(1);
				const int dist = wf->cigar->score;
				if (dist >= 0 && dist <= k)
					results.push_back({i, dist});
			}
			return results;
		};

		for (std::size_t i = 0; i < tc.queries.size(); ++i) {
			const auto [k, query] = tc.queries[i];
			const auto b = std::chrono::steady_clock::now();
			const auto results = report_candidates(query, k);
			const auto e = std::chrono::steady_clock::now();
			report_query(fd, i, k, results, b, e);
		}
		wavefront_aligner_delete(wf);
#endif
#ifdef HLP_GREP_WITH_DT_PATRICIA
	} else if (method == "dt_patricia") {
		// DT-Patricia: exact threshold edit distance search over a
		// Patricia tree with the diagonal-transition algorithm. The tree
		// is a shared-prefix index, so its construction cost is the
		// method's build time (unlike naive/wfa). Unit-cost model and
		// a supported alphabet (DNA ACGT or the 20 protein letters)
		// are validated in the parent process.
		using namespace dt_patricia;

		// The alphabet policy must cover every testcase character.
		const auto is_dna = [&] {
			for (const char c : tc.alphabet)
				if (std::string_view("ACGT").find(c) == std::string_view::npos)
					return false;
			return true;
		}();

		const auto build_b = std::chrono::steady_clock::now();
		if (is_dna) {
			PatriciaTree<DnaAlphabet> tree(tc.dict);
			DTPatricia<DnaAlphabet, UnitCost> aligner(tree);
			const auto build_e = std::chrono::steady_clock::now();
			write_all(fd, "BUILD " +
			                  std::to_string(elapsed_ms(build_b, build_e)) +
			                  "\n");

			for (std::size_t i = 0; i < tc.queries.size(); ++i) {
				const auto [k, query] = tc.queries[i];
				const auto b = std::chrono::steady_clock::now();
				auto results = aligner.ed_within_k(query, k);
				std::sort(results.begin(), results.end(),
				          [](const AlignmentResult &a,
				             const AlignmentResult &b) {
					          return a.string_id < b.string_id;
				          });
				const auto e = std::chrono::steady_clock::now();
				report_query(fd, i, k, from_dt_patricia(results), b, e);
			}
		} else {
			PatriciaTree<ProteinAlphabet> tree(tc.dict);
			DTPatricia<ProteinAlphabet, UnitCost> aligner(tree);
			const auto build_e = std::chrono::steady_clock::now();
			write_all(fd, "BUILD " +
			                  std::to_string(elapsed_ms(build_b, build_e)) +
			                  "\n");

			for (std::size_t i = 0; i < tc.queries.size(); ++i) {
				const auto [k, query] = tc.queries[i];
				const auto b = std::chrono::steady_clock::now();
				auto results = aligner.ed_within_k(query, k);
				std::sort(results.begin(), results.end(),
				          [](const AlignmentResult &a,
				             const AlignmentResult &b) {
					          return a.string_id < b.string_id;
				          });
				const auto e = std::chrono::steady_clock::now();
				report_query(fd, i, k, from_dt_patricia(results), b, e);
			}
		}
#endif
	} else if (method == "bed_tree") {
		// B^ed-tree gram-count range query (see BedTreeIndex): the q-gram
		// inverted index is the method's build cost; each query applies
		// the length and count filters and verifies survivors exactly.
		// Unit-cost model validated in the parent process.
		const auto build_b = std::chrono::steady_clock::now();
		const BedTreeIndex index(tc.dict);
		const auto build_e = std::chrono::steady_clock::now();
		write_all(fd, "BUILD " +
		                  std::to_string(elapsed_ms(build_b, build_e)) + "\n");

		for (std::size_t i = 0; i < tc.queries.size(); ++i) {
			const auto [k, query] = tc.queries[i];
			const auto b = std::chrono::steady_clock::now();
			const auto results =
			    index.query(tc.dict, tc.cost, tc.alphabet, query, k);
			const auto e = std::chrono::steady_clock::now();
			report_query(fd, i, k, results, b, e);
		}
	} else if (method == "vpsearch") {
		// Vantage-point tree range search (see VpSearchTree): the tree is
		// the method's build cost; queries prune with the triangle
		// inequality and report exact distances. Unit-cost model
		// validated in the parent process.
		const auto build_b = std::chrono::steady_clock::now();
		const VpSearchTree tree(tc.dict, tc.cost, tc.alphabet);
		const auto build_e = std::chrono::steady_clock::now();
		write_all(fd, "BUILD " +
		                  std::to_string(elapsed_ms(build_b, build_e)) + "\n");

		for (std::size_t i = 0; i < tc.queries.size(); ++i) {
			const auto [k, query] = tc.queries[i];
			const auto b = std::chrono::steady_clock::now();
			const auto results = tree.query(query, k);
			const auto e = std::chrono::steady_clock::now();
			report_query(fd, i, k, results, b, e);
		}
	} else { // naive: no index to build.
		write_all(fd, "BUILD 0\n");

		const NaiveSolver solver(tc.dict, tc.cost, tc.alphabet);
		for (std::size_t i = 0; i < tc.queries.size(); ++i) {
			const auto [k, query] = tc.queries[i];
			const auto b = std::chrono::steady_clock::now();
			const auto results = solver.query(query, k);
			const auto e = std::chrono::steady_clock::now();
			report_query(fd, i, k, results, b, e);
		}
	}

	// _exit: do not run the parent's inherited atexit/stdio flushes.
	::_exit(0);
}

/**
 * @brief Runs one method in a forked child and collects its records.
 *
 * @return Build time, peak memory in bytes and one QueryStat per query.
 *         Reports failure on stderr and returns nullopt on any error.
 */
std::optional<MethodRun> measure_method(const std::string &method,
                                         const Testcase &tc,
                                         bool compact_nodes,
                                         const std::string &index_file) {
	int pipefd[2];
	if (pipe(pipefd) != 0) {
		std::cerr << "pipe() failed\n";
		return std::nullopt;
	}
	const pid_t pid = fork();
	if (pid < 0) {
		std::cerr << "fork() failed\n";
		return std::nullopt;
	}
	if (pid == 0) {
		close(pipefd[0]);
		run_method(pipefd[1], method, tc, compact_nodes,
		           index_file); // never returns
	}
	close(pipefd[1]);

	std::string data;
	char buf[4096];
	long n;
	while ((n = read(pipefd[0], buf, sizeof(buf))) > 0)
		data.append(buf, static_cast<std::size_t>(n));
	close(pipefd[0]);

	int status;
	struct rusage ru;
	if (wait4(pid, &status, 0, &ru) < 0) {
		std::cerr << "wait4() failed\n";
		return std::nullopt;
	}
	if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
		std::cerr << method << " run failed\n";
		return std::nullopt;
	}
	// Linux reports ru_maxrss in KiB.
	const long peak_bytes = ru.ru_maxrss * 1024L;

	// Parse the child's protocol: first "BUILD <ms>", then one QUERY line
	// per query in testcase order.
	double build_ms = 0.0;
	std::vector<QueryStat> qstats;
	{
		std::size_t pos = 0;
		auto next_line = [&data, &pos]() {
			const std::size_t nl = data.find('\n', pos);
			if (nl == std::string::npos)
				return std::string();
			std::string s = data.substr(pos, nl - pos);
			pos = nl + 1;
			return s;
		};
		if (sscanf(data.c_str(), "BUILD %lf", &build_ms) != 1) {
			std::cerr << "missing BUILD record from " << method << " child\n";
			return std::nullopt;
		}
		while (pos < data.size()) {
			QueryStat q{};
			if (sscanf(next_line().c_str(), "QUERY %zu %d %lf %lu %lu",
			           &q.index, &q.k, &q.runtime_ms, &q.hash, &q.num_matches)
			    == 5)
				qstats.push_back(q);
		}
	}
	if (qstats.size() != tc.queries.size()) {
		std::cerr << method << ": received " << qstats.size() << " of "
		          << tc.queries.size() << " query records; child aborted early?\n";
		return std::nullopt;
	}
	std::sort(qstats.begin(), qstats.end(),
	          [](const QueryStat &a, const QueryStat &b) {
		          return a.index < b.index;
	          });
	return MethodRun{method, build_ms, peak_bytes, std::move(qstats)};
}

[[noreturn]] void usage(const char *argv0) {
	std::cerr << "usage: " << argv0
	          << " --method hlp_grep [--method naive|wfa|dt_patricia|bed_tree|vpsearch ...]"
	          << " <testcase-file> [--out out.json] [--no-compact]"
	          << " [--index graph.hlpg]\n";
	std::exit(1);
}

/** True when the model is unary unit-cost (WFA2 edit metric). */
bool unit_cost_model(const Testcase &tc) {
	const std::string &alpha = tc.alphabet;
	for (const char a : alpha) {
		if (tc.cost.ins() != 1 || tc.cost.del() != 1)
			return false;
		for (const char b : alpha) {
			const int g = tc.cost.consume(a, b);
			if (a == b ? g != 0 : g != 1)
				return false;
		}
	}
	return true;
}

/** True when every alphabet character maps in the chosen policy sets. */
bool supported_alphabet(const Testcase &tc) {
	for (const char c : tc.alphabet) {
		const char up = static_cast<char>(std::toupper(
		    static_cast<unsigned char>(c)));
		if (std::string_view("ACGT").find(up) == std::string_view::npos &&
		    std::string_view("ACDEFGHIKLMNPQRSTVWY").find(up) ==
		        std::string_view::npos)
			return false;
	}
	return true;
}

} // namespace

int main(int argc, char **argv) {
	std::vector<std::string> method_names;
	std::string testcase, out_path, index_file;
	bool compact_nodes = true;
	for (int i = 1; i < argc; ++i) {
		const std::string arg = argv[i];
		if (arg == "--method") {
			if (i + 1 >= argc)
				usage(argv[0]);
			method_names.push_back(argv[++i]);
		} else if (arg == "--out") {
			if (i + 1 >= argc)
				usage(argv[0]);
			out_path = argv[++i];
		} else if (arg == "--index") {
			if (i + 1 >= argc)
				usage(argv[0]);
			index_file = argv[++i];
		} else if (arg == "--no-compact") {
			compact_nodes = false; // POAGraph run compaction off
		} else if (!arg.empty() && arg[0] != '-' && testcase.empty()) {
			testcase = arg;
		} else {
			std::cerr << "unexpected argument: " << arg << '\n';
			usage(argv[0]);
		}
	}
	if (method_names.empty()) {
		std::cerr << "missing --method (repeatable)\n";
		usage(argv[0]);
	}
	for (const auto &method : method_names) {
		if (method != "hlp_grep" && method != "naive" && method != "wfa" &&
		    method != "dt_patricia" && method != "bed_tree" &&
		    method != "vpsearch") {
			std::cerr << "unknown method: '" << method << "'\n";
			usage(argv[0]);
		}
	}
#ifndef HLP_GREP_WITH_WFA2
	if (std::find(method_names.begin(), method_names.end(), "wfa")
	    != method_names.end()) {
		std::cerr << "method 'wfa' unavailable: bench was built without "
		          << "WFA2 (tests/scripts/setup_wfa2.sh or cmake with "
		          << "-DHLP_GREP_WITH_WFA2=OFF)\n";
		return 1;
	}
#endif
#ifndef HLP_GREP_WITH_DT_PATRICIA
	if (std::find(method_names.begin(), method_names.end(), "dt_patricia")
	    != method_names.end()) {
		std::cerr << "method 'dt_patricia' unavailable: bench was built "
		          << "without DT-Patricia (tests/scripts/setup_dt_patricia.sh "
		          << "or cmake with -DHLP_GREP_WITH_DT_PATRICIA=OFF)\n";
		return 1;
	}
#endif

	if (testcase.empty()) {
		std::cerr << "missing testcase file\n";
		usage(argv[0]);
	}
	const fs::path file(testcase);
	if (!fs::is_regular_file(file)) {
		std::cerr << "not a testcase file: " << file << '\n';
		return 1;
	}
	if (out_path.empty()) {
		fs::path out = file;
		out.replace_extension(".bench.json");
		out_path = out.string();
	}

	const Testcase tc = parse_testcase(file);

	// A prebuilt index (written by `save`) only feeds the hlp_grep method:
	// it reconstructs the solver without building the dictionary, so the
	// reported build time is zero.
	if (!index_file.empty()) {
		const bool wants_hlp =
		    std::find(method_names.begin(), method_names.end(),
		              "hlp_grep") != method_names.end();
		if (!wants_hlp) {
			std::cerr << "--index applies only to --method hlp_grep\n";
			return 1;
		}
		if (!fs::is_regular_file(fs::path(index_file))) {
			std::cerr << "not an index file: " << index_file << '\n';
			return 1;
		}
		if (!compact_nodes)
			std::cerr << "note: --no-compact has no effect with --index "
			          << "(compaction is stored in the index file)\n";
	}

	// WFA2, DT-Patricia, BED-tree and VP-search as used here implement
	// the unit-cost Levenshtein distance only.
	const auto wants = [&](const char *m) {
		return std::find(method_names.begin(), method_names.end(), m) !=
		       method_names.end();
	};
	const bool wants_uniform =
	    wants("wfa") || wants("dt_patricia") || wants("bed_tree") || wants("vpsearch");
	if (wants_uniform && !unit_cost_model(tc)) {
		std::cerr << "method(s) 'wfa'/'dt_patricia'/'bed_tree'/'vpsearch' "
		          << "require a unit-cost model (match 0, ins/del/mismatch 1) "
		          << "but " << file << " uses a different cost model\n";
		return 1;
	}
#ifdef HLP_GREP_WITH_DT_PATRICIA
	if ((wants("wfa") || wants("dt_patricia")) && !supported_alphabet(tc)) {
		std::cerr << "method(s) 'wfa'/'dt_patricia' require an alphabet "
		          << "of DNA (ACGT) or protein (20 letters) characters; "
		          << "testcase alphabet is '" << tc.alphabet << "'\n";
		return 1;
	}
#endif

	std::vector<MethodRun> runs;
	for (const auto &method : method_names) {
		auto run = measure_method(method, tc, compact_nodes, index_file);
		if (!run)
			return 1;
		runs.push_back(std::move(*run));
	}

	// Dictionary and query lengths, reported alongside the timings.
	std::vector<std::size_t> dict_lens, query_lens;
	dict_lens.reserve(tc.dict.size());
	for (const auto &seq : tc.dict)
		dict_lens.push_back(seq.size());
	query_lens.reserve(tc.queries.size());
	for (const auto &[k, query] : tc.queries) {
		(void)k;
		query_lens.push_back(query.size());
	}
	auto sorted = dict_lens;
	std::sort(sorted.begin(), sorted.end());

	std::ostringstream json;
	json << std::fixed << std::setprecision(6);
	json << "{\n"
	     << "  \"testcase\": \"" << file.string() << "\",\n"
	     << "  \"compact_nodes\": " << (compact_nodes ? "true" : "false" ) << ",\n"
	     << "  \"dict_size\": " << tc.dict.size() << ",\n"
	     << "  \"num_queries\": " << tc.queries.size() << ",\n";
	if (sorted.empty()) {
		json << "  \"dict_len_min\": 0, \"dict_len_max\": 0,\n"
		     << "  \"dict_len_median\": 0, \"dict_len_mean\": 0,\n";
	} else {
		const long sum = std::accumulate(sorted.begin(), sorted.end(), 0L);
		json << "  \"dict_len_min\": " << sorted.front() << ",\n"
		     << "  \"dict_len_max\": " << sorted.back() << ",\n"
		     << "  \"dict_len_median\": " << sorted[sorted.size() / 2] << ",\n"
		     << "  \"dict_len_mean\": "
		     << static_cast<double>(sum) / sorted.size() << ",\n";
	}
	json << "  \"dict_len\": [";
	for (std::size_t i = 0; i < dict_lens.size(); ++i)
		json << (i ? ", " : "") << dict_lens[i];
	json << "],\n"
	     << "  \"query_len\": [";
	for (std::size_t i = 0; i < query_lens.size(); ++i)
		json << (i ? ", " : "") << query_lens[i];
	json << "],\n";

	json << "  \"methods\": [\n";
	for (std::size_t m = 0; m < runs.size(); ++m) {
		const auto &r = runs[m];
		json << "    {\n"
		     << "      \"method\": \"" << r.name << "\",\n"
		     << "      \"build_time_ms\": " << r.build_ms << ",\n"
		     << "      \"peak_memory_bytes\": " << r.peak_bytes << ",\n"
		     << "      \"queries\": [\n";
		for (std::size_t i = 0; i < r.queries.size(); ++i) {
			const auto &q = r.queries[i];
			json << "        {\"index\": " << q.index << ", \"k\": " << q.k
			     << ", \"len\": " << query_lens[q.index]
			     << ", \"runtime_ms\": " << q.runtime_ms
			     << ", \"num_matches\": " << q.num_matches
			     << ", \"results_hash\": " << q.hash << "}"
			     << (i + 1 < r.queries.size() ? "," : "") << "\n";
		}
		json << "      ]\n"
		     << "    }" << (m + 1 < runs.size() ? "," : "") << "\n";
	}
	json << "  ]\n"
	     << "}\n";

	std::ofstream out(out_path);
	if (!out) {
		std::cerr << "cannot write " << out_path << '\n';
		return 1;
	}
	out << json.str() << std::flush;
	if (!out) {
		std::cerr << "write failed while writing " << out_path << '\n';
		return 1;
	}
	std::cout << "BENCH " << file << " -> " << out_path
	          << " (dict=" << tc.dict.size() << " queries=" << tc.queries.size()
	          << " methods=";
	for (std::size_t m = 0; m < runs.size(); ++m)
		std::cout << (m ? "," : "") << runs[m].name;
	std::cout << ")\n";
	return 0;
}
