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

// Vendored BED-tree headers (legacy code with `using namespace std` and
// generic global class names like Node/Page/Queue): include last, and keep
// this file's own helper names collision-free below.
#ifdef HLP_GREP_WITH_BED_TREE
#include <BufferManager.h>
#include <SBPTree.h>
#include <StrDictOrder.h>
#endif

// Vendored hstree header (legacy code with `using namespace std` and
// generic globals like dict/N/results): same collision discipline as above.
#ifdef HLP_GREP_WITH_HSTREE
#include <hstree.h>
#endif

#ifdef HLP_GREP_WITH_EDLIB
#include <edlib.h>
#endif

#ifdef HLP_GREP_WITH_PARASAIL
#include <parasail.h>
#endif

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <cctype>
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

#ifdef HLP_GREP_WITH_BED_TREE
/**
 * @brief Tuning for the vendored BED-tree backend (see run_method).
 *
 * The buffer table holds BED_BUF_ENTRIES nodes / BED_BUF_PAGES open pages;
 * the ctor preallocates entries * page bytes of backing file, so the entry
 * count stays modest (8 MiB at 8 KiB pages). Everything benchmark-sized
 * fits in memory and the disk file only sees evictions. Pages must fit the
 * longest headed entry, so the page size scales up from BED_PAGE_KB
 * minimum for long-sequence testcases (upstream's fixed MAXLEN/page
 * assumes short strings and cannot index entries larger than a page).
 * A page must hold at least three max-length entries: FindSplitPos
 * halves pages by bytes, so with room for fewer the split moves every
 * entry across and makes no progress (infinite split loop). Note the
 * library's per-node LowerBound is O(q^2), so long sequences stay slow
 * regardless — this backend targets the short-string regime.
 */
constexpr int BED_PAGE_KB = 8;
constexpr int BED_BUF_ENTRIES = 1024;
constexpr int BED_BUF_PAGES = 1024;

/** Longest raw string the BED-tree DP matrix must cover (dict + queries). */
std::size_t bed_tree_matrix_size(const Testcase &tc) {
	std::size_t m = 0;
	for (const auto &seq : tc.dict)
		m = std::max(m, seq.size());
	for (const auto &[k, query] : tc.queries) {
		(void)k;
		m = std::max(m, query.size());
	}
	return m;
}

/**
 * @brief Lowercases/prepares a headed entry for the BED-tree comparator.
 *
 * BED-tree strips a leading numeric head ("<id> <seq>") during
 * comparisons, but RemoveHead unconditionally scans for a blank, so every
 * string handed to it needs one; heads also carry the result ids back out
 * (see the bench patch in setup_bed_tree.sh). Returns the headless length.
 */
int bed_tree_prepare(StrDictOrder &comp, std::string entry,
                     std::vector<char> &buf) {
	if (buf.size() < entry.size() + 1)
		buf.resize(entry.size() + 1);
	return comp.OrgProcessing(entry.data(), buf.data());
}
#endif

#ifdef HLP_GREP_WITH_HSTREE
/**
 * @brief Exact threshold search with the vendored hstree library.
 *
 * hstree sorts its dictionary by (length, lexicographic) and reports
 * sorted-position ids, so the glue sorts the same way and maps hits back
 * to dictionary ids (ties need no care: every entry keeps its own id).
 * Empty queries bypass the index (strHash would read out of bounds on
 * them) and are answered directly from lengths: ed("", s) == |s|.
 * Thresholds beyond the indexed substring levels (lv(tau) >= maxLevel)
 * fall back to exact caled() verification over the dictionary.
 * Unit-cost model only (validated in main).
 */
class HSTreeIndex {
public:
	explicit HSTreeIndex(const std::vector<std::string> &dict) {
		order.resize(dict.size());
		std::iota(order.begin(), order.end(), 0);
		// Same order as hstree's sortData: by length, then lexicographic.
		std::sort(order.begin(), order.end(),
		          [&](std::size_t a, std::size_t b) {
			          if (dict[a].size() != dict[b].size())
				          return dict[a].size() < dict[b].size();
			          return dict[a] < dict[b];
		          });
		dict_lens.reserve(dict.size());
		for (std::size_t i = 0; i < dict.size(); ++i)
			dict_lens.push_back(static_cast<int>(dict[i].size()));
		::dict.clear();
		::dict.reserve(dict.size());
		::minDictlen = 0x7FFFFFFF;
		::maxDictlen = 0;
		for (const std::size_t id : order) {
			::dict.push_back(dict[id]);
			::minDictlen =
			    std::min(::minDictlen, static_cast<int>(dict[id].size()));
			::maxDictlen =
			    std::max(::maxDictlen, static_cast<int>(dict[id].size()));
		}
		::N = static_cast<int>(::dict.size());
		::initial();
		::createIndex();
		// hstree's substring tables only cover levels below maxLevel, so
		// larger thresholds would index out of bounds (lv(tau) grows with
		// tau); those queries fall back to exact verification below.
		max_tau = ::maxLevel >= 31 ? INT_MAX : ((1 << ::maxLevel) - 1);
	}

	std::vector<Result> query(const std::string &query_str, int k) const {
		std::vector<Result> results;
		// Empty queries bypass the index (see the class comment).
		if (query_str.empty()) {
			for (std::size_t i = 0; i < dict_lens.size(); ++i)
				if (dict_lens[i] <= k)
					results.push_back({i, dict_lens[i]});
			return results; // ascending ids: already sorted
		}
		if (k <= max_tau) {
			::presim(query_str, k);
			std::vector<std::pair<int, int>> hits;
			::hsearch(query_str, k, &hits);
			results.reserve(hits.size());
			for (const auto &[pos, dist] : hits)
				results.push_back({order[pos], dist});
		} else {
			for (std::size_t pos = 0; pos < order.size(); ++pos) {
				const int dist = ::caled(query_str, ::dict[pos], k);
				if (dist <= k)
					results.push_back({order[pos], dist});
			}
		}
		std::sort(results.begin(), results.end(),
		          [](const Result &a, const Result &b) {
			          return a.id < b.id;
		          });
		return results;
	}

private:
	std::vector<std::size_t> order; ///< Sorted position -> dictionary id.
	std::vector<int> dict_lens;     ///< Dictionary sequence lengths.
	int max_tau = 0; ///< Largest threshold the substring levels cover.
};
#endif

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
#ifdef HLP_GREP_WITH_BED_TREE
	} else if (method == "bed_tree") {
		// Real BED-tree (SBPTree, dictionary order from the SIGMOD 2010
		// paper): the tree over the headed dictionary is the method's
		// build cost; each query runs RangeQuery with the integer
		// threshold and converts the collected (id, dist) hits. The disk
		// backing file is sandboxed in a per-child temp dir via the
		// bench patch (BufferManager::diskPath). Unit-cost model
		// validated in the parent process.
		if (tc.dict.empty()) {
			write_all(fd, "BUILD 0\n");
			for (std::size_t i = 0; i < tc.queries.size(); ++i) {
				const auto [k, query] = tc.queries[i];
				(void)query;
				const auto now = std::chrono::steady_clock::now();
				report_query(fd, i, k, {}, now, now);
			}
		} else {
			char dir_template[] = "/tmp/hlp_bedtree_XXXXXX";
			if (mkdtemp(dir_template) == nullptr)
				::_exit(1);
			const std::string backing =
			    std::string(dir_template) + "/disk.dat";
			BufferManager::setDiskPath(backing.c_str());

			const auto build_b = std::chrono::steady_clock::now();
			const std::size_t max_len = bed_tree_matrix_size(tc);
			const int page_kb =
			    std::max(BED_PAGE_KB,
			             (3 * static_cast<int>(max_len) + 4096 + 1023) / 1024);
			StrDictOrder comp(static_cast<int>(max_len));
			BufferManager bm(BED_BUF_ENTRIES, BED_BUF_PAGES, page_kb);
			SBPTree tree(page_kb * 1024, &bm, &comp);
			std::vector<char> buf;
			for (std::size_t i = 0; i < tc.dict.size(); ++i) {
				bed_tree_prepare(comp,
				                 std::to_string(i) + " " + tc.dict[i], buf);
				tree.InsertString(buf.data());
			}
			const auto build_e = std::chrono::steady_clock::now();
			write_all(fd, "BUILD " +
			                  std::to_string(elapsed_ms(build_b, build_e)) +
			                  "\n");

			for (std::size_t i = 0; i < tc.queries.size(); ++i) {
				const auto [k, query] = tc.queries[i];
				const auto b = std::chrono::steady_clock::now();
				const int qlen = bed_tree_prepare(
				    comp, std::string("0 ") + query, buf);
				std::vector<std::pair<int, float>> hits;
				tree.RangeQuery(buf.data(), qlen, static_cast<float>(k),
				                &hits);
				std::vector<Result> results;
				results.reserve(hits.size());
				for (const auto &[id, dist] : hits)
					results.push_back({static_cast<std::size_t>(id),
					                   static_cast<int>(dist)});
				std::sort(results.begin(), results.end(),
				          [](const Result &a, const Result &b) {
					          return a.id < b.id;
				          });
				const auto e = std::chrono::steady_clock::now();
				report_query(fd, i, k, results, b, e);
			}

			unlink(backing.c_str());
			rmdir(dir_template);
		}
#endif
#ifdef HLP_GREP_WITH_HSTREE
	} else if (method == "hstree") {
		// Real hstree (threshold driver): the length-sorted dictionary +
		// substring index is the method's build cost; each query runs
		// presim/hsearch with the integer threshold and maps the
		// collected (sorted-pos, dist) hits back to dictionary ids.
		// Unit-cost model validated in the parent process.
		if (tc.dict.empty()) {
			write_all(fd, "BUILD 0\n");
			for (std::size_t i = 0; i < tc.queries.size(); ++i) {
				const auto [k, query] = tc.queries[i];
				(void)query;
				const auto now = std::chrono::steady_clock::now();
				report_query(fd, i, k, {}, now, now);
			}
		} else {
			const auto build_b = std::chrono::steady_clock::now();
			const HSTreeIndex index(tc.dict);
			const auto build_e = std::chrono::steady_clock::now();
			write_all(fd, "BUILD " +
			                  std::to_string(elapsed_ms(build_b, build_e)) +
			                  "\n");

			for (std::size_t i = 0; i < tc.queries.size(); ++i) {
				const auto [k, query] = tc.queries[i];
				const auto b = std::chrono::steady_clock::now();
				const auto results = index.query(query, k);
				const auto e = std::chrono::steady_clock::now();
				report_query(fd, i, k, results, b, e);
			}
		}
#endif
#ifdef HLP_GREP_WITH_EDLIB
	} else if (method == "edlib") {
		// Threshold edit distance search with edlib (Myers' bit-vector):
		// align the query against every dictionary sequence in NW global
		// mode and keep distances <= k (edlib reports -1 above the k
		// bound). No index to build. Unit-cost model validated in the
		// parent process.
		write_all(fd, "BUILD 0\n");

		EdlibAlignConfig config = edlibDefaultAlignConfig();
		config.mode = EDLIB_MODE_NW;
		config.task = EDLIB_TASK_DISTANCE;
		for (std::size_t i = 0; i < tc.queries.size(); ++i) {
			const auto [k, query] = tc.queries[i];
			const auto b = std::chrono::steady_clock::now();
			config.k = k;
			std::vector<Result> results;
			for (std::size_t j = 0; j < tc.dict.size(); ++j) {
				const std::string &seq = tc.dict[j];
				if (std::abs(static_cast<long>(seq.size()) -
				             static_cast<long>(query.size())) > k)
					continue;
				EdlibAlignResult result = edlibAlign(
				    query.data(), static_cast<int>(query.size()),
				    seq.data(), static_cast<int>(seq.size()), config);
				if (result.status != EDLIB_STATUS_OK)
					::_exit(1);
				if (result.editDistance >= 0 && result.editDistance <= k)
					results.push_back(
					    {j, static_cast<int>(result.editDistance)});
				edlibFreeAlignResult(result);
			}
			const auto e = std::chrono::steady_clock::now();
			report_query(fd, i, k, results, b, e);
		}
#endif
#ifdef HLP_GREP_WITH_PARASAIL
	} else if (method == "parasail") {
		// Threshold edit distance search with parasail's banded global
		// alignment: the unit-cost matrix (match 0, mismatch -1) with
		// affine gaps (open 1, extend 1, i.e. 1 per gapped base) turns
		// the banded NW score into exactly -ed, and the band of
		// half-width k keeps every alignment of cost <= k. No index to
		// build. The banded routine rejects empty strings and k == 0, so
		// those are answered directly (ed is the other length /
		// equality). Unit-cost model validated in the parent process.
		write_all(fd, "BUILD 0\n");

		parasail_matrix_t *matrix = parasail_matrix_create(
		    tc.alphabet.c_str(), 0, -1);
		if (matrix == nullptr)
			::_exit(1);
		for (std::size_t i = 0; i < tc.queries.size(); ++i) {
			const auto [k, query] = tc.queries[i];
			const auto b = std::chrono::steady_clock::now();
			std::vector<Result> results;
			for (std::size_t j = 0; j < tc.dict.size(); ++j) {
				const std::string &seq = tc.dict[j];
				if (std::abs(static_cast<long>(seq.size()) -
				             static_cast<long>(query.size())) > k)
					continue;
				if (query.empty() || seq.empty()) {
					// ed is the non-empty length here (length filter
					// above already guarantees it is <= k when one
					// side is empty... unless both are empty: dist 0).
					results.push_back(
					    {j, static_cast<int>(std::max(query.size(),
					                                  seq.size()))});
					continue;
				}
				if (k == 0) {
					if (query == seq)
						results.push_back({j, 0});
					continue;
				}
				parasail_result_t *result = parasail_nw_banded(
				    query.data(), static_cast<int>(query.size()),
				    seq.data(), static_cast<int>(seq.size()), 1, 1, k,
				    matrix);
				if (result == nullptr)
					::_exit(1);
				const int dist = -parasail_result_get_score(result);
				parasail_result_free(result);
				if (dist >= 0 && dist <= k)
					results.push_back({j, dist});
			}
			const auto e = std::chrono::steady_clock::now();
			report_query(fd, i, k, results, b, e);
		}
		parasail_matrix_free(matrix);
#endif
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
	          << " --method hlp_grep [--method naive|wfa|dt_patricia|bed_tree|hstree|edlib|parasail ...]"
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
		    method != "hstree" && method != "edlib" &&
		    method != "parasail") {
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
#ifndef HLP_GREP_WITH_BED_TREE
	if (std::find(method_names.begin(), method_names.end(), "bed_tree")
	    != method_names.end()) {
		std::cerr << "method 'bed_tree' unavailable: bench was built "
		          << "without BED-tree (tests/scripts/setup_bed_tree.sh "
		          << "or cmake with -DHLP_GREP_WITH_BED_TREE=OFF)\n";
		return 1;
	}
#endif
#ifndef HLP_GREP_WITH_HSTREE
	if (std::find(method_names.begin(), method_names.end(), "hstree")
	    != method_names.end()) {
		std::cerr << "method 'hstree' unavailable: bench was built "
		          << "without hstree (tests/scripts/setup_hstree.sh "
		          << "or cmake with -DHLP_GREP_WITH_HSTREE=OFF)\n";
		return 1;
	}
#endif
#ifndef HLP_GREP_WITH_EDLIB
	if (std::find(method_names.begin(), method_names.end(), "edlib")
	    != method_names.end()) {
		std::cerr << "method 'edlib' unavailable: bench was built "
		          << "without edlib (tests/scripts/setup_edlib.sh "
		          << "or cmake with -DHLP_GREP_WITH_EDLIB=OFF)\n";
		return 1;
	}
#endif
#ifndef HLP_GREP_WITH_PARASAIL
	if (std::find(method_names.begin(), method_names.end(), "parasail")
	    != method_names.end()) {
		std::cerr << "method 'parasail' unavailable: bench was built "
		          << "without parasail (tests/scripts/setup_parasail.sh "
		          << "or cmake with -DHLP_GREP_WITH_PARASAIL=OFF)\n";
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

	// WFA2, DT-Patricia, BED-tree, hstree, edlib and parasail as used here
	// implement the unit-cost Levenshtein distance only.
	const auto wants = [&](const char *m) {
		return std::find(method_names.begin(), method_names.end(), m) !=
		       method_names.end();
	};
	const bool wants_uniform = wants("wfa") || wants("dt_patricia") ||
	                           wants("bed_tree") || wants("hstree") ||
	                           wants("edlib") || wants("parasail");
	if (wants_uniform && !unit_cost_model(tc)) {
		std::cerr << "method(s) 'wfa'/'dt_patricia'/'bed_tree'/'hstree'/"
		          << "'edlib'/'parasail' require a unit-cost model (match 0, "
		          << "ins/del/mismatch 1) but " << file
		          << " uses a different cost model\n";
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
