/**
 * @file hlp_grep.hpp
 * @brief Edit distance search queries on a DNA sequence dictionary using
 *        heavy-light decomposition on POA graph paths.
 *
 * HLP-grep indexes a dictionary of DNA sequences (represented as paths in a
 * partial order alignment (POA) graph) and answers edit distance search
 * queries: given a query string and a threshold k, it finds all dictionary
 * sequences within edit distance k of the query string.
 */
#pragma once
#include <hlp_grep/binary_lifter.hpp>
#include <hlp_grep/cost_model.hpp>
#include <hlp_grep/dist_matrix.hpp>
#include <hlp_grep/poa_graph.hpp>
#include <hlp_grep/result.hpp>

#include <algorithm>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace hlp_grep {

/**
 * @brief Edit distance search engine over a sequence dictionary using
 *        heavy-light decomposition on POA graph paths.
 *
 * Indexes the dictionary sequences (represented as paths in a POA graph) and
 * answers edit distance search queries.
 *
 * Usage: construct the solver with the dictionary, then call query() to
 * retrieve all dictionary sequences within a given edit distance threshold
 * of a query string.
 */
class Solver {
public:
	/**
	 * @brief Constructs the solver over the given sequence dictionary.
	 *
	 * @param dict Dictionary of DNA sequences to search. The position of each
	 *             sequence in this vector defines the `id` reported in Result.
	 * @param cost Cost model defining the costs of the basic edit operations;
	 *             defaults to the unit-cost model. The referenced model must
	 *             outlive the solver.
	 * @param compact_nodes Whether to merge single-in/single-out node runs
	 *             into multi-character nodes after graph construction.
	 * @param num_threads Worker threads for the per-query BinaryLifter
	 *             chain table build; 0 means use the hardware concurrency.
	 */
	explicit Solver(std::vector<std::string> dict,
	                const CostModel &cost = DEFAULT_COST_MODEL,
	                bool compact_nodes = true,
	                std::size_t num_threads = 0)
	    : dict(std::move(dict)), cost(cost), num_threads(num_threads),
	      graph(this->dict, this->cost, compact_nodes) {
		build_compressed_paths();
	}

	/**
	 * @brief Finds all dictionary sequences within edit distance k of the
	 *        query.
	 *
	 * Builds a BinaryLifter over the graph for this query, then scores
	 * every dictionary sequence's compressed path: the DP row starts at
	 * the path's start node and each compressed step advances it across a
	 * heavy chain (jump) or a single light edge (step). The last entry of
	 * the final row is the edit distance between the sequence's path and
	 * the query; it is reported when it does not exceed k.
	 *
	 * Scoring aborts a path early once every entry of its DP row exceeds
	 * k: with nonnegative operation costs each transition only adds cost,
	 * so the row minimum never decreases and no later state can fall back
	 * to <= k. Costs are nonnegative in every sane model; negative-cost
	 * models make the check unsound and must not rely on it.
	 *
	 * Sequences whose length differs from the query by more than k are
	 * skipped without scoring: their length difference is a lower bound
	 * on the edit distance, and their DP states need not fit the windows
	 * the lifter precomputes.
	 *
	 * @param query The query DNA sequence.
	 * @param k     Maximum allowed edit distance (threshold).
	 * @return A vector of Result, one per matched dictionary sequence,
	 *         each containing the sequence's index in the dictionary and its
	 *         edit distance to the query. Results are sorted by `id` in
	 *         ascending order.
	 */
	std::vector<Result> query(const std::string &query, int k) const {
		std::vector<Result> results;
		if (dict.empty())
			return results;
		const long m = static_cast<long>(query.size());
		BinaryLifter lifter(graph, query, k, num_threads);
		// Every path is scored independently: each iteration reads const
		// graph/lifter state and its own compressed path, and writes only
		// its own slot of dist.
		std::vector<int> dist(dict.size(), DistMatrix::INF);
#pragma omp parallel for schedule(dynamic)
		for (std::size_t i = 0; i < dict.size(); ++i) {
			if (std::abs(static_cast<long>(dict[i].size()) - m) > k)
				continue;
			const auto &cp = compressed_paths[i];
			std::vector<int> row = initial_row(cp.start, query, k);
			POAGraph::node_id cur = cp.start;
			for (const auto &st : cp.steps) {
				if (st.type == POAGraph::EdgeType::HEAVY) {
					cur = lifter.jump(cur, st.length, row);
				} else {
					row = lifter.step(cur, st.next, row);
					cur = st.next;
				}
				if (*std::min_element(row.begin(), row.end()) > k)
					break;
			}
			dist[i] = row.back();
		}
		for (std::size_t i = 0; i < dict.size(); ++i)
			if (dist[i] <= k)
				results.push_back({i, dist[i]});
		return results;
	}

private:
	/**
	 * @brief Initial DP row at the start node of a compressed path.
	 *
	 * Row entries are the node's window positions: the number of query
	 * characters consumed before the start node's label begins. The
	 * label is aligned as a block against some query substring ending at
	 * b (matches, substitutions, and internal deletions via the ED
	 * recurrence below), with every earlier query character inserted.
	 * The window must match BinaryLifter::window (the start node's label
	 * is never crossed by an edge, so no matrix can produce this row).
	 *
	 * @param start First node of the compressed path.
	 * @param query Query string (may be empty).
	 * @param k     Edit distance threshold widening the window.
	 * @return One entry per position of start's window.
	 */
	std::vector<int> initial_row(POAGraph::node_id start,
	                             const std::string &query, int k) const {
		const auto [jmin, jmax] = graph.pos_range(start);
		const long m = static_cast<long>(query.size());
		const long mk = static_cast<long>(k);
		const std::string &label = graph.seq(start);
		const long label_len = static_cast<long>(label.size());
		// After-label band of the start node: the output window, kept
		// non-empty to match BinaryLifter::window (start nodes have
		// jmin = 0, so this is [max(0, |label| - k), min(|query| + 1,
		// jmax + |label| + k + 1))).
		const long lo = std::max(0L, static_cast<long>(jmin) + label_len - mk);
		const long hL =
		    std::min(m + 1, static_cast<long>(jmax) + label_len + mk + 1);
		const long hi = std::max(hL, lo + 1);
		std::vector<int> row(static_cast<std::size_t>(hi - lo),
		                     DistMatrix::INF);
		// Semiglobal DP over (label prefix, query prefix): f[i][b] is the
		// min cost of consuming the first i label characters against a
		// query substring ending at b. Row 0 is the pure insertion prefix
		// from the virtual source at position 0; each label character
		// sweeps over its own band [max(0, jmin + i - k), min(m + 1,
		// jmax + i + k + 1)) — mirroring BinaryLifter::sweep_label — so
		// the work stays O(|label| * (spread + 2k)). Entries saturate at
		// INF; when the after-label band is empty every band is, and the
		// row stays all-INF.
		if (lo >= hL)
			return row;
		const long l0 = std::max(0L, static_cast<long>(jmin) - mk);
		const long h0 = std::min(m + 1, static_cast<long>(jmax) + mk + 1);
		std::vector<long> f(static_cast<std::size_t>(hL - l0),
		                    DistMatrix::INF);
		// Virtual source at position 0 (start nodes have jmin = 0, so the
		// seed sits at the entry band's bottom); chain insertions within
		// the entry band [l0, h0).
		f[0] = 0;
		for (long b = 1; b < h0; ++b)
			f[b - l0] = std::min(f[b - 1 - l0] + cost.ins(),
			                     static_cast<long>(DistMatrix::INF));
		for (long i = 1; i <= label_len; ++i) {
			const long li =
			    std::max(0L, static_cast<long>(jmin) + i - mk);
			const long hi_i =
			    std::min(m + 1, static_cast<long>(jmax) + i + mk + 1);
			const char c = label[i - 1];
			const int del_c = cost.del();
			long diag = li - 1 >= l0 ? f[li - 1 - l0] : DistMatrix::INF;
			long left = DistMatrix::INF; // f[i][b - 1]: below-band at li
			for (long b = li; b < hi_i; ++b) {
				const long up = f[b - l0]; // f[i - 1][b]
				long cur = up + del_c;
				if (left < DistMatrix::INF)
					cur = std::min(cur, left + cost.ins());
				if (diag < DistMatrix::INF)
					cur = std::min(cur,
					               diag + cost.consume(c,
					                                   query[b - 1]));
				if (cur > DistMatrix::INF)
					cur = DistMatrix::INF;
				f[b - l0] = cur;
				left = cur;
				diag = up;
			}
		}
		for (long b = lo; b < hi; ++b)
			row[static_cast<std::size_t>(b - lo)] =
			    static_cast<int>(f[b - l0]);
		return row;
	}

	std::vector<std::string> dict; ///< Dictionary of DNA sequences to search.
	/// Cost model used for the edit distance computations; must outlive
	/// the solver.
	const CostModel &cost;
	/// Worker threads for the per-query BinaryLifter build (0 = auto).
	std::size_t num_threads = 0;
	POAGraph graph;                ///< POA graph built from the dictionary.
	/// Compressed (heavy-chain / light-step) representation of each dict path.
	std::vector<POAGraph::CompressedPath> compressed_paths;

	/** Builds the compressed representation of every dictionary path. */
	void build_compressed_paths() {
		compressed_paths.clear();
		compressed_paths.reserve(dict.size());
		for (std::size_t i = 0; i < dict.size(); ++i)
			compressed_paths.push_back(graph.compressed_path(i));
	}
};

} // namespace hlp_grep
