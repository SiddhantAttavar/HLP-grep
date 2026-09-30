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
#include <cassert>
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
	 * @param num_threads Worker threads for per-query chain-table building
	 *             and path scoring; 0 means use the OpenMP maximum for scoring
	 *             (the chain-table build caps automatic parallelism at 8).
	 */
	explicit Solver(std::vector<std::string> dict,
	                const CostModel &cost = DEFAULT_COST_MODEL,
	                bool compact_nodes = true,
	                std::size_t num_threads = 0)
	    : dict(std::move(dict)), cost(cost), num_threads(num_threads),
	      graph(this->dict, this->cost, compact_nodes),
	      lifter(this->graph, num_threads) {
		build_compressed_paths();
		lifter.decompose(compressed_paths);
		sort_compressed_paths();
	}

	/**
	 * @brief Finds all dictionary sequences within edit distance k of the
	 *        query.
	 *
	 * Rebuilds the per-query chain tables of the pre-marked BinaryLifter
	 * (build(query, k) over the blocks decompose() touched — the touched set is
	 * query-independent, so every query fills the same blocks), then scores
	 * the length-filtered compressed paths in lexicographic order, divided
	 * into contiguous batches scored in parallel. Each batch keeps its own
	 * stack of DP rows for shared prefixes and computes only each path's
	 * unmatched suffix across heavy blocks (jump) or light edges (step).
	 * Prefixes at batch boundaries are recomputed. The final row's entry at
	 * query position |query| is the edit distance between the sequence's path
	 * and the query; it is reported when retained and no greater than k.
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
		if (dict.empty())
			return {};
		const std::vector<std::size_t> eligible = eligible_paths(query.size(), k);
		if (!eligible.empty())
			lifter.build(query, k);
		return score_paths(query, k, eligible);
	}

private:
	/** A reusable DP row at a prefix boundary of a sorted path. */
	struct QueryState {
		POAGraph::node_id node;
		std::vector<int> row;
		std::pair<std::size_t, std::size_t> pos_range;
	};

	/** Lexicographic order over the canonical path-operation tokens. */
	bool path_less(std::size_t lhs, std::size_t rhs) const {
		const auto &a = compressed_paths[lhs];
		const auto &b = compressed_paths[rhs];
		if (a.start != b.start)
			return a.start < b.start;
		const std::size_t common = std::min(a.steps.size(), b.steps.size());
		for (std::size_t i = 0; i < common; ++i) {
			const auto &x = a.steps[i];
			const auto &y = b.steps[i];
			if (x.type != y.type)
				return static_cast<unsigned char>(x.type) <
				       static_cast<unsigned char>(y.type);
			if (x.type == POAGraph::EdgeType::HEAVY) {
				if (x.length != y.length)
					return x.length < y.length;
			} else if (x.next != y.next) {
				return x.next < y.next;
			}
		}
		if (a.steps.size() != b.steps.size())
			return a.steps.size() < b.steps.size();
		return lhs < rhs;
	}

	/** Number of saved states shared by two canonical path prefixes. */
	static std::size_t common_prefix_states(
	    const POAGraph::CompressedPath &a,
	    const POAGraph::CompressedPath &b) {
		if (a.start != b.start)
			return 0;
		std::size_t common = 0;
		while (common < a.steps.size() && common < b.steps.size()) {
			const auto &x = a.steps[common];
			const auto &y = b.steps[common];
			if (x.type != y.type)
				break;
			if (x.type == POAGraph::EdgeType::HEAVY) {
				if (x.length != y.length)
					break;
			} else if (x.next != y.next) {
				break;
			}
			++common;
		}
		// State zero is the row after the starting node's label.
		return common + 1;
	}

	/** Filters by the length lower bound while preserving sorted path order. */
	std::vector<std::size_t> eligible_paths(std::size_t query_length,
	                                        int k) const {
		std::vector<std::size_t> eligible;
		eligible.reserve(sorted_path_ids.size());
		const long m = static_cast<long>(query_length);
		for (const std::size_t id : sorted_path_ids)
			if (std::abs(static_cast<long>(dict[id].size()) - m) <= k)
				eligible.push_back(id);
		return eligible;
	}

	/** Scores contiguous batches in parallel, preserving prefix sharing per batch. */
	std::vector<Result> score_paths(
	    const std::string &query, int k,
	    const std::vector<std::size_t> &eligible) const {
		std::vector<int> dist(dict.size(), DistMatrix::INF);
		const long m = static_cast<long>(query.size());
		for (std::size_t id = 0; id < dict.size(); ++id)
			if (dict[id].empty() &&
			    std::abs(static_cast<long>(dict[id].size()) - m) <= k) {
				const long empty_distance = m * cost.ins();
				dist[id] = static_cast<int>(
				    std::min(empty_distance,
				             static_cast<long>(DistMatrix::INF)));
			}
		if (!eligible.empty()) {
			const std::size_t thread_count = static_cast<std::size_t>(
			    num_threads == 0 ? omp_get_max_threads() : num_threads);
			const std::size_t batch_count =
			    thread_count > eligible.size() / 4
			        ? eligible.size()
			        : std::min(eligible.size(), thread_count * 4);
			const std::size_t batch_size = eligible.size() / batch_count;
			const std::size_t larger_batches = eligible.size() % batch_count;
			const std::size_t active_threads =
			    std::min(thread_count, batch_count);
#pragma omp parallel for schedule(dynamic) num_threads(active_threads)
			for (std::size_t batch = 0; batch < batch_count; ++batch) {
				const std::size_t begin =
				    batch * batch_size + std::min(batch, larger_batches);
				const std::size_t end =
				    begin + batch_size + (batch < larger_batches ? 1 : 0);
				score_batch(query, k, eligible, begin, end, dist);
			}
		}

		return collect_results(dist, k);
	}

	/** Scores one contiguous range, sharing DP states only within the batch. */
	void score_batch(const std::string &query, int k,
	                 const std::vector<std::size_t> &eligible,
	                 std::size_t begin, std::size_t end,
	                 std::vector<int> &dist) const {
		std::vector<QueryState> stack;
		const POAGraph::CompressedPath *previous = nullptr;

		for (std::size_t i = begin; i < end; ++i) {
			const std::size_t id = eligible[i];
			const auto &cp = compressed_paths[id];
			std::size_t common =
			    previous ? common_prefix_states(*previous, cp) : 0;
			// An earlier path may have stopped after its row became
			// impossible, so only states actually present can be reused.
			common = std::min(common, stack.size());
			while (stack.size() > common)
				stack.pop_back();

			if (stack.empty()) {
				std::vector<int> row = initial_row(cp.start, query, k);
				auto pos_range = lifter.window(cp.start);
				lifter.clip_row(row, pos_range);
				stack.push_back(
				    {cp.start, std::move(row), pos_range});
			}

			std::size_t step_index = stack.size() - 1;
			bool dead = stack.back().row.empty();
			while (!dead && step_index < cp.steps.size()) {
				const auto &step = cp.steps[step_index];
				const QueryState &state = stack.back();
				POAGraph::node_id cur = state.node;
				const auto full_range = lifter.window(cur);
				assert(state.pos_range.first >= full_range.first &&
				       state.pos_range.second <= full_range.second &&
				       state.pos_range.first <= state.pos_range.second);
				assert(state.row.size() ==
				       state.pos_range.second - state.pos_range.first);
				QueryState next_state = state;
				if (step.type == POAGraph::EdgeType::HEAVY) {
					if (step.length == 1) {
						const auto destination =
						    *graph.node(cur).heavy_neighbour;
						lifter.step(cur, destination, next_state.row,
						            next_state.pos_range);
						next_state.node = destination;
					} else {
						next_state.node = lifter.jump(
						    cur, step.length, next_state.row,
						    next_state.pos_range);
					}
				} else {
					lifter.step(cur, step.next, next_state.row,
					            next_state.pos_range);
					next_state.node = step.next;
				}
				stack.push_back(std::move(next_state));
				++step_index;
				dead = stack.back().row.empty();
			}
			if (!dead) {
				const QueryState &final_state = stack.back();
				const std::size_t query_end = query.size();
				if (query_end >= final_state.pos_range.first &&
				    query_end < final_state.pos_range.second)
					dist[id] = final_state.row[
					    query_end - final_state.pos_range.first];
			}
			previous = &cp;
		}
	}

	/** Returns result records in dictionary-id order. */
	std::vector<Result> collect_results(const std::vector<int> &dist,
	                                    int k) const {
		std::vector<Result> results;
		// Iterate dictionary ids in increasing order, independent of batch
		// completion order, so the public result order is deterministic.
		for (std::size_t id = 0; id < dict.size(); ++id)
			if (dist[id] <= k)
				results.push_back({id, dist[id]});
		return results;
	}

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
	/// Query-independent lifter: up table + decomposed-block set filled
	/// at construction; build(current_query, k) runs per query() call.
	mutable BinaryLifter lifter;
	/// Compressed (heavy-chain / light-step) representation of each dict path.
	std::vector<POAGraph::CompressedPath> compressed_paths;
	/// Dictionary ids in lexicographic compressed-path order.
	std::vector<std::size_t> sorted_path_ids;

	/** Builds the merged compressed representation of every dictionary path. */
	void build_compressed_paths() {
		compressed_paths.clear();
		compressed_paths.reserve(dict.size());
		for (std::size_t i = 0; i < dict.size(); ++i) {
			if (dict[i].empty())
				compressed_paths.emplace_back();
			else
				compressed_paths.push_back(graph.compressed_path(i));
		}
	}

	/** Sorts nonempty dictionary paths after decompose() rewrites their steps. */
	void sort_compressed_paths() {
		sorted_path_ids.clear();
		sorted_path_ids.reserve(dict.size());
		for (std::size_t i = 0; i < dict.size(); ++i)
			if (!dict[i].empty())
				sorted_path_ids.push_back(i);
		std::sort(sorted_path_ids.begin(), sorted_path_ids.end(),
		          [this](std::size_t a, std::size_t b) {
			          return path_less(a, b);
		          });
	}
};

} // namespace hlp_grep
