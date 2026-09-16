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
#include <memory>
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
	 *             defaults to the unit-cost model. The solver clones the
	 *             model, so it owns its own copy and there is no lifetime
	 *             requirement on @p cost (temporaries are fine).
	 */
	explicit Solver(std::vector<std::string> dict,
	                const CostModel &cost = DEFAULT_COST_MODEL)
	    : dict(std::move(dict)), cost(cost.clone()),
	      graph(this->dict, *this->cost) {
		build_compressed_paths();
	}

	/**
	 * @brief Copy constructor.
	 *
	 * The cost member is a unique_ptr holding a polymorphic model, so the
	 * compiler-generated copy machinery cannot be used; the model is
	 * deep-cloned and the graph is rebuilt against the cloned model.
	 */
	Solver(const Solver &other)
	    : dict(other.dict), cost(other.cost->clone()),
	      graph(this->dict, *this->cost) {
		build_compressed_paths();
	}

	Solver &operator=(const Solver &other) = delete;

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
		BinaryLifter lifter(graph, query, k);
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
			}
			const int dist = row.back();
			if (dist <= k)
				results.push_back({i, dist});
		}
		return results;
	}

private:
	/**
	 * @brief Initial DP row at the start node of a compressed path.
	 *
	 * Row entries are the node's window positions, i.e. the number of
	 * query characters already consumed arriving with the start node's
	 * base already handled: that base is matched against some consumed
	 * character or deleted, and every other consumed character is
	 * inserted. Same recurrence style as BinaryLifter's single-edge
	 * matrices. The window must match BinaryLifter::window (the start
	 * node's base is never crossed by an edge, so no matrix can produce
	 * this row).
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
		const long lo = std::max(0L, static_cast<long>(jmin) - mk);
		const long hi = std::min(m + 1, static_cast<long>(jmax) + mk + 2);

		std::vector<int> row;
		row.reserve(hi - lo);
		// Running prefix data over the whole query; windows only decide
		// which entries get recorded. The seed is INF (no consume seen
		// yet); entries stay true transition costs and saturate at INF.
		long best_alt =
		    DistMatrix::INF; // min over j <= i of (consume_j - ins_j)
		long sum_ins = 0;    // sum of insertion costs over query[0..i)
		for (long i = 0; i < hi; ++i) {
			if (i > 0) {
				sum_ins += cost->ins(query[i - 1]);
				best_alt = std::min(
				    best_alt,
				    static_cast<long>(
				        cost->consume(graph.base(start),
				                      query[i - 1])) -
				        static_cast<long>(cost->ins(query[i - 1])));
			}
			long cur = static_cast<long>(cost->del(graph.base(start))) +
			           sum_ins;
			if (i > 0)
				cur = std::min(cur, best_alt + sum_ins);
			if (cur > DistMatrix::INF)
				cur = DistMatrix::INF;
			if (i >= lo)
				row.push_back(static_cast<int>(cur));
		}
		return row;
	}

	std::vector<std::string> dict; ///< Dictionary of DNA sequences to search.
	/// Owned clone of the cost model; polymorphic storage (Solver is the
	/// only unique_ptr user) so it needs the manual deep-copy ctor above.
	std::unique_ptr<CostModel> cost;
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
