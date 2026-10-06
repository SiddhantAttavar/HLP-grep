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
#include <cstdlib>
#include <iostream>
#include <stdexcept>
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
	 *             Every sequence must be non-empty; an empty sequence is
	 *             rejected with std::runtime_error by POAGraph.
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
		initialize_paths();
	}

	/**
	 * @brief Loads a solver from a previously saved POAGraph file.
	 *
	 * Skips the dictionary build step: the graph (nodes, edges, paths,
	 * topological order, heavy-chain metadata) is read from @p graph_file
	 * instead of being rebuilt from dictionary sequences. The supplied
	 * cost model must match the one stored in the file. Solver query
	 * helpers are reconstructed from the loaded graph, while
	 * query-specific tables are still built per query() call.
	 *
	 * @param graph_file Path to a file written by save_graph().
	 * @param cost Cost model the graph was saved with; must outlive the
	 *             solver.
	 * @param num_threads Worker threads, as in the dictionary constructor.
	 */
	static Solver load_graph(const std::string &graph_file,
	                         const CostModel &cost = DEFAULT_COST_MODEL,
	                         std::size_t num_threads = 0) {
		return Solver(POAGraph::load_file(graph_file, cost), cost,
		              num_threads);
	}

	/** Saves the complete graph structure used by this solver. */
	void save_graph(const std::string &filename) const {
		graph.save_file(filename);
	}

	/**
	 * @brief Finds all dictionary sequences within edit distance k of the
	 *        query.
	 *
	 * Rebuilds the per-query chain tables of the pre-marked BinaryLifter, then
	 * evaluates each length-eligible path with bounded A*. Integer f-score
	 * buckets 0..k are scanned in order. The heuristic is the unavoidable
	 * net-indel cost: min(ins, del) times the distance from the current to the
	 * terminal diagonal. Each search offset indexes a
	 * sparse path coordinate whose actual sequence position is built from
	 * CompressedEdge::label_shift; query_pos is that position plus the diagonal.
	 * Character transitions advance through retained label positions; a heavy
	 * edge uses its matrix to skip every destination label in its compressed
	 * block. A path's distance is finalized when its terminal state is removed
	 * from its f-score bucket.
	 *
	 * Scoring ignores relaxations above k. All operation costs are
	 * nonnegative, and insertion/deletion costs are at least 1, so an
	 * alignment within k stays in the represented diagonal band.
	 *
	 * Sequences whose length differs from the query by more than k are
	 * skipped without scoring: their length difference is a lower bound
	 * on the edit distance, and their DP states need not fit the windows
	 * the lifter precomputes and the diagonal band remains sufficient.
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
	/** Builds solver state around an already-built graph (load path). */
	explicit Solver(POAGraph graph, const CostModel &cost,
	                std::size_t num_threads)
	    : cost(cost), num_threads(num_threads), graph(std::move(graph)),
	      lifter(this->graph, num_threads) {
		initialize_paths();
	}

	void initialize_paths() {
		if (dict.empty() && graph.num_sequences() != 0) {
			dict.reserve(graph.num_sequences());
			for (std::size_t i = 0; i < graph.num_sequences(); ++i) {
				std::string sequence;
				for (POAGraph::node_id u : graph.path(i))
					sequence += graph.seq(u);
				dict.push_back(std::move(sequence));
			}
		}
		build_compressed_paths();
		lifter.decompose(compressed_paths);
		build_search_paths();
	}

	// -- implementation ----------------------------------------------------
	/** Filters paths by the length-difference lower bound. */
	std::vector<std::size_t> eligible_paths(std::size_t query_length,
	                                        int k) const {
		std::vector<std::size_t> eligible;
		eligible.reserve(dict.size());
		const long m = static_cast<long>(query_length);
		for (std::size_t id = 0; id < dict.size(); ++id)
			if (std::abs(static_cast<long>(dict[id].size()) - m) <= k)
				eligible.push_back(id);
		return eligible;
	}

	/** Scores independent compressed paths in parallel using bounded A*. */
	std::vector<Result> score_paths(
	    const std::string &query, int k,
	    const std::vector<std::size_t> &eligible) const {
		std::vector<int> dist(dict.size(), DistMatrix::INF);
		if (!eligible.empty()) {
			const std::size_t thread_count = static_cast<std::size_t>(
			    num_threads == 0 ? std::max(omp_get_max_threads() / 2, 1)
				: num_threads);
			const std::size_t active_threads =
			    std::min(thread_count, eligible.size());
#pragma omp parallel for schedule(dynamic) num_threads(active_threads)
			for (std::size_t i = 0; i < eligible.size(); ++i)
				dist[eligible[i]] = score_path(query, k, eligible[i]);
		}

		return collect_results(dist, k);
	}

	/** One queued state at a sparse path-coordinate index and diagonal. */
	struct SearchState {
		std::size_t offset;
		std::size_t diagonal_index;
		int distance;
	};

	/** Sparse per-path coordinates used by the distance-bucket search. */
	struct PathSearchIndex {
		std::vector<std::size_t> step_len;
		std::vector<std::size_t> label_pos;
		std::vector<std::size_t> seq_pos;
		std::vector<POAGraph::node_id> node;
	};

	/** Computes one path's exact distance with buckets indexed by f-score. */
	int score_path(const std::string &query, int k, std::size_t id) const {
		const auto &path_index = search_paths[id];
		const std::size_t sequence_length = path_index.seq_pos.back();
		const long diag_delta = static_cast<long>(query.size()) - 
			static_cast<long>(sequence_length);
		if (std::abs(diag_delta) > k) {
			return DistMatrix::INF;
		}

		const long max_shift = (k - std::abs(diag_delta)) / 2;
		const long min_diag = std::min(0L, diag_delta) - max_shift;
		const long max_diag = std::max(0L, diag_delta) + max_shift;

		const std::size_t state_count = path_index.label_pos.size();
		const std::size_t width = max_diag - min_diag + 1;
		std::vector<int> best(state_count * width, DistMatrix::INF);

		// A net diagonal change requires at least this many indels. Using the
		// minimum indel cost makes the lower bound consistent for either sign.
		std::vector<std::vector<SearchState>> buckets(
		    static_cast<std::size_t>(k) + 1);

		auto relax = [&](std::size_t offset, int diagonal, int candidate) {
			const std::size_t query_pos = diagonal + path_index.seq_pos[offset];
			if (query_pos < 0 || query_pos > query.size())
				return;
			const int delta = std::abs(diagonal - static_cast<int>(diag_delta));
			const int f_score = candidate + delta;
			if (f_score > k)
				return;
			const std::size_t d = diagonal - min_diag;
			int &known = best[offset * width + d];
			if (candidate < known) {
				known = candidate;
				buckets[static_cast<std::size_t>(f_score)].push_back(
				    {offset, d, candidate});
			}
		};
		relax(0, 0, 0);

		for (int f_score = 0; f_score <= k; ++f_score) {
			auto &bucket = buckets[static_cast<std::size_t>(f_score)];
			for (std::size_t cursor = 0; cursor < bucket.size(); ++cursor) {
				const auto [offset, diagonal_index, distance] = bucket[cursor];
				if (best[offset * width + diagonal_index] != distance)
					continue;

				const POAGraph::node_id node = path_index.node[offset];
				const std::size_t seq_pos = path_index.seq_pos[offset];
				const std::size_t label_pos = path_index.label_pos[offset];
				const std::size_t step_len = path_index.step_len[offset];

				const int diagonal = static_cast<int>(diagonal_index + min_diag);
				const std::size_t query_pos = diagonal + seq_pos;

				if (query_pos == query.size() && offset == state_count - 1)
					return distance;

				// Insert a query character without advancing in the dictionary.
				relax(offset, diagonal + 1, distance + cost.ins());

				if (offset == state_count - 1)
					continue;

				const std::size_t next_offset = offset + 1;
				const std::size_t next_seq_pos = path_index.seq_pos[next_offset];
				if (step_len == 0) {
					relax(next_offset, diagonal - 1, distance + cost.del());
					if (query_pos < query.size()) {
						const int consume_cost =
							cost.consume(graph.seq(node)[label_pos], query[query_pos]);
						relax(next_offset, diagonal, distance + consume_cost);
					}
					continue;
				}

				if (diagonal > min_diag && distance - cost.ins() ==
					best[offset * width + diagonal_index - 1]) {
					continue;
				}

				const std::size_t label_shift = next_seq_pos - seq_pos;
				if (diagonal < max_diag - 1 && distance - cost.ins() ==
					best[offset * width + diagonal_index + 1]) {
					relax(next_offset, diagonal - static_cast<int>(label_shift),
						static_cast<int>(distance + label_shift));
					continue;
				}

				// lifter.jump_point(
				// 	node, step_len, query_pos, distance,
				// 	diag_delta - diagonal, label_shift, row, range
				// );

				const DistMatrix &mat = lifter.jump_mat(node, step_len);

				const auto [col_lo, col_hi] = mat.col_range();
				const std::size_t target_pos = query_pos + label_shift;
				const int curr_delta = static_cast<int>(diag_delta) - diagonal;
				const std::size_t shifted_target_pos = 
					std::max(0L, static_cast<long>(target_pos + curr_delta));
				const long max_shift = (k - distance - std::abs(curr_delta)) / 2;
				const std::size_t lo = std::max(static_cast<long>(col_lo),
					static_cast<long>(std::min(shifted_target_pos, target_pos)) -
					max_shift);
				const std::size_t hi = std::min(col_hi,
					std::max(shifted_target_pos, target_pos) + max_shift + 1);

				bool change = false;
				const std::size_t matrix_row = query_pos - mat.row_range().first;
				int output_diagonal = static_cast<int>(lo) - 
					static_cast<int>(next_seq_pos);
				for (std::size_t b = lo; b < hi; ++b) {
					int next_distance = distance + mat(matrix_row, b - col_lo);
					const int delta = std::abs(output_diagonal -
						 static_cast<int>(diag_delta));
					const int f_score = next_distance + delta;
					if (f_score <= k) {
						const std::size_t d = output_diagonal - min_diag;
						int &known = best[next_offset * width + d];
						if (next_distance < known) {
							known = next_distance;
							buckets[static_cast<std::size_t>(f_score)].push_back(
								{next_offset, d, next_distance});
							change = true;
						}
						else if (change) {
							break;
						}
					}
					output_diagonal++;
				}
			}
		}
		return DistMatrix::INF;
	}

	/** Returns result records in dictionary-id order. */
	std::vector<Result> collect_results(const std::vector<int> &dist,
	                                    int k) const {
		std::vector<Result> results;
		// Iterate dictionary ids in increasing order, independent of parallel
		// path completion order, so the public result order is deterministic.
		for (std::size_t id = 0; id < dict.size(); ++id)
			if (dist[id] <= k)
				results.push_back({id, dist[id]});
		return results;
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
	/// Sparse character and heavy-landing coordinates for each dictionary path.
	std::vector<PathSearchIndex> search_paths;
	/** Builds the merged compressed representation of every dictionary path. */
	void build_compressed_paths() {
		compressed_paths.clear();
		compressed_paths.reserve(dict.size());
		for (std::size_t i = 0; i < dict.size(); ++i)
			compressed_paths.push_back(graph.compressed_path(i));
	}

	/** Builds sparse coordinates, omitting all interior positions of heavy blocks. */
	void build_search_paths() {
		search_paths.clear();
		search_paths.resize(compressed_paths.size());
		for (std::size_t i = 0; i < compressed_paths.size(); ++i) {
			const auto &cp = compressed_paths[i];
			auto &index = search_paths[i];

			POAGraph::node_id current_node = cp.start;
			std::size_t sequence_offset = 0;

			auto append_light = [&]() {
				const std::size_t label_size = graph.seq(current_node).size();
				for (std::size_t label_pos = 0; label_pos < label_size;
					 ++label_pos) {
					index.step_len.push_back(0);
					index.label_pos.push_back(label_pos);
					index.seq_pos.push_back(sequence_offset);
					index.node.push_back(current_node);
					sequence_offset++;
				}
			};

			append_light();
			for (std::size_t step_index = 0;
			     step_index < cp.steps.size(); ++step_index) {
				const auto &step = cp.steps[step_index];
				if (step.type == POAGraph::EdgeType::LIGHT) {
					current_node = step.next;
					append_light();
				}
				else if (step.label_shift == 1) {
					current_node = *graph.node(current_node).heavy_neighbour;
					append_light();
				}
				else {
					index.step_len.push_back(step.length);
					index.label_pos.push_back(0);
					index.seq_pos.push_back(sequence_offset);
					index.node.push_back(current_node);
					sequence_offset += step.label_shift;
					current_node = lifter.jump_node(current_node, step.length);
				}
			}
			index.step_len.push_back(0);
			index.label_pos.push_back(0);
			index.seq_pos.push_back(sequence_offset);
			index.node.push_back(current_node);
		}
	}
};

} // namespace hlp_grep
