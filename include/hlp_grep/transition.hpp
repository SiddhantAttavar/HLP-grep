/**
 * @file transition.hpp
 * @brief Hash and mismatch summary for a power-of-two heavy-chain block.
 */
#pragma once

#include <hlp_grep/poa_graph.hpp>
#include <hlp_grep/query.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>

namespace hlp_grep {

struct Transition {
	std::uint64_t label_hash = 0;
	std::size_t label_len = 0;
	int mismatches = 0;
	Transition() = default;

	/** Summary of one heavy edge, from u to v. */
	Transition(POAGraph::node_id, POAGraph::node_id v,
	           const POAGraph &graph, const Query &query) {
		label_hash = query.hash(std::string(1, graph.seq(v)));
		label_len = 1;
	}

	/** Composes consecutive power-of-two heavy blocks and applies the filter. */
	static Transition compose(const Transition &left, const Transition &right,
	                          const POAGraph &graph, POAGraph::node_id start,
	                          const Query &query) {
		Transition result;
		result.label_len = left.label_len + right.label_len;
		std::uint64_t power = 1;
		for (std::size_t i = 0; i < left.label_len; ++i)
			power = (power * Query::HASH_BASE) % Query::HASH_MOD;
		result.label_hash =
		    (left.label_hash + right.label_hash * power) % Query::HASH_MOD;

		const auto [jmin, jmax] = graph.pos_range(start);
		const std::size_t k = static_cast<std::size_t>(query.k);
		const std::size_t lo = jmin + 1 > k ? jmin + 1 - k : 0;
		const std::size_t hi = std::min(query.sequence.size() + 1,
		                                jmax + 1 + k + 1);
		bool found = false;
		if (result.label_len <= query.sequence.size()) {
			const std::size_t level =
			    63 - static_cast<std::size_t>(__builtin_clzll(result.label_len));
			const std::size_t max_start =
				query.sequence.size() - result.label_len + 1;
			if (level < query.hashes.size())
				for (std::size_t pos = lo; pos < std::min(hi, max_start); ++pos)
					if (query.hashes[level][pos] == result.label_hash) {
						found = true;
						break;
					}
		}
		result.mismatches = found ? 0 : std::max(
		    1, left.mismatches + right.mismatches);
		return result;
	}
};

} // namespace hlp_grep
