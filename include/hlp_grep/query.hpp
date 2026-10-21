/**
 * @file query.hpp
 * @brief Query text and threshold used by the edit-distance tables.
 */
#pragma once

#include <hlp_grep/string_hash.hpp>
#include <string>
#include <unordered_set>

namespace hlp_grep {

/** Immutable query data shared by edge-transition construction. */
struct Query {
	std::string sequence;
	int k;
	int max_level;
	StringHash hash;
	std::vector<std::unordered_set<StringHash::hash_type>> hash_sets;

	explicit Query(const std::string &sequence, int k)
	    : sequence(sequence), k(k), max_level(63 - __builtin_clzll(sequence.size())), hash(sequence) {
		hash_sets.resize(max_level + 1);
		for (int i = 0; i <= max_level; i++) {
			for (std::size_t j = 0; j + (1 << i) <= sequence.size(); j++)
				hash_sets[i].insert(hash.hash(j, 1 << i));
		}
	}
};

} // namespace hlp_grep
