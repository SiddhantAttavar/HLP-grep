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
	StringHash hash;

	explicit Query(const std::string &sequence, int k)
	    : sequence(sequence), k(k), hash(sequence) {}
};

} // namespace hlp_grep
