/**
 * @file query.hpp
 * @brief Query text, threshold, and rolling hashes of power-of-two substrings.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace hlp_grep {

/** Immutable query data shared by edge-transition construction. */
struct Query {
	static constexpr std::uint64_t HASH_BASE = 256;
	static constexpr std::uint64_t HASH_MOD = 1000000007ULL;

	std::string sequence;
	int k;
	std::vector<std::vector<std::uint64_t>> hashes;

	explicit Query(const std::string &sequence, int k)
	    : sequence(sequence), k(k) {
		const std::size_t n = sequence.size();
		std::size_t levels = 0;
		while (levels < 63 && (std::size_t{1} << levels) <= n)
			++levels;
		hashes.resize(levels);
		std::vector<std::uint64_t> powers(levels, 1);
		if (levels == 0)
			return;

		powers[0] = HASH_BASE;
		hashes[0].reserve(n);
		for (unsigned char c : sequence)
			hashes[0].push_back(c);

		for (std::size_t level = 1; level < levels; ++level) {
			const std::size_t half = std::size_t{1} << (level - 1);
			powers[level] = (powers[level - 1] * powers[level - 1]) % HASH_MOD;
			hashes[level].resize(n - (std::size_t{1} << level) + 1);
			for (std::size_t i = 0; i < hashes[level].size(); ++i)
				hashes[level][i] = (hashes[level - 1][i] +
				    hashes[level - 1][i + half] * powers[level - 1]) % HASH_MOD;
		}
	}

	std::uint64_t hash(const std::string &text) const {
		std::uint64_t result = 0, power = 1;
		for (unsigned char c : text) {
			result = (result + static_cast<std::uint64_t>(c) * power) % HASH_MOD;
			power = (power * HASH_BASE) % HASH_MOD;
		}
		return result;
	}
};

} // namespace hlp_grep
