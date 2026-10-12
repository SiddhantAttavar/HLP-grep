/**
 * @file query.hpp
 * @brief Query text, threshold, and rolling hashes of power-of-two substrings.
 */
#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace hlp_grep {

/** Immutable query data shared by edge-transition construction. */
struct Query {
	static constexpr std::uint64_t HASH_BASE = 256;
	static constexpr std::uint64_t HASH_MOD = 1000000007ULL;
	static constexpr std::size_t HASH_POWER_LEVELS =
	    std::numeric_limits<std::size_t>::digits;

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
		if (levels == 0)
			return;

		hashes[0].reserve(n);
		for (unsigned char c : sequence)
			hashes[0].push_back(c);

		for (std::size_t level = 1; level < levels; ++level) {
			const std::size_t half = std::size_t{1} << (level - 1);
			hashes[level].resize(n - (std::size_t{1} << level) + 1);
			for (std::size_t i = 0; i < hashes[level].size(); ++i)
				hashes[level][i] = (hashes[level - 1][i] +
				    hashes[level - 1][i + half] * power_of_base(half)) % HASH_MOD;
		}
	}

	/** Returns HASH_BASE^exponent mod HASH_MOD for a power-of-two exponent. */
	static std::uint64_t power_of_base(std::size_t exponent) {
		assert(exponent > 0 && (exponent & (exponent - 1)) == 0);
		static const auto powers = [] {
			std::array<std::uint64_t, HASH_POWER_LEVELS> result{};
			result[0] = HASH_BASE;
			for (std::size_t level = 1; level < result.size(); ++level)
				result[level] = (result[level - 1] * result[level - 1]) % HASH_MOD;
			return result;
		}();
		const std::size_t level =
		    static_cast<std::size_t>(__builtin_ctzll(exponent));
		assert(level < powers.size());
		return powers[level];
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
