/**
 * @file string_hash.hpp
 * @brief O(1) polynomial-hash queries for arbitrary string substrings.
 */
#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace hlp_grep {

/** Immutable string with O(1) hash queries for arbitrary substrings. */
struct StringHash {
	using hash_type = std::uint64_t;

	explicit StringHash(const std::string &text) : text_size(text.size()) {
		assert(text.size() <= MAX_LEN);
		prefix_hash.resize(text_size + 1, 0);
		hash_type z = 1;
		for (std::size_t i = 0; i < text_size; ++i) {
			prefix_hash[i + 1] = prefix_hash[i] + static_cast<unsigned char>(text[i]) * z;
			z = z * HASH_BASE;
		}

		static std::size_t max_len = 0;
		if (max_len < text_size) {
			for (std::size_t i = max_len; i < text_size; ++i) {
				inv_pow.push_back(inv_pow.back() * INV_HASH_BASE);
			}
			max_len = text_size;
		}
	}

	/** Hashes text[start..start + length) in O(1). */
	hash_type hash(std::size_t start, std::size_t length) const {
		assert(start <= text_size && length <= text_size - start);
		const std::size_t end = start + length;
		return (prefix_hash[end] - prefix_hash[start]) * inv_pow[start];
	}

	std::size_t size() const { return text_size; }

private:
	static constexpr hash_type HASH_BASE = 257ULL;
	static constexpr hash_type INV_HASH_BASE = 18374966859414961921ULL;

	std::size_t text_size;
	std::vector<hash_type> prefix_hash;
	inline static std::vector<hash_type> inv_pow{1};
};

} // namespace hlp_grep
