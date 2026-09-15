/**
 * @file cost_model.hpp
 * @brief Cost model for edit distance computations.
 */
#pragma once
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace hlp_grep {

/**
 * @brief Cost model for edit distance computations.
 *
 * Defines the costs of the basic edit operations.
 */
struct CostModel {
	/**
	 * @brief Constructs a cost model with the given insertion/deletion costs
	 *        and substitution matrix.
	 *
	 * @param ins      Cost of inserting a character.
	 * @param del      Cost of deleting a character.
	 * @param alphabet Characters of the alphabet, defining the row/column
	 *                 order of @p matrix.
	 * @param matrix   Substitution cost matrix: `matrix[i][j]` is the cost of
	 *                 matching `alphabet[i]` with `alphabet[j]`. Characters not
	 *                 present in @p alphabet are looked up at index 0.
	 *                 Defaults to the unit-cost model over the DNA alphabet
	 *                 (0 on the diagonal, 1 elsewhere).
	 */
	CostModel(int ins = 1, int del = 1, std::string alphabet = "AGCT",
	          std::vector<std::vector<int>> matrix = {})
	    : ins(ins), del(del) {
		this->alphabet = std::move(alphabet);
		if (matrix.empty())
			this->matrix = unit_matrix(this->alphabet.size());
		else
			this->matrix = std::move(matrix);
		for (std::size_t i = 0; i < this->alphabet.size(); ++i)
			index[static_cast<unsigned char>(this->alphabet[i])] = i;
	}

	int ins; ///< Cost of inserting a character.
	int del; ///< Cost of deleting a character.

	/**
	 * @brief Returns the cost of matching two characters.
	 *
	 * @param a First character.
	 * @param b Second character.
	 * @return The cost of matching @p a with @p b (0 if they match).
	 */
	int match(char a, char b) const {
		return matrix[index_of(a)][index_of(b)];
	}

private:
	std::string alphabet;                 ///< Alphabet, defining matrix order.
	std::vector<std::vector<int>> matrix; ///< Substitution cost matrix.
	std::size_t index[256] = {};          ///< char -> alphabet index.

	static std::vector<std::vector<int>> unit_matrix(std::size_t sigma) {
		std::vector<std::vector<int>> m(sigma, std::vector<int>(sigma, 1));
		for (std::size_t i = 0; i < sigma; ++i)
			m[i][i] = 0;
		return m;
	}

	std::size_t index_of(char c) const {
		return index[static_cast<unsigned char>(c)];
	}
};

} // namespace hlp_grep
