/**
 * @file hlp_grep.hpp
 * @brief Edit distance search queries on a DNA sequence dictionary using
 *        heavy-light decomposition on pangenome paths.
 *
 * HLP-grep indexes a dictionary of DNA sequences (represented as paths in a
 * pangenome) and answers edit distance search queries: given a query string
 * and a threshold k, it finds all dictionary sequences within edit distance
 * k of the query string.
 */
#pragma once
#include <string>
#include <utility>
#include <vector>

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
	    : ins(ins), del(del), alphabet_(std::move(alphabet)),
	      matrix_(std::move(matrix)) {
		if (matrix_.empty())
			matrix_ = unit_matrix(alphabet_.size());
		for (std::size_t i = 0; i < alphabet_.size(); ++i)
			index_[static_cast<unsigned char>(alphabet_[i])] = i;
	}

	int ins; ///< Cost of inserting a character.
	int del; ///< Cost of deleting a character.

	/**
	 * @brief Returns the cost of matching two characters.
	 *
	 * Looks the characters up in the substitution matrix provided at
	 * construction time.
	 *
	 * @param a First character.
	 * @param b Second character.
	 * @return The cost of matching @p a with @p b (0 if they match).
	 */
	int match(char a, char b) const {
		return matrix_[index_of(a)][index_of(b)];
	}

private:
	std::string alphabet_;                ///< Alphabet, defining matrix order.
	std::vector<std::vector<int>> matrix_; ///< Substitution cost matrix.
	std::size_t index_[256] = {};          ///< char -> alphabet index.

	static std::vector<std::vector<int>> unit_matrix(std::size_t sigma) {
		std::vector<std::vector<int>> m(sigma, std::vector<int>(sigma, 1));
		for (std::size_t i = 0; i < sigma; ++i)
			m[i][i] = 0;
		return m;
	}

	std::size_t index_of(char c) const {
		return index_[static_cast<unsigned char>(c)];
	}
};

/**
 * @brief A single approximate match returned by Solver::query().
 */
struct Result {
	std::size_t id;   ///< Index of the matched sequence in the dictionary (as passed to the constructor).
	int dist; ///< Edit distance between the matched sequence and the query string.
};

/**
 * @brief Index and query engine for edit distance search over a sequence dictionary.
 *
 * Usage: construct a concrete Solver implementation with the dictionary, then
 * call query() to retrieve all dictionary sequences within a given edit
 * distance threshold of a query string.
 */
class Solver {
public:
	/**
	 * @brief Constructs the solver over the given dictionary with the given
	 *        edit cost model.
	 *
	 * @param dict Dictionary of DNA sequences to search. The position of each
	 *             sequence in this vector defines the `id` reported in Result.
	 * @param cost Cost model defining the costs of the basic edit operations.
	 *             Defaults to the unit-cost model (`CostModel{}`).
	 */
	explicit Solver(std::vector<std::string> dict, CostModel cost = {})
	    : dict(std::move(dict)), cost(cost) {}

	virtual ~Solver() = default;

	/**
	 * @brief Finds all dictionary sequences within edit distance k of the query.
	 *
	 * @param query The query DNA sequence.
	 * @param k     Maximum allowed edit distance (threshold).
	 * @return A vector of Result, one per matched dictionary sequence,
	 *         each containing the sequence's index in the dictionary and its
	 *         edit distance to the query. Results are sorted by `id` in
	 *         ascending order.
	 */
	virtual std::vector<Result> query(const std::string &query, int k) = 0;

protected:
	std::vector<std::string> dict; ///< Dictionary of DNA sequences to search.
	CostModel cost; ///< Cost model used for the edit distance computations.
};
