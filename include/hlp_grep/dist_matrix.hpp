/**
 * @file dist_matrix.hpp
 * @brief Distance matrices over the (min, +) semiring.
 *
 * DistMatrix stores a rectangular matrix of integer costs and provides the
 * min-plus matrix and matrix-vector products used to compose edit-distance
 * transition tables along heavy chains. No infinity
 * sentinel is built in; callers represent unreachable states with their own
 * large values and must keep them small enough that sums do not overflow.
 */
#pragma once
#include <hlp_grep/cost_model.hpp>

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace hlp_grep {

/**
 * @brief A rectangular matrix over the (min, +) semiring.
 *
 * Entry (i, j) is a cost. Storage is flat row-major.
 */
class DistMatrix {
public:
	/**
	 * @brief Constructs a rows x cols matrix filled with @p fill.
	 *
	 * @param rows Number of rows.
	 * @param cols Number of columns.
	 * @param fill Initial value for every entry (defaults to 0).
	 */
	DistMatrix(std::size_t rows, std::size_t cols, int fill = 0)
	    : rows(rows), cols(cols), data(rows * cols, fill) {}

	/// Number of rows.
	std::size_t num_rows() const { return rows; }

	/// Number of columns.
	std::size_t num_cols() const { return cols; }

	/**
	 * @brief Mutable access to entry (i, j).
	 *
	 * @param i Row index.
	 * @param j Column index.
	 */
	int &operator()(std::size_t i, std::size_t j) {
		return data[i * cols + j];
	}

	/**
	 * @brief Read-only access to entry (i, j).
	 *
	 * @param i Row index.
	 * @param j Column index.
	 */
	int operator()(std::size_t i, std::size_t j) const {
		return data[i * cols + j];
	}

	/**
	 * @brief Min-plus product of two matrices.
	 *
	 * Computes `out(i, j) = min_c (a(i, c) + b(c, j))` with plain integer
	 * addition, in O(a.rows * a.cols * b.cols) time.
	 *
	 * @param a Left factor.
	 * @param b Right factor.
	 * @return The min-plus product a (x) b.
	 * @throws std::invalid_argument if a.cols != b.rows or a.cols == 0
	 *         (an empty inner dimension leaves the product undefined).
	 */
	static DistMatrix min_plus_product(const DistMatrix &a,
	                                   const DistMatrix &b) {
		if (a.cols != b.rows)
			throw std::invalid_argument(
			    "DistMatrix::min_plus_product: dimension mismatch");
		if (a.cols == 0)
			throw std::invalid_argument(
			    "DistMatrix::min_plus_product: empty inner dimension");

		DistMatrix out(a.rows, b.cols);
		for (std::size_t i = 0; i < a.rows; ++i) {
			const int *arow = &a.data[i * a.cols];
			int *orow = &out.data[i * b.cols];
			for (std::size_t j = 0; j < b.cols; ++j)
				orow[j] = arow[0] + b.data[j];
			for (std::size_t c = 1; c < a.cols; ++c) {
				const int left = arow[c];
				const int *brow = &b.data[c * b.cols];
				for (std::size_t j = 0; j < b.cols; ++j) {
					const int cand = left + brow[j];
					if (cand < orow[j])
						orow[j] = cand;
				}
			}
		}
		return out;
	}

	/**
	 * @brief Min-plus product of a row vector with this matrix.
	 *
	 * Computes `out(j) = min_i (vec(i) + (*this)(i, j))` with plain integer
	 * addition.
	 *
	 * @param vec Row vector of costs, one per row of this matrix.
	 * @return The resulting row vector, one entry per column.
	 * @throws std::invalid_argument if vec.size() != rows or rows == 0
	 *         (an empty row vector leaves the result undefined).
	 */
	std::vector<int> min_plus_apply(const std::vector<int> &vec) const {
		if (vec.size() != rows)
			throw std::invalid_argument(
			    "DistMatrix::min_plus_apply: dimension mismatch");
		if (rows == 0)
			throw std::invalid_argument(
			    "DistMatrix::min_plus_apply: empty row vector");

		std::vector<int> out(cols);
		for (std::size_t j = 0; j < cols; ++j)
			out[j] = vec[0] + data[j];
		for (std::size_t i = 1; i < rows; ++i) {
			const int left = vec[i];
			const int *row = &data[i * cols];
			for (std::size_t j = 0; j < cols; ++j) {
				const int cand = left + row[j];
				if (cand < out[j])
					out[j] = cand;
			}
		}
		return out;
	}

private:
	std::size_t rows;      ///< Number of rows.
	std::size_t cols;      ///< Number of columns.
	std::vector<int> data; ///< Flat row-major storage.
};

} // namespace hlp_grep
