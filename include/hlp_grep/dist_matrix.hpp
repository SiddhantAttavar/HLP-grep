/**
 * @file dist_matrix.hpp
 * @brief Distance matrices over the (min, +) semiring.
 *
 * DistMatrix stores a rectangular matrix of integer costs and provides the
 * min-plus matrix and matrix-vector products used to compose edit-distance
 * transition tables along heavy chains. DistMatrix::INF is the infinity
 * sentinel for unreachable states; every product saturates candidates at
 * INF, so entries never exceed INF and plain integer sums (at most
 * 2 * INF = INT_MAX - 1) cannot overflow.
 */
#pragma once
#include <hlp_grep/cost_model.hpp>

#include <cstddef>
#include <limits>
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
	/// Infinity sentinel for unreachable states; saturated by every
	/// product, so entries never exceed it.
	static constexpr int INF = std::numeric_limits<int>::max() / 2;
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
	 * addition, in O(a.rows * a.cols * b.cols) time. Candidates saturate
	 * at INF, so entries never exceed INF.
	 *
	 * @param a Left factor; entries must not exceed INF.
	 * @param b Right factor; entries must not exceed INF.
	 * @param monge When true, both factors are Monge: for each row of a,
	 *        the column argmins of arow(c) + b(c, j) are monotone in j,
	 *        so each output row is found with divide-and-conquer in
	 *        O(a.cols + b.cols) argmin steps instead of the cubic scan.
	 *        The caller guarantees Monge-ness; no verification is done.
	 * @return The min-plus product a (x) b.
	 * @throws std::invalid_argument if a.cols != b.rows or a.cols == 0
	 *         (an empty inner dimension leaves the product undefined).
	 */
	static DistMatrix min_plus_product(const DistMatrix &a,
	                                   const DistMatrix &b,
	                                   bool monge = false) {
		if (a.cols != b.rows)
			throw std::invalid_argument(
			    "DistMatrix::min_plus_product: dimension mismatch");
		if (a.cols == 0)
			throw std::invalid_argument(
			    "DistMatrix::min_plus_product: empty inner dimension");

		DistMatrix out(a.rows, b.cols);
		if (a.rows == 0 || b.cols == 0)
			return out;
		if (monge) {
			for (std::size_t i = 0; i < a.rows; ++i) {
				const int *arow = &a.data[i * a.cols];
				int *orow = &out.data[i * b.cols];
				// An all-INF input row saturates every candidate at INF:
				// skip the divide-and-conquer and emit the row directly.
				bool all_inf = true;
				for (std::size_t c = 0; c < a.cols; ++c)
					if (arow[c] < INF) {
						all_inf = false;
						break;
					}
				if (all_inf) {
					for (std::size_t j = 0; j < b.cols; ++j)
						orow[j] = INF;
					continue;
				}
				monge_row(arow, b, orow, 0, b.cols - 1, 0, a.cols - 1);
			}
			return out;
		}
		for (std::size_t i = 0; i < a.rows; ++i) {
			const int *arow = &a.data[i * a.cols];
			int *orow = &out.data[i * b.cols];
			for (std::size_t j = 0; j < b.cols; ++j) {
				const int cand = arow[0] + b.data[j];
				orow[j] = cand > INF ? INF : cand;
			}
			for (std::size_t c = 1; c < a.cols; ++c) {
				const int left = arow[c];
				const int *brow = &b.data[c * b.cols];
				for (std::size_t j = 0; j < b.cols; ++j) {
					int cand = left + brow[j];
					if (cand > INF)
						cand = INF;
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
	 * addition. Candidates saturate at INF, so entries never exceed INF.
	 *
	 * @param vec Row vector of costs, one per row of this matrix; entries
	 *            must not exceed INF.
	 * @param monge When true, this matrix is Monge: the column argmins
	 *        of vec(i) + (*this)(i, j) are monotone in j, so the result
	 *        is found with divide-and-conquer in O(rows + cols) argmin
	 *        steps instead of the full rows x cols scan. The caller
	 *        guarantees Monge-ness; no verification is done.
	 * @return The resulting row vector, one entry per column.
	 * @throws std::invalid_argument if vec.size() != rows or rows == 0
	 *         (an empty row vector leaves the result undefined).
	 */
	std::vector<int> min_plus_apply(const std::vector<int> &vec,
	                                bool monge = false) const {
		if (vec.size() != rows)
			throw std::invalid_argument(
			    "DistMatrix::min_plus_apply: dimension mismatch");
		if (rows == 0)
			throw std::invalid_argument(
			    "DistMatrix::min_plus_apply: empty row vector");

		std::vector<int> out(cols);
		if (monge) {
			if (cols == 0)
				return out;
			monge_row(&vec[0], *this, &out[0], 0, cols - 1, 0,
			          rows - 1);
			return out;
		}
		for (std::size_t j = 0; j < cols; ++j) {
			const int cand = vec[0] + data[j];
			out[j] = cand > INF ? INF : cand;
		}
		for (std::size_t i = 1; i < rows; ++i) {
			const int left = vec[i];
			const int *row = &data[i * cols];
			for (std::size_t j = 0; j < cols; ++j) {
				int cand = left + row[j];
				if (cand > INF)
					cand = INF;
				if (cand < out[j])
					out[j] = cand;
			}
		}
		return out;
	}

private:
	/**
	 * @brief Divide-and-conquer column minima of M(c, j) = arow(c) + b(c, j).
	 *
	 * With b Monge, M is Monge and the column argmins are monotone in j,
	 * so the argmin of the middle column bounds the argmin ranges of both
	 * halves. Ties resolve to the leftmost argmin, preserving monotonicity.
	 *
	 * @param arow One row of the left factor, length b.rows.
	 * @param b Right factor.
	 * @param orow Output row to fill, length b.cols.
	 * @param j_lo First column of the range to solve.
	 * @param j_hi Last column of the range to solve (inclusive).
	 * @param k_lo First candidate row of the argmin range.
	 * @param k_hi Last candidate row of the argmin range (inclusive).
	 */
	static void monge_row(const int *arow, const DistMatrix &b, int *orow,
	                      std::size_t j_lo, std::size_t j_hi,
	                      std::size_t k_lo, std::size_t k_hi) {
		if (j_lo > j_hi)
			return;
		const std::size_t j_mid = (j_lo + j_hi) / 2;
		std::size_t best_k = k_lo;
		int best = arow[k_lo] + b.data[k_lo * b.cols + j_mid];
		if (best > INF)
			best = INF;
		for (std::size_t k = k_lo + 1; k <= k_hi; ++k) {
			int cand = arow[k] + b.data[k * b.cols + j_mid];
			if (cand > INF)
				cand = INF;
			if (cand < best) {
				best = cand;
				best_k = k;
			}
		}
		orow[j_mid] = best;
		if (j_mid > j_lo)
			monge_row(arow, b, orow, j_lo, j_mid - 1, k_lo, best_k);
		if (j_mid < j_hi)
			monge_row(arow, b, orow, j_mid + 1, j_hi, best_k, k_hi);
	}

	std::size_t rows;      ///< Number of rows.
	std::size_t cols;      ///< Number of columns.
	std::vector<int> data; ///< Flat row-major storage.
};

} // namespace hlp_grep
