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

#include <cassert>
#include <cstddef>
#include <limits>
#include <utility>
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
	    : DistMatrix(rows, cols, fill, {0, rows}, {0, cols}) {}

	/**
	 * @brief Constructs a matrix with explicit query-position ranges.
	 *
	 * @param rows Number of rows.
	 * @param cols Number of columns.
	 * @param fill Initial value for every entry.
	 * @param row_range Half-open query-position range represented by rows.
	 * @param col_range Half-open query-position range represented by columns.
	 */
	DistMatrix(std::size_t rows, std::size_t cols, int fill,
	           std::pair<std::size_t, std::size_t> row_range,
	           std::pair<std::size_t, std::size_t> col_range)
	    : rows(rows), cols(cols), row_range_value(row_range),
	      col_range_value(col_range),
	      data(rows * cols, fill) {
		assert(row_range.second >= row_range.first &&
		       row_range.second - row_range.first == rows);
		assert(col_range.second >= col_range.first &&
		       col_range.second - col_range.first == cols);
	}

	/// Number of rows.
	std::size_t num_rows() const { return rows; }

	/// Number of columns.
	std::size_t num_cols() const { return cols; }

	/// Half-open query-position range represented by matrix rows.
	std::pair<std::size_t, std::size_t> row_range() const {
		return row_range_value;
	}

	/// Half-open query-position range represented by matrix columns.
	std::pair<std::size_t, std::size_t> col_range() const {
		return col_range_value;
	}

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
	 * addition. Candidates saturate at INF, so entries never exceed INF.
	 *
	 * Both factors must be layers of edit-distance DAGs driven by a
	 * CostModel-validated cost assignment (match == 0, nonnegative
	 * ins/del, ins + del >= mismatch): mat(a, b) is the lowest ED cost to
	 * go from (row1, a) to (row3, b) across the shared DAG, so the
	 * shortest paths of (row3, b) and (row3, b - 1) cannot cross and the
	 * leftmost argmins over the candidate layer c obey the staircase
	 *   opt(a, b - 1) <= opt(a, b) <= opt(a + 1, b)
	 * (the right inequality follows the same non-crossing argument read
	 * from the reverse direction, paths starting at (row3, b)). Every
	 * candidate scan in [opt(a, b - 1), opt(a + 1, b)] therefore
	 * telescopes: rows a are processed in decreasing order, columns in
	 * increasing order, for O(a.rows * a.cols) cell operations instead
	 * of the cubic scan (plus the divide-and-conquer solve of the
	 * bottom row). The monotone-argmin behaviour of every single
	 * output row (opt(a, b - 1) <= opt(a, b)) is what the divide-and-
	 * conquer of min_plus_apply relies on as well.
	 *
	 * @param a Left factor; entries must not exceed INF.
	 * @param b Right factor; entries must not exceed INF.
	 * @return The min-plus product a (x) b.
	 * @pre a.cols == b.rows with a.cols > 0, and the shared position
	 *      ranges agree (a.col_range() == b.row_range()).
	 */
	static DistMatrix min_plus_product(const DistMatrix &a,
	                                   const DistMatrix &b) {
		assert(a.cols == b.rows && a.cols > 0);
		assert(a.col_range() == b.row_range());

		DistMatrix out(a.rows, b.cols, 0, a.row_range(), b.col_range());
		if (b.cols == 0)
			return out;
		argmin_sweep(a, b, out);
		return out;
	}

	/**
	 * @brief Min-plus product of a row vector with this matrix.
	 *
	 * Computes `out(j) = min_i (vec(i) + (*this)(i, j))` with plain integer
	 * addition. Candidates saturate at INF, so entries never exceed INF.
	 *
	 * This matrix must be a layer of an edit-distance DAG driven by a
	 * CostModel-validated cost assignment, so for a fixed vec the
	 * leftmost argmins of vec(i) + (*this)(i, j) are monotone in j (the
	 * same non-crossing shortest-path argument as
	 * min_plus_product): the result is found with divide-and-conquer in
	 * O((rows + cols) log cols) argmin steps instead of the full
	 * rows x cols scan.
	 *
	 * @param vec Row vector of costs, one per row of this matrix; entries
	 *            must not exceed INF.
	 * @return The resulting row vector, one entry per column.
	 * @pre vec.size() == rows.
	 */
	std::vector<int> min_plus_apply(const std::vector<int> &vec) const {
		return min_plus_apply(vec, row_range_value);
	}

	/**
	 * @brief Applies this matrix to a contiguous subrange of its row range.
	 *
	 * The first element of @p vec corresponds to input_range.first. Rows
	 * outside the input range are excluded. This lets callers omit source
	 * states already known to exceed their threshold.
	 *
	 * @param vec Costs for the selected contiguous range of matrix rows.
	 * @param input_range Query-position range represented by @p vec.
	 * @return One entry per matrix column.
	 * @pre input_range lies within the matrix row range and vec.size()
	 *      matches its width.
	 */
	std::vector<int> min_plus_apply(
	    const std::vector<int> &vec,
	    std::pair<std::size_t, std::size_t> input_range) const {
		assert(input_range.first >= row_range_value.first &&
		       input_range.second <= row_range_value.second &&
		       input_range.first <= input_range.second);
		const std::size_t offset = input_range.first - row_range_value.first;
		const std::size_t expected = input_range.second - input_range.first;
		assert(vec.size() == expected);

		std::vector<int> out(cols, INF);
		if (cols == 0)
			return out;
		if (vec.empty())
			return out;
		argmin_row(vec.data(), *this, out.data(), 0, cols - 1,
		           offset, offset + vec.size() - 1, nullptr, offset);
		return out;
	}

private:
	/**
	 * @brief Argmin-window computation of the min-plus product.
	 *
	 * Implements the staircase described in min_plus_product. The bottom
	 * row (a = rows - 1) has no (a + 1) upper bound yet, so it is solved
	 * by the divide-and-conquer once, capturing its leftmost argmins for
	 * the rows above. Rows then run in decreasing a, columns increasing:
	 * each cell scans candidates in [opt(a, j - 1), opt(a + 1, j)].
	 *
	 * Cheap INF emission: per call the sweep computes
	 * first_finite1(a) = first c with in1(a, c) < INF and
	 * first_finite2(c) = first j with in2(c, j) < INF — both with a
	 * forward-only cursor (copy the previous row's value, advance while
	 * the cell is INF), linear in the matrix since the values are
	 * monotone across consecutive ED-layer rows. With the row's first
	 * finite candidate at first_finite1(a), the earliest finite landing
	 * over all remaining candidates is the suffix minimum of
	 * first_finite2 at that index (a suffix-minimum prefix is needed
	 * because real window-clamped layer blocks can have non-monotone
	 * first_finite2 values). For b below that minimum every candidate
	 * saturates at INF, so the cell is emitted in O(1) without
	 * scanning; scans also start at left = first_finite1(a), since no
	 * argmin can be smaller. An all-INF left row (first_finite1 == B)
	 * whole-row INF-emits with a loose bound for the row above.
	 *
	 * @param a Left factor.
	 * @param b Right factor.
	 * @param out Pre-sized product matrix to fill (a.rows x b.cols).
	 */
	static void argmin_sweep(const DistMatrix &a, const DistMatrix &b,
	                         DistMatrix &out) {
		const std::size_t E = a.rows;
		const std::size_t F = b.cols;
		const std::size_t B = a.cols;
		// first_finite1(ra): first candidate c with in1(ra, c) < INF
		// (or B for an all-INF row); first_finite2(c): first column j
		// with in2(c, j) < INF (or F). For ED-DAG layers all-INF rows
		// are reserved pedantically: delete below and the reaches band
		// shifts right monotonically with the index, so both are swept
		// with a forward-only cursor: copy the previous row's value,
		// then advance while the cell is INF. Amortized
		// O(a.rows + B) and O(b.rows + F) over the whole
		// precomputation; an all-INF row records itself as B/F and
		// hands its predecessor's cursor to the next row. Nothing below
		// relies on the monotonicity packing the cursor with (a wrong
		// too-small cursor only costs extra scans).
		std::vector<std::size_t> ff1(a.rows, B);
		{
			std::size_t cur = 0;
			for (std::size_t ra = 0; ra < E; ++ra) {
				std::size_t p = cur;
				while (p < B && a.data[ra * B + p] >= INF)
					++p;
				ff1[ra] = p;
				if (p < B)
					cur = p; // all-INF rows retry at cur
			}
		}
		std::vector<std::size_t> ff2(b.rows, F);
		{
			std::size_t cur = 0;
			for (std::size_t c = 0; c < B; ++c) {
				std::size_t p = cur;
				while (p < F && b.data[c * F + p] >= INF)
					++p;
				ff2[c] = p;
				if (p < F)
					cur = p;
			}
		}
		// suffix_min_ff2[c] = min over c' >= c of first_finite2(c').
		// Columns below this are INF for every candidate from first
		// finite on — the gate needs no monotonicity of ff2 itself
		// (real window-clamped layer blocks can violate it).
		std::vector<std::size_t> smin_ff2(b.rows + 1, F);
		for (std::size_t c = b.rows; c-- > 0;)
			smin_ff2[c] = std::min(ff2[c], smin_ff2[c + 1]);
		// Argmins of the row below the one being processed (loose B - 1
		// until the bottom row is computed for real). Argmins are
		// candidate indices, so any loose bound must be the last
		// candidate row; using an output-row index would read out of
		// range. Processing rows in decreasing a keeps every needed
		// bound one swap behind.
		std::vector<std::size_t> opt_next(F, B - 1);
		std::vector<std::size_t> opt_cur(F);
		for (std::size_t ra = E; ra-- > 0;) {
			const int *arow = &a.data[ra * a.cols];
			int *orow = &out.data[ra * F];
			const std::size_t ff = ff1[ra];
			// First column that can be finite at all for this row:
			// with the row's first finite candidate at ff, every
			// candidate's first finite landing is the suffix
			// minimum from ff, so columns below that are INF
			// outright (smin_ff2[B] == F covers all-INF rows).
			const std::size_t thr = smin_ff2[ff];
			if (ra + 1 == E) {
				// Bottom row: solve with the divide-and-conquer
				// and capture its leftmost argmins, seeded at the
				// row's first finite candidate.
				argmin_row(arow, b, orow, 0, F - 1,
				           ff < B ? ff : 0, a.cols - 1, &opt_cur);
			} else {
				// Columns below the first finite landing are INF
				// without scanning; the loose argmin bound keeps
				// the staircase valid for the row above.
				std::size_t j = 0;
				for (; j < thr; ++j) {
					orow[j] = INF;
					opt_cur[j] = B - 1;
				}
				std::size_t left = ff;
				for (; j < F; ++j) {
					const std::size_t hi = opt_next[j];
					std::size_t best_k = left;
					int best =
					    arow[left] + b.data[left * F + j];
					if (best > INF)
						best = INF;
					for (std::size_t k = left + 1;
					     k <= hi; ++k) {
						int cand =
						    arow[k] +
						    b.data[k * F + j];
						if (cand > INF)
							cand = INF;
						if (cand < best) {
							best = cand;
							best_k = k;
						}
					}
					orow[j] = best;
					// A saturated cell has no
					// meaningful argmin; the
					// loose whole-range bound
					// keeps the staircase valid
					// as the upper bound for the
					// row above.
					opt_cur[j] =
					    best < INF ? best_k : B - 1;
					left = best_k;
				}
			}
			opt_next.swap(opt_cur);
		}
	}

	/**
	 * @brief Divide-and-conquer column minima of M(c, j) = arow(c) + b(c, j).
	 *
	 * With both factors ED-DAG layers, M's leftmost column argmins are
	 * monotone in j (the non-crossing shortest-path staircase), so the
	 * argmin of the middle column bounds the argmin ranges of both
	 * halves. Ties resolve to the leftmost argmin, preserving monotonicity.
	 *
	 * @param arow One row of the left factor, length b.rows.
	 * @param b Right factor.
	 * @param orow Output row to fill, length b.cols.
	 * @param j_lo First column of the range to solve.
	 * @param j_hi Last column of the range to solve (inclusive).
	 * @param k_lo First candidate row of the argmin range.
	 * @param k_hi Last candidate row of the argmin range (inclusive).
	 * @param found Optional buffer of length b.cols — when set, receives
	 *        the leftmost argmin of each solved column (indexed by j).
	 * @param a_offset Candidate-row index represented by arow[0].
	 * @pre j_lo <= j_hi and k_lo <= k_hi.
	 */
	static void argmin_row(const int *arow, const DistMatrix &b, int *orow,
	                       std::size_t j_lo, std::size_t j_hi,
	                       std::size_t k_lo, std::size_t k_hi,
	                       std::vector<std::size_t> *found = nullptr,
	                       std::size_t a_offset = 0) {
		assert(j_lo <= j_hi && k_lo <= k_hi);
		const std::size_t j_mid = (j_lo + j_hi) / 2;
		std::size_t best_k = k_lo;
		int best = arow[k_lo - a_offset] + b.data[k_lo * b.cols + j_mid];
		if (best > INF)
			best = INF;
		for (std::size_t k = k_lo + 1; k <= k_hi; ++k) {
			int cand = arow[k - a_offset] + b.data[k * b.cols + j_mid];
			if (cand > INF)
				cand = INF;
			if (cand < best) {
				best = cand;
				best_k = k;
			}
		}
		orow[j_mid] = best;
		const bool has_argmin = best < INF;
		// A saturated cell has no meaningful argmin. Store the loose upper
		// bound for callers, and retain the full candidate range on both
		// recursive sides instead of narrowing around the arbitrary best_k.
		if (found)
			(*found)[j_mid] = has_argmin ? best_k : k_hi;
		const std::size_t left_k_hi = has_argmin ? best_k : k_hi;
		const std::size_t right_k_lo = has_argmin ? best_k : k_lo;
		if (j_mid > j_lo)
			argmin_row(arow, b, orow, j_lo, j_mid - 1, k_lo, left_k_hi,
			           found, a_offset);
		if (j_mid < j_hi)
			argmin_row(arow, b, orow, j_mid + 1, j_hi, right_k_lo, k_hi,
			           found, a_offset);
	}

	std::size_t rows;      ///< Number of rows.
	std::size_t cols;      ///< Number of columns.
	/// Half-open query-position ranges represented by rows and columns.
	std::pair<std::size_t, std::size_t> row_range_value;
	std::pair<std::size_t, std::size_t> col_range_value;
	std::vector<int> data; ///< Flat row-major storage.
};

} // namespace hlp_grep
