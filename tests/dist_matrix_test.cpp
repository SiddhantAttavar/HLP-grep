/**
 * @file dist_matrix_test.cpp
 * @brief Validates DistMatrix min-plus products against brute force.
 *
 * Checks hand-computed matrix-matrix products, the matrix-vector product,
 * dimension-mismatch and empty-inner-dimension errors, and degenerate
 * empty matrices.
 */

#include <hlp_grep/dist_matrix.hpp>

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>


using namespace hlp_grep;

namespace {

/** Exits with an error message on failure, like the testcase helpers. */
void check(bool ok, const std::string &msg) {
	if (!ok) {
		std::cerr << "FAIL: " << msg << '\n';
		std::exit(1);
	}
}

/** Builds a matrix from row strings; each cell is one digit value. */
DistMatrix make(std::size_t rows, std::size_t cols,
                const std::vector<std::string> &cells) {
	DistMatrix m(rows, cols);
	for (std::size_t i = 0; i < rows; ++i)
		for (std::size_t j = 0; j < cols; ++j)
			m(i, j) = cells[i][j] - '0';
	return m;
}

/** Builds a matrix from row strings of multi-digit values. */
DistMatrix make_values(std::size_t rows, std::size_t cols,
                       const std::vector<std::vector<int>> &cells) {
	DistMatrix m(rows, cols);
	for (std::size_t i = 0; i < rows; ++i)
		for (std::size_t j = 0; j < cols; ++j)
			m(i, j) = cells[i][j];
	return m;
}

/** Compares two matrices entry by entry. */
void check_equal(const DistMatrix &got, const DistMatrix &want,
                 const std::string &name) {
	check(got.num_rows() == want.num_rows() &&
	          got.num_cols() == want.num_cols(),
	      name + ": shape mismatch");
	for (std::size_t i = 0; i < want.num_rows(); ++i)
		for (std::size_t j = 0; j < want.num_cols(); ++j)
			check(got(i, j) == want(i, j),
			      name + ": entry (" + std::to_string(i) + "," +
			          std::to_string(j) + ") is " +
			          std::to_string(got(i, j)) + ", expected " +
			          std::to_string(want(i, j)));
}

void test_product() {
	// Hand-computed 2x2 (x) 2x2.
	const DistMatrix a = make(2, 2, {"02", "39"});
	const DistMatrix b = make(2, 2, {"14", "05"});
	const DistMatrix got = DistMatrix::min_plus_product(a, b);
	// out(0,0) = min(0+1, 2+0) = 1; out(0,1) = min(0+4, 2+5) = 4
	// out(1,0) = min(3+1, 9+0) = 4; out(1,1) = min(3+4, 9+5) = 7
	check_equal(got, make(2, 2, {"14", "47"}), "product");

	// Rectangular shapes: 1x3 (x) 3x2.
	const DistMatrix ra = make(1, 3, {"501"});
	const DistMatrix rb = make(3, 2, {"12", "34", "56"});
	const DistMatrix rgot = DistMatrix::min_plus_product(ra, rb);
	// out(0,0) = min(5+1, 0+3, 1+5) = 3; out(0,1) = min(5+2, 0+4, 1+6) = 4
	check_equal(rgot, make(1, 2, {"34"}), "rectangular product");

	// Non-negative and negative costs alike; large-but-finite values must
	// not be special-cased.
	const DistMatrix la = make_values(2, 1, {{1000000}, {2000000}});
	const DistMatrix lb = make_values(1, 2, {{3, 4}});
	const DistMatrix lgot = DistMatrix::min_plus_product(la, lb);
	check_equal(lgot, make_values(2, 2, {{1000003, 1000004},
	                                     {2000003, 2000004}}),
	            "large values");
}

void test_monge_product() {
	// Both factors Monge: a(i, c) = i + c, b(c, j) = |c - j|.
	// Hand-computed 2x3 (x) 3x3 with b(c, j) = |c - j|:
	// b rows are "012", "101", "210".
	const DistMatrix a = make(2, 3, {"021", "310"});
	const DistMatrix b = make(3, 3, {"012", "101", "210"});
	const DistMatrix want = make(2, 3, {"011", "210"});
	check_equal(DistMatrix::min_plus_product(a, b, true), want,
	            "monge product");
	// Rectangular Monge factors must match the cubic path exactly.
	DistMatrix ra(3, 4);
	for (std::size_t i = 0; i < 3; ++i)
		for (std::size_t c = 0; c < 4; ++c)
			ra(i, c) = static_cast<int>(i + c);
	DistMatrix rb(4, 5);
	for (std::size_t c = 0; c < 4; ++c)
		for (std::size_t j = 0; j < 5; ++j)
			rb(c, j) =
			    static_cast<int>(c > j ? c - j : j - c);
	check_equal(DistMatrix::min_plus_product(ra, rb, true),
	            DistMatrix::min_plus_product(ra, rb, false),
	            "monge product matches cubic");

	// Zero outer dimensions take the early return on the monge path too.
	const DistMatrix ea(0, 2);
	const DistMatrix eb(2, 3);
	const DistMatrix egot = DistMatrix::min_plus_product(ea, eb, true);
	check(egot.num_rows() == 0 && egot.num_cols() == 3,
	      "monge empty: zero-row product shape wrong");
}

void test_monge_inf_fills() {
	// Lifter-style tables: realized span costs with INF below the
	// staircase (the exact FRFF corner that k+1 fills broke).
	// A(5x6): row i has span costs from column i, INF before.
	DistMatrix a(5, 6, DistMatrix::INF);
	for (std::size_t i = 0; i < 5; ++i)
		for (std::size_t k = i; k < 6; ++k)
			a(i, k) = static_cast<int>(k - i);
	// B(6x7): same shape; note B(1,1) = 1 < old k+1-style fills.
	DistMatrix b(6, 7, DistMatrix::INF);
	for (std::size_t k = 0; k < 6; ++k)
		for (std::size_t j = k; j < 7; ++j)
			b(k, j) = static_cast<int>(j - k);
	const DistMatrix cubic = DistMatrix::min_plus_product(a, b, false);
	const DistMatrix fast = DistMatrix::min_plus_product(a, b, true);
	check_equal(fast, cubic, "monge INF fills match cubic");
	// Hand checks: out(0,0) = a(0,0)+b(0,0) = 0; out(0,6) = 0+6 = 6;
	// column 0 of b has only b(0,0) = 0 realized while a(1,0) is INF,
	// so out(1,0) saturates to INF instead of missing the minimum.
	check(fast(0, 0) == 0, "monge INF fills: entry (0,0) wrong");
	check(fast(0, 6) == 6, "monge INF fills: entry (0,6) wrong");
	check(fast(1, 0) == DistMatrix::INF,
	      "monge INF fills: entry (1,0) wrong");
	check(fast(4, 0) == DistMatrix::INF,
	      "monge INF fills: saturation wrong");
}

void test_apply() {
	const DistMatrix m = make(2, 3, {"092", "831"});
	const std::vector<int> vec = {1, 2};
	const std::vector<int> got = m.min_plus_apply(vec);
	// out(0) = min(1+0, 2+8) = 1; out(1) = min(1+9, 2+3) = 5
	// out(2) = min(1+2, 2+1) = 3
	check(got.size() == 3, "apply: result size");
	check(got[0] == 1 && got[1] == 5 && got[2] == 3, "apply: values wrong");
}

void test_monge_apply() {
	// Monge matrix b(c, j) = |c - j|; vec {2, 0, 3}:
	// out(0) = min(2+0, 0+1, 3+2) = 1; out(1) = min(2+1, 0+0, 3+1) = 0
	// out(2) = min(2+2, 0+1, 3+0) = 1
	const DistMatrix m = make(3, 3, {"012", "101", "210"});
	const std::vector<int> got = m.min_plus_apply({2, 0, 3}, true);
	check(got.size() == 3, "monge apply: result size");
	check(got[0] == 1 && got[1] == 0 && got[2] == 1,
	      "monge apply: values wrong");

	// Rectangular Monge matrix must match the full scan exactly.
	DistMatrix big(4, 5);
	for (std::size_t c = 0; c < 4; ++c)
		for (std::size_t j = 0; j < 5; ++j)
			big(c, j) = static_cast<int>(c > j ? c - j : j - c);
	const std::vector<int> fast = big.min_plus_apply({3, 0, 2, 5}, true);
	const std::vector<int> slow = big.min_plus_apply({3, 0, 2, 5}, false);
	check(fast == slow, "monge apply: mismatch against full scan");

	// Lifter-style INF-below-the-staircase table with an INF input
	// entry: saturation and the DnC argmin monotonicity must survive.
	DistMatrix stair(5, 6, DistMatrix::INF);
	for (std::size_t i = 0; i < 5; ++i)
		for (std::size_t j = i; j < 6; ++j)
			stair(i, j) = static_cast<int>(j - i);
	const std::vector<int> sf =
	    stair.min_plus_apply({0, DistMatrix::INF, 1, DistMatrix::INF, 0},
	                         true);
	const std::vector<int> ss =
	    stair.min_plus_apply({0, DistMatrix::INF, 1, DistMatrix::INF, 0},
	                         false);
	check(sf == ss, "monge apply INF: mismatch against full scan");
	check(sf[0] == 0, "monge apply INF: entry 0 wrong");
	check(sf[5] == 1, "monge apply INF: entry 5 wrong");

	// Zero columns take the early return on the monge path too.
	const DistMatrix wide(2, 0);
	check(wide.min_plus_apply({0, 0}, true).empty(),
	      "monge apply: zero-column result not empty");
}

void test_errors() {
	const DistMatrix a(2, 3);
	const DistMatrix b(2, 2);
	bool threw = false;
	try {
		DistMatrix::min_plus_product(a, b);
	} catch (const std::invalid_argument &) {
		threw = true;
	}
	check(threw, "errors: product dimension mismatch not thrown");

	const DistMatrix e(2, 0);
	const DistMatrix f(0, 2);
	threw = false;
	try {
		DistMatrix::min_plus_product(e, f);
	} catch (const std::invalid_argument &) {
		threw = true;
	}
	check(threw, "errors: empty inner dimension not thrown");

	const DistMatrix m(3, 2);
	threw = false;
	try {
		m.min_plus_apply({0, 0});
	} catch (const std::invalid_argument &) {
		threw = true;
	}
	check(threw, "errors: apply dimension mismatch not thrown");

	const DistMatrix z(0, 2);
	threw = false;
	try {
		z.min_plus_apply({});
	} catch (const std::invalid_argument &) {
		threw = true;
	}
	check(threw, "errors: apply empty row vector not thrown");
}

void test_empty() {
	// Zero outer dimensions are fine: the loops simply do not run.
	const DistMatrix a(0, 2);
	const DistMatrix b(2, 3);
	const DistMatrix got = DistMatrix::min_plus_product(a, b);
	check(got.num_rows() == 0 && got.num_cols() == 3,
	      "empty: zero-row product shape wrong");

	const DistMatrix m(2, 0);
	const std::vector<int> got_vec = m.min_plus_apply({0, 0});
	check(got_vec.empty(), "empty: apply result not empty");
}

} // namespace

int main() {
	test_product();
	test_monge_product();
	test_monge_inf_fills();
	test_apply();
	test_monge_apply();
	test_errors();
	test_empty();
	std::cout << "dist_matrix_test: all checks passed\n";
	return 0;
}
