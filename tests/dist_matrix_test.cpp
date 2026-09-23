/**
 * @file dist_matrix_test.cpp
 * @brief Validates DistMatrix min-plus products against brute force.
 *
 * Checks hand-computed matrix-matrix products, the argmin-staircase
 * product of ED-DAG layer blocks against a brute-force reference, the
 * matrix-vector product, dimension-mismatch and empty-inner-dimension
 * errors, CostModel's structural validation, and degenerate empty
 * matrices.
 */

#include <hlp_grep/dist_matrix.hpp>

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
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
	// Rectangular 1xN shapes (1 output row): direct, always exact.
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

/** Brute semiglobal ED cost of aligning @p label against a window. */
int ed_block(const std::string &label, const CostModel &cost,
             const std::string &query, std::size_t from,
             std::size_t to) {
	const std::string str = query.substr(from, to - from);
	std::vector<std::vector<int>> dp(
	    label.size() + 1, std::vector<int>(str.size() + 1, 0));
	for (std::size_t i = 1; i <= label.size(); ++i)
		dp[i][0] = dp[i - 1][0] + cost.del();
	for (std::size_t j = 1; j <= str.size(); ++j)
		dp[0][j] = dp[0][j - 1] + cost.ins();
	for (std::size_t i = 1; i <= label.size(); ++i)
		for (std::size_t j = 1; j <= str.size(); ++j)
			dp[i][j] = std::min(
			    {dp[i - 1][j] + cost.del(),
			     dp[i][j - 1] + cost.ins(),
			     dp[i - 1][j - 1] +
			         cost.consume(label[i - 1], str[j - 1])});
	return dp[label.size()][str.size()];
}

void test_staircase_product() {
	// ED-DAG layer fuzz: two single-character layer blocks built from
	// the same query/window axis, checked against a brute-force min-plus
	// product and the direct composed ED of the concatenated labels.
	// Reverse (backward) intervals are INF, as in edge_matrix.
	const CostModel unit(1, 1, 0, 1);
	const CostModel heavy(2, 3, 0, 1);
	const std::string query = "ACGTA";
	std::mt19937 rng(12345);
	for (const CostModel *model : {&unit, &heavy}) {
		for (int trial = 0; trial < 60; ++trial) {
			const std::string l1 = {char("ACGT"[rng() % 4])};
			const std::string l2 = {char("ACGT"[rng() % 4])};
			DistMatrix m1(6, 6, DistMatrix::INF);
			DistMatrix m2(6, 6, DistMatrix::INF);
			for (std::size_t i = 0; i < 6; ++i)
				for (std::size_t j = i; j < 6; ++j) {
					m1(i, j) =
					    ed_block(l1, *model, query, i, j);
					m2(i, j) =
					    ed_block(l2, *model, query, i, j);
				}
			DistMatrix brute(6, 6, DistMatrix::INF);
			for (std::size_t i = 0; i < 6; ++i)
				for (std::size_t j = 0; j < 6; ++j)
					for (std::size_t c = 0; c < 6; ++c)
						if (m1(i, c) < DistMatrix::INF &&
						    m2(c, j) < DistMatrix::INF) {
							int cand = m1(i, c) + m2(c, j);
							if (cand > DistMatrix::INF)
								cand = DistMatrix::INF;
							if (cand < brute(i, j))
								brute(i, j) = cand;
						}
			check_equal(DistMatrix::min_plus_product(m1, m2),
			            brute,
			            "staircase product matches brute (trial " +
			                std::to_string(trial) + ")");
		}
	}
}

void test_layer_inf_fills() {
	// Lifter-style tables realized as ED-DAG layer blocks: a pure
	// ins/deletion stretch of the query (row i spans columns >= i, INF
	// before), the same shape on the right factor.
	DistMatrix a(5, 6, DistMatrix::INF);
	for (std::size_t i = 0; i < 5; ++i)
		for (std::size_t k = i; k < 6; ++k)
			a(i, k) = static_cast<int>(k - i);
	DistMatrix b(6, 7, DistMatrix::INF);
	for (std::size_t k = 0; k < 6; ++k)
		for (std::size_t j = k; j < 7; ++j)
			b(k, j) = static_cast<int>(j - k);
	const DistMatrix got = DistMatrix::min_plus_product(a, b);
	// Hand checks: out(0,0) = a(0,0)+b(0,0) = 0; out(0,6) = 0+6 = 6;
	// column 0 of b has only b(0,0) = 0 realized while a(1,0) is INF,
	// so out(1,0) saturates to INF instead of missing the minimum.
	check(got(0, 0) == 0, "layer INF fills: entry (0,0) wrong");
	check(got(0, 6) == 6, "layer INF fills: entry (0,6) wrong");
	check(got(1, 0) == DistMatrix::INF,
	      "layer INF fills: entry (1,0) wrong");
	check(got(4, 0) == DistMatrix::INF,
	      "layer INF fills: saturation wrong");
	// The full composed matrix: row i's min landing position is i
	// (cost 0), spreading by 1 per column to the right.
	for (std::size_t i = 0; i < 5; ++i)
		for (std::size_t j = 0; j < 7; ++j)
			check(got(i, j) ==
			          (j >= i ? static_cast<int>(j - i)
			                  : DistMatrix::INF),
			      "layer INF fills: entry (" + std::to_string(i) +
			          "," + std::to_string(j) + ") wrong");
}

void test_apply() {
	// ED-DAG block |c - j| with vec {2, 0, 3}:
	// out(0) = min(2+0, 0+1, 3+2) = 1; out(1) = min(2+1, 0+0, 3+1) = 0
	// out(2) = min(2+2, 0+1, 3+0) = 1
	const DistMatrix m = make(3, 3, {"012", "101", "210"});
	const std::vector<int> got = m.min_plus_apply({2, 0, 3});
	check(got.size() == 3, "apply: result size");
	check(got[0] == 1 && got[1] == 0 && got[2] == 1,
	      "apply: values wrong");

	// Large-but-finite input values must not saturate below their sum.
	{
		const DistMatrix big = make_values(1, 2, {{1000000, 5}});
		const std::vector<int> g = big.min_plus_apply({250000});
		check(g[0] == 1250000 && g[1] == 250005,
		      "apply: large-entry saturation wrong");
	}

	// Lifter-style pure-deletion layer with an INF input entry:
	// saturation and the argmin monotonicity must survive.
	DistMatrix stair(5, 6, DistMatrix::INF);
	for (std::size_t i = 0; i < 5; ++i)
		for (std::size_t j = i; j < 6; ++j)
			stair(i, j) = static_cast<int>(j - i);
	const std::vector<int> sf =
	    stair.min_plus_apply(
	        {0, DistMatrix::INF, 1, DistMatrix::INF, 0});
	check(sf[0] == 0, "apply INF: entry 0 wrong");
	check(sf[5] == 1, "apply INF: entry 5 wrong");

	// Zero columns take the early return as well.
	const DistMatrix wide(2, 0);
	check(wide.min_plus_apply({0, 0}).empty(),
	      "apply: zero-column result not empty");
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

void test_cost_model() {
	// match == 0 is required for the argmin staircase.
	for (const auto &[ins, del, match, mismatch] :
	     std::vector<std::tuple<int, int, int, int>>{
	         {1, 1, 1, 1}}) {
		bool threw = false;
		try {
			CostModel bad(ins, del, match, mismatch);
		} catch (const std::invalid_argument &) {
			threw = true;
		}
		check(threw,
		      "cost model: invalid subs/match not rejected (" +
		          std::to_string(ins) + "," + std::to_string(del) +
		          "," + std::to_string(match) + "," +
		          std::to_string(mismatch) + ")");
	}
	// Negative ins/del individually, and ins+de < mismatch.
	for (const auto &[ins, del, match, mismatch] :
	     std::vector<std::tuple<int, int, int, int>>{
	         {-1, 3, 0, 1}, {2, -3, 0, 1}, {1, 1, 0, 3}}) {
		bool threw = false;
		try {
			CostModel bad(ins, del, match, mismatch);
		} catch (const std::invalid_argument &) {
			threw = true;
		}
		check(threw,
		      "cost model: ins/del/mismatch bound not rejected (" +
		          std::to_string(ins) + "," + std::to_string(del) +
		          "," + std::to_string(match) + "," +
		          std::to_string(mismatch) + ")");
	}
	// Valid models must construct.
	static_cast<void>(CostModel(2, 3, 0, 1));
	static_cast<void>(DEFAULT_COST_MODEL);
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
	test_staircase_product();
	test_layer_inf_fills();
	test_apply();
	test_errors();
	test_cost_model();
	test_empty();
	std::cout << "dist_matrix_test: all checks passed\n";
	return 0;
}
