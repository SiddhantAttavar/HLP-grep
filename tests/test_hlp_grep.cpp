#include <hlp_grep/hlp_grep.hpp>

#include <brute_dp.hpp>

#include <algorithm>
#include <cassert>

static bool contains_dist_at_most(const std::vector<Result> &results, int id, int dist) {
	return std::any_of(results.begin(), results.end(), [&](const Result &r) {
		return r.id == id && r.dist <= dist;
	});
}

static void test(Solver &solver, const std::string &query) {
	for (int k = 0; k < 5; ++k) {
		const auto results = solver.query(query, k);
		for (const Result &r : results)
			assert(r.dist <= k); // every reported match meets the threshold
	}
}

int main() {
	BruteSolver solver({"ACGT", "ACGA", "TTTT", "ACG"});

	// Exact distances when k covers them.
	{
		test(solver, "ACGT");
		const auto results = solver.query("ACGT", 10);
		assert(results.size() == 4);
		assert(contains_dist_at_most(results, 0, 0));
		assert(contains_dist_at_most(results, 1, 1));
		assert(contains_dist_at_most(results, 2, 3));
		assert(contains_dist_at_most(results, 3, 1));
	}

	// Empty query: distance equals the length of the dictionary entry.
	{
		test(solver, "");
		const auto results = solver.query("", 3);
		assert(results.size() == 1);
		assert(contains_dist_at_most(results, 3, 3));
	}

	return 0;
}
