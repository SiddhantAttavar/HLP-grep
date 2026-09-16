/**
 * @file testcase.hpp
 * @brief Shared parsing of HLP-grep testcase files and collection of testcase
 *        file paths from CLI arguments.
 *
 * The testcase format is described in tests/testcases/README.md.
 */
#pragma once

#include <hlp_grep/hlp_grep.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace hlp_grep;

namespace fs = std::filesystem;

/**
 * @brief A parsed testcase: dictionary, cost model, queries and their thresholds.
 *
 * Follows the format described in tests/testcases/README.md.
 */
struct Testcase {
	std::string alphabet; ///< Characters of the alphabet, in matrix order.
	/// One of the two concrete models below is populated by the parse, and
	/// cost points at the selected one.
	UnitCostModel unit_cost{1, 1};
	MatrixCostModel matrix_cost;
	const CostModel *cost = nullptr; ///< Selected parsed model.
	std::vector<std::string> dict;   ///< Dictionary sequences.
	std::vector<std::pair<int, std::string>> queries; ///< (k, query) pairs.
};

[[noreturn]] inline void fail(const fs::path &file, int line_no,
                              const std::string &msg) {
	std::cerr << file << ':' << line_no << ": " << msg << '\n';
	std::exit(1);
}

inline Testcase parse_testcase(const fs::path &file) {
	std::ifstream in(file);
	if (!in) {
		std::cerr << "cannot open " << file << '\n';
		std::exit(1);
	}

	std::string line;
	int line_no = 0;
	auto next = [&](const char *what) -> std::string {
		if (!std::getline(in, line)) {
			std::cerr << file << ": unexpected EOF while reading " << what << '\n';
			std::exit(1);
		}
		++line_no;
		return line;
	};

	Testcase tc;

	// 1. Alphabet.
	tc.alphabet = next("alphabet");
	if (tc.alphabet.empty())
		fail(file, line_no, "alphabet must not be empty");
	const std::size_t sigma = tc.alphabet.size();

	// 2. Cost model: per-character insertion costs, per-character deletion
	//    costs, then the |Sigma| x |Sigma| substitution matrix.
	auto read_costs = [&](const char *what) {
		std::vector<int> costs(sigma);
		std::istringstream iss(next(what));
		for (std::size_t i = 0; i < sigma; ++i)
			if (!(iss >> costs[i]))
				fail(file, line_no, "expected " + std::to_string(sigma) +
				                       " cost values");
		return costs;
	};
	std::vector<int> ins_costs = read_costs("insertion costs");
	std::vector<int> del_costs = read_costs("deletion costs");
	std::vector<std::vector<int>> matrix(sigma, std::vector<int>(sigma));
	for (std::size_t i = 0; i < sigma; ++i) {
		std::istringstream iss(next("cost matrix row"));
		for (std::size_t j = 0; j < sigma; ++j) {
			if (!(iss >> matrix[i][j]))
				fail(file, line_no, "cost matrix row too short");
		}
	}
	// Uniform match/mismatch matrices map to the matrix-free model; any
	// other matrix keeps the explicit matrix model.
	bool uniform = true;
	const int diag0 = matrix[0][0];
	const int off0 = matrix[0][1];
	for (std::size_t i = 0; i < sigma && uniform; ++i)
		for (std::size_t j = 0; j < sigma; ++j) {
			const int expect = (i == j) ? diag0 : off0;
			if (matrix[i][j] != expect) {
				uniform = false;
				break;
			}
		}
	if (uniform) {
		tc.unit_cost = UnitCostModel(ins_costs[0], del_costs[0], diag0, off0);
		tc.cost = &tc.unit_cost;
	} else {
		tc.matrix_cost = MatrixCostModel(
		    tc.alphabet, std::move(matrix), std::move(ins_costs),
		    std::move(del_costs));
		tc.cost = &tc.matrix_cost;
	}

	// 3. Dictionary header + sequences.
	{
		const std::string header = next("dictionary size");
		std::size_t n = 0;
		try {
			n = std::stoul(header);
		} catch (...) {
			fail(file, line_no, "dictionary size is not an integer: " + header);
		}
		tc.dict.reserve(n);
		for (std::size_t i = 0; i < n; ++i)
			tc.dict.push_back(next("dictionary sequence"));
	}

	// 4. Query header + queries ("<k> <seq>" per line).
	{
		const std::string header = next("query count");
		std::size_t q = 0;
		try {
			q = std::stoul(header);
		} catch (...) {
			fail(file, line_no, "query count is not an integer: " + header);
		}
		tc.queries.reserve(q);
		for (std::size_t i = 0; i < q; ++i) {
			std::string qline = next("query");
			const auto sp = qline.find(' ');
			if (sp == std::string::npos)
				fail(file, line_no, "expected '<k> <query>': " + qline);
			int k = 0;
			try {
				k = std::stoi(qline.substr(0, sp));
			} catch (...) {
				fail(file, line_no, "k is not an integer: " + qline);
			}
			tc.queries.emplace_back(k, qline.substr(sp + 1));
		}
	}
	return tc;
}

/** Collects testcase files from CLI args: each arg is either a file or a folder. */
inline std::vector<fs::path> collect_files(int argc, char **argv, int first = 1) {
	std::vector<fs::path> files;
	for (int i = first; i < argc; ++i) {
		const fs::path p(argv[i]);
		std::error_code ec;
		if (fs::is_regular_file(p, ec)) {
			files.push_back(p);
		} else if (fs::is_directory(p, ec)) {
			for (const auto &entry : fs::directory_iterator(p, ec))
				if (entry.is_regular_file(ec) && entry.path().extension() == ".txt")
					files.push_back(entry.path());
		} else {
			std::cerr << "not a file or directory: " << p << '\n';
			std::exit(1);
		}
	}
	std::sort(files.begin(), files.end());
	return files;
}

/** Returns the solution file path for a testcase file (`<name>.sol`). */
inline fs::path solution_path(const fs::path &testcase) {
	fs::path sol = testcase;
	sol.replace_extension(".sol");
	return sol;
}
