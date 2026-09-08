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

namespace fs = std::filesystem;

/**
 * @brief A parsed testcase: dictionary, cost model, queries and their thresholds.
 *
 * Follows the format described in tests/testcases/README.md.
 */
struct Testcase {
	std::string alphabet;           ///< Characters of the alphabet, in matrix order.
	CostModel cost;                 ///< Cost model parsed from the file.
	std::vector<std::string> dict;  ///< Dictionary sequences.
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

	// 2. Cost model: insertion/deletion costs, then the |Sigma| x |Sigma|
	//    substitution matrix.
	int ins = 0, del = 0;
	{
		std::istringstream iss(next("insertion/deletion costs"));
		if (!(iss >> ins >> del))
			fail(file, line_no, "expected two integers: <ins> <del>");
	}
	std::vector<std::vector<int>> matrix(sigma, std::vector<int>(sigma));
	for (std::size_t i = 0; i < sigma; ++i) {
		std::istringstream iss(next("cost matrix row"));
		for (std::size_t j = 0; j < sigma; ++j) {
			if (!(iss >> matrix[i][j]))
				fail(file, line_no, "cost matrix row too short");
		}
	}
	tc.cost = CostModel(ins, del, tc.alphabet, std::move(matrix));

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
