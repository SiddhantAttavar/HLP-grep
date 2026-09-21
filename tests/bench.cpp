/**
 * @file bench.cpp
 * @brief Benchmarks a single testcase with one or more methods: build time,
 *        peak memory and per-query runtime for each method, written as one
 *        JSON file.
 *
 * Usage: bench --method hlp_grep [--method naive ...] <testcase-file>
 *              [--out out.json]
 *
 * Each measured method runs in its own forked child process; the parent
 * recovers the child's rusage (ru_maxrss = peak RSS of the whole run, parse
 * included) via wait4 and assembles the JSON. Per query, the child reports
 * the index, threshold k, runtime, match count and an FNV-1a hash of the
 * result vector, so runs of different methods can be cross-checked for
 * identical results.
 */

#include <hlp_grep/hlp_grep.hpp>

#include "naive_solver.hpp"
#include "testcase.hpp"

#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <iostream>
#include <optional>
#include <string>
#include <utility>
#include <vector>

using namespace hlp_grep;

namespace {

/** FNV-1a 64-bit hash of a results vector: bytes of each (id, dist) pair. */
std::uint64_t results_hash(const std::vector<Result> &results) {
	std::uint64_t h = 0xcbf29ce484222325ULL;
	for (const auto &r : results) {
		const auto feed = [&h](auto value) {
			unsigned char bytes[sizeof(value)];
			std::memcpy(bytes, &value, sizeof(value));
			for (const unsigned char b : bytes) {
				h ^= b;
				h *= 0x100000001b3ULL;
			}
		};
		feed(r.id);
		feed(r.dist);
	}
	return h;
}

double elapsed_ms(const std::chrono::steady_clock::time_point begin,
                  const std::chrono::steady_clock::time_point end) {
	return std::chrono::duration<double, std::milli>(end - begin).count();
}

/** Serialised per-query record sent from the child to the parent. */
struct QueryStat {
	std::size_t index;
	int k;
	double runtime_ms;
	std::uint64_t num_matches;
	std::uint64_t hash;
};

/** Collected measurements of one method run. */
struct MethodRun {
	std::string name;
	double build_ms;
	long peak_bytes;
	std::vector<QueryStat> queries;
};

void write_all(int fd, const std::string &s) {
	std::size_t done = 0;
	while (done < s.size()) {
		const long n = write(fd, s.data() + done, s.size() - done);
		if (n <= 0)
			::_exit(1);
		done += static_cast<std::size_t>(n);
	}
}

void report_query(int fd, std::size_t index, int k,
                  const std::vector<Result> &results,
                  const std::chrono::steady_clock::time_point begin,
                  const std::chrono::steady_clock::time_point end) {
	const double ms = elapsed_ms(begin, end);
	write_all(fd, "QUERY " + std::to_string(index) + " " +
	                  std::to_string(k) + " " + std::to_string(ms) + " " +
	                  std::to_string(results_hash(results)) + " " +
	                  std::to_string(results.size()) + "\n");
}

/**
 * @brief Runs the selected method and streams per-query records to @p fd.
 *
 * First line: `BUILD <build_ms>`; then one `QUERY <index> <k> <ms> <hash>
 * <matches>` line per query, in testcase order.
 */
void run_method(int fd, const std::string &method, const Testcase &tc) {
	if (method == "hlp_grep") {
		const auto build_b = std::chrono::steady_clock::now();
		const Solver solver(tc.dict, *tc.cost);
		const auto build_e = std::chrono::steady_clock::now();
		write_all(fd, "BUILD " + std::to_string(elapsed_ms(build_b, build_e)) +
		                  "\n");

		for (std::size_t i = 0; i < tc.queries.size(); ++i) {
			const auto [k, query] = tc.queries[i];
			const auto b = std::chrono::steady_clock::now();
			const auto results = solver.query(query, k);
			const auto e = std::chrono::steady_clock::now();
			report_query(fd, i, k, results, b, e);
		}
	} else { // naive: no index to build.
		write_all(fd, "BUILD 0\n");

		const NaiveSolver solver(tc.dict, *tc.cost, tc.alphabet);
		for (std::size_t i = 0; i < tc.queries.size(); ++i) {
			const auto [k, query] = tc.queries[i];
			const auto b = std::chrono::steady_clock::now();
			const auto results = solver.query(query, k);
			const auto e = std::chrono::steady_clock::now();
			report_query(fd, i, k, results, b, e);
		}
	}

	// _exit: do not run the parent's inherited atexit/stdio flushes.
	::_exit(0);
}

/**
 * @brief Runs one method in a forked child and collects its records.
 *
 * @return Build time, peak memory in bytes and one QueryStat per query.
 *         Reports failure on stderr and returns nullopt on any error.
 */
std::optional<MethodRun> measure_method(const std::string &method,
                                        const Testcase &tc) {
	int pipefd[2];
	if (pipe(pipefd) != 0) {
		std::cerr << "pipe() failed\n";
		return std::nullopt;
	}
	const pid_t pid = fork();
	if (pid < 0) {
		std::cerr << "fork() failed\n";
		return std::nullopt;
	}
	if (pid == 0) {
		close(pipefd[0]);
		run_method(pipefd[1], method, tc); // never returns
	}
	close(pipefd[1]);

	std::string data;
	char buf[4096];
	long n;
	while ((n = read(pipefd[0], buf, sizeof(buf))) > 0)
		data.append(buf, static_cast<std::size_t>(n));
	close(pipefd[0]);

	int status;
	struct rusage ru;
	if (wait4(pid, &status, 0, &ru) < 0) {
		std::cerr << "wait4() failed\n";
		return std::nullopt;
	}
	if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
		std::cerr << method << " run failed\n";
		return std::nullopt;
	}
	// Linux reports ru_maxrss in KiB.
	const long peak_bytes = ru.ru_maxrss * 1024L;

	// Parse the child's protocol: first "BUILD <ms>", then one QUERY line
	// per query in testcase order.
	double build_ms = 0.0;
	std::vector<QueryStat> qstats;
	{
		std::size_t pos = 0;
		auto next_line = [&data, &pos]() {
			const std::size_t nl = data.find('\n', pos);
			if (nl == std::string::npos)
				return std::string();
			std::string s = data.substr(pos, nl - pos);
			pos = nl + 1;
			return s;
		};
		if (sscanf(data.c_str(), "BUILD %lf", &build_ms) != 1) {
			std::cerr << "missing BUILD record from " << method << " child\n";
			return std::nullopt;
		}
		while (pos < data.size()) {
			QueryStat q{};
			if (sscanf(next_line().c_str(), "QUERY %zu %d %lf %lu %lu",
			           &q.index, &q.k, &q.runtime_ms, &q.hash, &q.num_matches)
			    == 5)
				qstats.push_back(q);
		}
	}
	if (qstats.size() != tc.queries.size()) {
		std::cerr << method << ": received " << qstats.size() << " of "
		          << tc.queries.size() << " query records; child aborted early?\n";
		return std::nullopt;
	}
	std::sort(qstats.begin(), qstats.end(),
	          [](const QueryStat &a, const QueryStat &b) {
		          return a.index < b.index;
	          });
	return MethodRun{method, build_ms, peak_bytes, std::move(qstats)};
}

[[noreturn]] void usage(const char *argv0) {
	std::cerr << "usage: " << argv0
	          << " --method hlp_grep [--method naive ...] <testcase-file>"
	          << " [--out out.json]\n";
	std::exit(1);
}

} // namespace

int main(int argc, char **argv) {
	std::vector<std::string> method_names;
	std::string testcase, out_path;
	for (int i = 1; i < argc; ++i) {
		const std::string arg = argv[i];
		if (arg == "--method") {
			if (i + 1 >= argc)
				usage(argv[0]);
			method_names.push_back(argv[++i]);
		} else if (arg == "--out") {
			if (i + 1 >= argc)
				usage(argv[0]);
			out_path = argv[++i];
		} else if (!arg.empty() && arg[0] != '-' && testcase.empty()) {
			testcase = arg;
		} else {
			std::cerr << "unexpected argument: " << arg << '\n';
			usage(argv[0]);
		}
	}
	if (method_names.empty()) {
		std::cerr << "missing --method (repeatable)\n";
		usage(argv[0]);
	}
	for (const auto &method : method_names) {
		if (method != "hlp_grep" && method != "naive") {
			std::cerr << "unknown method: '" << method << "'\n";
			usage(argv[0]);
		}
	}

	if (testcase.empty()) {
		std::cerr << "missing testcase file\n";
		usage(argv[0]);
	}
	const fs::path file(testcase);
	if (!fs::is_regular_file(file)) {
		std::cerr << "not a testcase file: " << file << '\n';
		return 1;
	}
	if (out_path.empty()) {
		fs::path out = file;
		out.replace_extension(".bench.json");
		out_path = out.string();
	}

	const Testcase tc = parse_testcase(file);

	std::vector<MethodRun> runs;
	for (const auto &method : method_names) {
		auto run = measure_method(method, tc);
		if (!run)
			return 1;
		runs.push_back(std::move(*run));
	}

	// Dictionary and query lengths, reported alongside the timings.
	std::vector<std::size_t> dict_lens, query_lens;
	dict_lens.reserve(tc.dict.size());
	for (const auto &seq : tc.dict)
		dict_lens.push_back(seq.size());
	query_lens.reserve(tc.queries.size());
	for (const auto &[k, query] : tc.queries) {
		(void)k;
		query_lens.push_back(query.size());
	}
	auto sorted = dict_lens;
	std::sort(sorted.begin(), sorted.end());

	std::ostringstream json;
	json << std::fixed << std::setprecision(6);
	json << "{\n"
	     << "  \"testcase\": \"" << file.string() << "\",\n"
	     << "  \"dict_size\": " << tc.dict.size() << ",\n"
	     << "  \"num_queries\": " << tc.queries.size() << ",\n";
	if (sorted.empty()) {
		json << "  \"dict_len_min\": 0, \"dict_len_max\": 0,\n"
		     << "  \"dict_len_median\": 0, \"dict_len_mean\": 0,\n";
	} else {
		const long sum = std::accumulate(sorted.begin(), sorted.end(), 0L);
		json << "  \"dict_len_min\": " << sorted.front() << ",\n"
		     << "  \"dict_len_max\": " << sorted.back() << ",\n"
		     << "  \"dict_len_median\": " << sorted[sorted.size() / 2] << ",\n"
		     << "  \"dict_len_mean\": "
		     << static_cast<double>(sum) / sorted.size() << ",\n";
	}
	json << "  \"dict_len\": [";
	for (std::size_t i = 0; i < dict_lens.size(); ++i)
		json << (i ? ", " : "") << dict_lens[i];
	json << "],\n"
	     << "  \"query_len\": [";
	for (std::size_t i = 0; i < query_lens.size(); ++i)
		json << (i ? ", " : "") << query_lens[i];
	json << "],\n";

	json << "  \"methods\": [\n";
	for (std::size_t m = 0; m < runs.size(); ++m) {
		const auto &r = runs[m];
		json << "    {\n"
		     << "      \"method\": \"" << r.name << "\",\n"
		     << "      \"build_time_ms\": " << r.build_ms << ",\n"
		     << "      \"peak_memory_bytes\": " << r.peak_bytes << ",\n"
		     << "      \"queries\": [\n";
		for (std::size_t i = 0; i < r.queries.size(); ++i) {
			const auto &q = r.queries[i];
			json << "        {\"index\": " << q.index << ", \"k\": " << q.k
			     << ", \"len\": " << query_lens[q.index]
			     << ", \"runtime_ms\": " << q.runtime_ms
			     << ", \"num_matches\": " << q.num_matches
			     << ", \"results_hash\": " << q.hash << "}"
			     << (i + 1 < r.queries.size() ? "," : "") << "\n";
		}
		json << "      ]\n"
		     << "    }" << (m + 1 < runs.size() ? "," : "") << "\n";
	}
	json << "  ]\n"
	     << "}\n";

	std::ofstream out(out_path);
	if (!out) {
		std::cerr << "cannot write " << out_path << '\n';
		return 1;
	}
	out << json.str() << std::flush;
	if (!out) {
		std::cerr << "write failed while writing " << out_path << '\n';
		return 1;
	}
	std::cout << "BENCH " << file << " -> " << out_path
	          << " (dict=" << tc.dict.size() << " queries=" << tc.queries.size()
	          << " methods=";
	for (std::size_t m = 0; m < runs.size(); ++m)
		std::cout << (m ? "," : "") << runs[m].name;
	std::cout << ")\n";
	return 0;
}
