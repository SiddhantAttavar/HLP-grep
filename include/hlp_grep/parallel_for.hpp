/**
 * @file parallel_for.hpp
 * @brief Minimal header-only parallel loop used to speed up construction.
 *
 * The helpers here are deliberately dependency-free: they use only the C++17
 * standard library (std::thread), so the header-only INTERFACE target needs
 * no OpenMP/TBB toolchain flags. Work is distributed dynamically through an
 * atomic counter, which keeps load balanced when per-iteration cost varies
 * (e.g. edge matrices whose windows differ in width).
 */
#pragma once
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace hlp_grep {

/**
 * @brief Number of worker threads to use when the caller passes 0.
 *
 * Falls back to 1 when the hardware concurrency cannot be determined.
 */
inline std::size_t default_num_threads() {
	const unsigned int hc = std::thread::hardware_concurrency();
	return hc == 0 ? 1 : static_cast<std::size_t>(hc);
}

/**
 * @brief Runs `body(i)` for every i in [begin, end) across up to
 *        @p num_threads workers.
 *
 * Iterations are handed out dynamically, so an uneven cost distribution does
 * not strand workers. @p body must be safe to invoke concurrently on distinct
 * indices; it is the caller's responsibility to avoid data races. A single
 * thread (or an empty range) runs inline without spawning any threads.
 *
 * If a body invocation throws, the first exception is captured and rethrown
 * on the calling thread after every worker has joined.
 *
 * @param begin       First index (inclusive).
 * @param end         Last index (exclusive).
 * @param body        Callable invoked as body(i).
 * @param num_threads Worker count; 0 means default_num_threads().
 */
template <typename Body>
void parallel_for(std::size_t begin, std::size_t end, Body body,
                  std::size_t num_threads = 0) {
	if (num_threads == 0)
		num_threads = default_num_threads();
	if (end <= begin)
		return;
	const std::size_t total = end - begin;
	if (num_threads > total)
		num_threads = total;
	if (num_threads <= 1) {
		for (std::size_t i = begin; i < end; ++i)
			body(i);
		return;
	}

	std::atomic<std::size_t> next{begin};
	std::mutex err_mutex;
	std::exception_ptr error;
	auto worker = [&]() {
		for (;;) {
			const std::size_t i =
			    next.fetch_add(1, std::memory_order_relaxed);
			if (i >= end)
				break;
			try {
				body(i);
			} catch (...) {
				std::lock_guard<std::mutex> lock(err_mutex);
				if (!error)
					error = std::current_exception();
				break;
			}
		}
	};

	std::vector<std::thread> pool;
	pool.reserve(num_threads - 1);
	for (std::size_t t = 1; t < num_threads; ++t)
		pool.emplace_back(worker);
	worker();
	for (std::thread &th : pool)
		th.join();
	if (error)
		std::rethrow_exception(error);
}

} // namespace hlp_grep
