//
//  ThreadPool.cpp
//  DESPOT
//
//  Created by Tobias Wood on 29/07/2013.
//
//

#include "ThreadPool.h"

ThreadPool::ThreadPool(const size_t num_threads) {
	m_size = num_threads;
	m_pool.reserve(num_threads);
}

void ThreadPool::resize(const size_t num_threads) {
	m_size = num_threads;
}

void ThreadPool::for_loop(const function<void(size_t)> f,
                          const size_t start, const size_t stop, const size_t step) {
	size_t thread_step = static_cast<size_t>(m_size) * step;
	
	function<void(size_t)> worker = [&](size_t local) {
		while (m_run && (local < stop)) {
			f(local);
			local += thread_step;
		}
	};
	
	m_run = true;
	for (size_t thread_start = start; thread_start < thread_step; thread_start += step) {
		m_pool.emplace_back(worker, thread_start);
	}
	for (auto &t: m_pool) {
		t.join();
	}
	// Delete worker threads
	m_pool.resize(0);
}

void ThreadPool::for_loop(const function<void(size_t)> f, const size_t stop) {
	for_loop(f, 0, stop, 1);
}

void ThreadPool::stop() {
	m_run = false;
}