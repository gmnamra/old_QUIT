/*
 *  ThreadPool.cpp
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2013, 2014 Tobias Wood. All rights reserved.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "ThreadPool.h"

using namespace std;

namespace QUIT {

ThreadPool::ThreadPool(const size_t num_threads) :
	m_size{num_threads},
	m_continue{false},
	m_finished{false}
{
	// cout << __PRETTY_FUNCTION__ << " started" << endl;
	m_pool.reserve(num_threads);
	// cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

void ThreadPool::resize(const size_t num_threads) {
	// cout << __PRETTY_FUNCTION__ << " started" << endl;
	m_size = num_threads;
	// cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

void ThreadPool::for_loop(const function<void(size_t)> f,
                          const size_t start, const size_t stop, const size_t step) {
	// cout << __PRETTY_FUNCTION__ << " started" << endl;
	size_t thread_step = m_size * step;
	function<void(size_t)> worker = [&](size_t local) {
		while (m_continue && (local < stop)) {
			f(local);
			local += thread_step;
		}
	};
	m_finished = false;
	m_continue = true;
	registerInterrupt();
	for (size_t thread_start = start; thread_start < thread_step; thread_start += step) {
		m_pool.emplace_back(worker, thread_start);
	}
	for (auto &t: m_pool) {
		t.join();
	}
	deregisterInterrupt();
	// Delete worker threads
	m_pool.resize(0);
	m_finished = true;
	m_continue = false;
	// cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

void ThreadPool::for_loop(const function<void(size_t)> f, const size_t stop) {
	// cout << __PRETTY_FUNCTION__ << " started" << endl;
	for_loop(f, 0, stop, 1);
	// cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

void ThreadPool::for_loop2(const function<void(const size_t, const size_t)> f,
                           const size_t starti, const size_t stopi, const size_t stepi,
                           const size_t startj, const size_t stopj, const size_t stepj) {
	// cout << __PRETTY_FUNCTION__ << " started" << endl;
	size_t thread_stepi = m_size * stepi;

	function<void(const size_t, const size_t)> worker = [&](const size_t wi_start, const size_t wj_start) {
		size_t j = wj_start;
		while (m_continue && (j < stopj)) {
			size_t i = wi_start;
			while (m_continue && (i < stopi)) {
				f(i, j);
				i += thread_stepi;
			}
			j += stepj;
		}
	};

	m_finished = false;
	m_continue = true;
	registerInterrupt();
	for (size_t is = starti; is < thread_stepi; is += thread_stepi) {
		m_pool.emplace_back(worker, is, startj);
	}
	for (auto &t: m_pool) {
		t.join();
	}
	deregisterInterrupt();
	// Delete worker threads
	m_pool.resize(0);
	m_finished = true;
	m_continue = false;
	// cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

void ThreadPool::for_loop2(const function<void(const size_t, const size_t)> f, const size_t stopi, const size_t stopj) {
	// cout << __PRETTY_FUNCTION__ << " started" << endl;
	for_loop2(f, 0, stopi, 1, 0, stopj, 1);
	// cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

void ThreadPool::stop() {
	// cout << __PRETTY_FUNCTION__ << " started" << endl;
	m_continue = false;
	// cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

bool ThreadPool::finished() {
	// cout << __PRETTY_FUNCTION__ << " started" << endl;
	return m_finished;
	// cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

ThreadPool *ThreadPool::InterruptPool;
void ThreadPool::Interrupt(int) {
	// cout << __PRETTY_FUNCTION__ << " started" << endl;
	if (InterruptPool && InterruptPool->m_continue) {
		cout << "Interrupt received, halting threads." << endl;
		cout << "Pressing CTRL+C again will terminate the program." << endl;
		InterruptPool->stop();
	} else if (InterruptPool && !InterruptPool->m_continue) {
		throw(std::runtime_error("Registered ThreadPool is not running, nothing to interrupt."));
	} else if (!InterruptPool) {
		throw(std::runtime_error("No ThreadPool registered to interrupt."));
	}
	// cout << __PRETTY_FUNCTION__ << " finished" << endl;
}
void ThreadPool::registerInterrupt() {
	// cout << __PRETTY_FUNCTION__ << " started" << endl;
	InterruptPool = this;
	signal(SIGINT, ThreadPool::Interrupt);
	// cout << __PRETTY_FUNCTION__ << " finished" << endl;
}
void ThreadPool::deregisterInterrupt() {
	// cout << __PRETTY_FUNCTION__ << " started" << endl;
	InterruptPool = nullptr;
	signal(SIGINT, SIG_DFL);
	// cout << __PRETTY_FUNCTION__ << " finished" << endl;
}


} // End namespace QUIT
