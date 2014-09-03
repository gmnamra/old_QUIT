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

bool ThreadPool::EnableDebug = false;

ThreadPool::ThreadPool(const size_t num_threads) :
	m_size{num_threads},
	m_interrupted{false},
	m_finished{false}
{
	if (EnableDebug) {
		cout << __PRETTY_FUNCTION__ << " started" << endl;
		cout << "m_size: " << m_size << endl;
	}
	m_pool.reserve(m_size);
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

void ThreadPool::resize(const size_t num_threads) {
	if (EnableDebug) {
		cout << __PRETTY_FUNCTION__ << " started" << endl;
		cout << "setting number of threads to: " << num_threads << " was: " << m_size << endl;
	}
	m_size = num_threads;
	m_pool.reserve(m_size);
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

void ThreadPool::for_loop(const function<void(size_t)> f,
                          const size_t start, const size_t stop, const size_t step) {
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " started" << endl;
	size_t thread_step = m_size * step;
	function<void(size_t)> worker = [&](size_t local) {
		while (!m_interrupted && (local < stop)) {
			f(local);
			local += thread_step;
		}
	};
	m_finished = false;
	m_interrupted = false;
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
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

void ThreadPool::for_loop(const function<void(size_t)> f, const size_t stop) {
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " started" << endl;
	for_loop(f, 0, stop, 1);
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

void ThreadPool::for_loop2(const function<void(const size_t, const size_t)> f,
                           const size_t starti, const size_t stopi, const size_t stepi,
                           const size_t startj, const size_t stopj, const size_t stepj) {
	size_t stepj_per_thread = m_size * stepj;
	if (EnableDebug) {
		cout << __PRETTY_FUNCTION__ << " started" << endl;
		cout << "Function: " << &f << endl;
		cout << "starti " << starti << " stopi " << stopi << " stepi " << stepi << endl;
		cout << "starti " << startj << " stopi " << stopj << " stepi " << stepj << endl;
	}

	function<void(const size_t, const size_t)> worker = [&](const size_t wi_start, const size_t wj_start) {
		size_t j = wj_start;
		while (!m_interrupted && (j < stopj)) {
			size_t i = wi_start;
			while (!m_interrupted && (i < stopi)) {
				f(i, j);
				i += stepi;
			}
			j += stepj_per_thread;
		}
	};

	m_finished = false;
	m_interrupted = false;
	registerInterrupt();
	for (size_t startj_thread = startj; startj_thread < startj + stepj_per_thread; startj_thread += stepj) {
		m_pool.emplace_back(worker, starti, startj_thread);
		if (EnableDebug) {
			cout << "Started worker thread " << &(m_pool.back()) << " with starti: " << starti << " startj: " << startj_thread << endl;
		}
	}
	for (auto &t: m_pool) {
		t.join();
		if (EnableDebug) {
			cout << "Joined by thread " << &t << endl;
		}
	}
	deregisterInterrupt();
	// Delete worker threads
	m_pool.resize(0);
	m_finished = true;
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

void ThreadPool::for_loop2(const function<void(const size_t, const size_t)> f, const size_t stopi, const size_t stopj) {
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " started" << endl;
	if (EnableDebug) cout << "Function: " << &f << endl << "stopi " << stopi << " stopj " << stopj << endl;
	for_loop2(f, 0, stopi, 1, 0, stopj, 1);
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

void ThreadPool::stop() {
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " started" << endl;
	m_interrupted = true;
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

bool ThreadPool::finished() {
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " started" << endl;
	return m_finished;
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

bool ThreadPool::interrupted() {
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " started" << endl;
	return m_interrupted;
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " finished" << endl;
}

ThreadPool *ThreadPool::InterruptPool;
void ThreadPool::Interrupt(int) {
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " started" << endl;
	if (InterruptPool && !InterruptPool->m_interrupted) {
		cout << "Interrupt received, halting threads." << endl;
		cout << "Pressing CTRL+C again will terminate the program." << endl;
		InterruptPool->stop();
	} else if (InterruptPool && InterruptPool->m_interrupted) {
		throw(std::runtime_error("Registered ThreadPool is not running, nothing to interrupt."));
	} else if (!InterruptPool) {
		throw(std::runtime_error("No ThreadPool registered to interrupt."));
	}
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " finished" << endl;
}
void ThreadPool::registerInterrupt() {
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " started" << endl;
	InterruptPool = this;
	signal(SIGINT, ThreadPool::Interrupt);
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " finished" << endl;
}
void ThreadPool::deregisterInterrupt() {
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " started" << endl;
	InterruptPool = nullptr;
	signal(SIGINT, SIG_DFL);
	if (EnableDebug) cout << __PRETTY_FUNCTION__ << " finished" << endl;
}


} // End namespace QUIT
