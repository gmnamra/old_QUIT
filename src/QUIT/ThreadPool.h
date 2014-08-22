/*
 *  ThreadPool.h
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2013, 2014 Tobias Wood. All rights reserved.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_THREAD_POOL
#define QUIT_THREAD_POOL

#include <memory>
#include <iostream>
#include <thread>
#include <functional>
#include <vector>

namespace QUIT {

class ThreadPool {
	private:
		static void Interrupt(int);
		static ThreadPool *InterruptPool;

		std::vector<std::thread> m_pool;
		size_t m_size;
		bool m_continue, m_finished;

		void registerInterrupt();
		void deregisterInterrupt();

	public:
		ThreadPool(const size_t num_threads = std::thread::hardware_concurrency());
		ThreadPool(ThreadPool &) = delete;
		ThreadPool(ThreadPool &&) = delete;
		
		void resize(const size_t num_threads);
		void for_loop(const std::function<void(size_t)> f, const size_t start, const size_t stop, const size_t step);
		void for_loop(const std::function<void(size_t)> f, const size_t stop);
		void for_loop2(const std::function<void(const size_t, const size_t)> f,
		               const size_t starti, const size_t stopi, const size_t stepi,
		                const size_t startj, const size_t stopj, const size_t stepj);
		void for_loop2(const std::function<void(const size_t, const size_t)> f, const size_t stopi, const size_t stopj);
		bool finished();
		void stop();
};

} // End namespace QUIT

#endif // QUIT_THREAD_POOL
