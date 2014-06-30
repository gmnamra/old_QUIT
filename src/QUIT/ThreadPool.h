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

#include <iostream>
#include <thread>
#include <functional>
#include <vector>

namespace QUIT {

class ThreadPool {
	private:
		std::vector<std::thread> m_pool;
		size_t m_size;
		bool m_run;

	public:
		ThreadPool(const size_t num_threads = std::thread::hardware_concurrency());
		ThreadPool(ThreadPool &) = delete;
		ThreadPool(ThreadPool &&) = delete;
		
		void resize(const size_t num_threads);
		void for_loop(const std::function<void(size_t)> f, const size_t start, const size_t stop, const size_t step);
		void for_loop(const std::function<void(size_t)> f, const size_t stop);
		void stop();
};

} // End namespace QUIT

#endif // QUIT_THREAD_POOL
