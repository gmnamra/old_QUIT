/*
 *  Util.h
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2014 Tobias Wood. All rights reserved.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef QUIT_UTIL_H
#define QUIT_UTIL_H

#include <string>
#include <map>
#include <vector>
#include <random>
#include <functional>
#include <time.h>

#include <Eigen/Dense>

#include "Agilent/procpar.h"
#include "Nifti/Nifti.h"
#include "Nifti/ExtensionCodes.h"

namespace QUIT {

const std::string &OutExt(); //!< Return the extension stored in $QUIT_EXT
bool ReadPP(const Nifti::File &nii, Agilent::ProcPar &pp);
time_t printStartTime();
time_t printElapsedTime(const time_t &start);
void printElapsedClock(const clock_t &clockStart, const int voxCount);
void printLoopTime(const clock_t &loopStart, const int voxCount);
void checkHeaders(const Nifti::Header &n1, std::vector<Nifti::File> n_other); //!< Throws an exception if the passed in Nifti files do not share same matrix size and transform
Eigen::Vector3f parse_vector(char *str);

template<typename T>
T randNorm(double sigma)
{
  static std::mt19937_64 twister(time(NULL));
  static std::normal_distribution<T> nd(0., sigma);
  return nd(twister);
}


template<typename T> class Read;

template<typename T> class Read {
	public:
	static void FromLine(std::istream &in, T &val) {
		std::string line;
		if (!std::getline(in, line)) {
			throw(std::runtime_error("Failed to read input."));
		}
		FromString(line, val);
	}
	static void FromString(const std::string &s, T &val) {
		std::istringstream stream(s);
		if (!(stream >> val)) {
			throw(std::runtime_error("Failed to parse input: " + s));
		}
	}
};

template<typename T, long N> class Read<Eigen::Array<T, N, 1>> {
	public:
	static void FromLine(std::istream & in, Eigen::Ref<Eigen::Array<T, N, 1>> vals) {
		std::string line;
		if (!std::getline(in, line)) {
			throw(std::runtime_error("Failed to read input."));
		}
		FromString(line, vals);
	}

	static void FromString(const std::string &s, Eigen::Ref<Eigen::Array<T, N, 1>> vals) {
		std::istringstream stream(s);
		for (typename Eigen::Array<T, Eigen::Dynamic, 1>::Index i = 0; i < vals.size(); i++) {
			if (!(stream >> vals[i])) {
				throw(std::runtime_error("Failed to parse input: " + s));
			}
		}
	}
};

} // End namespace QUIT

#endif // QUIT_UTIL_H
