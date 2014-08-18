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

template<typename T>
T randNorm(double sigma)
{
  static std::mt19937_64 twister(time(NULL));
  static std::normal_distribution<T> nd(0., sigma);
  return nd(twister);
}

} // End namespace QUIT

#endif // QUIT_UTIL_H
