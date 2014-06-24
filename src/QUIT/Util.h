#ifndef UTIL_H
#define UTIL_H

#include <time.h>

#include "Agilent/procpar.h"
#include "Nifti/Nifti.h"
#include "Nifti/ExtensionCodes.h"

bool ReadPP(const Nifti &nii, Agilent::ProcPar &pp);
time_t printStartTime();
time_t printElapsedTime(const time_t &start);
void printElapsedClock(const clock_t &clockStart, const int voxCount);
void printLoopTime(const clock_t &loopStart, const int voxCount);
#endif // UTIL_H
