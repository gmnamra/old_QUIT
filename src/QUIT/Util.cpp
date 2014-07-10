/*
 *  Util.cpp
 *  Part of the QUantitative Image Toolbox
 *
 *  Copyright (c) 2014 Tobias Wood. All rights reserved.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "QUIT/Util.h"

using namespace std;

namespace QUIT {

bool ReadPP(const Nifti::File &nii, Agilent::ProcPar &pp) {
	const list<Nifti::Extension> &exts = nii.extensions();
	for (auto &e : exts) {
		if (e.code() == NIFTI_ECODE_COMMENT) {
			string s(e.data().begin(), e.data().end());
			stringstream ss(s);
			ss >> pp;
			return true;
		}
	}
	// If we got to here there are no procpar extensions, try the old method
	string path = nii.basePath() + ".procpar";
	ifstream pp_file(path);
	if (pp_file) {
		pp_file >> pp;
		return true;
	}
	return false;
}

time_t printStartTime() {
	time_t theTime = time(NULL);
	char timeStr[256];
	strftime(timeStr, 256, "%H:%M:%S", localtime(&theTime));
	cout << "Started at " << timeStr << endl;
	return theTime;
}

time_t printElapsedTime(const time_t &startTime) {
    time_t theTime = time(NULL);
    char timeStr[256];
    strftime(timeStr, 512, "%H:%M:%S", localtime(&theTime));
    double elapsed = difftime(theTime, startTime);
	cout << "Finished at " << timeStr << ". Elapsed time was " << elapsed << " s." << endl;
	return theTime;
}

void printElapsedClock(const clock_t &startClock, const int voxCount) {
	clock_t endClock = clock();
	float totalMilliseconds = (endClock - startClock) * (1.e3 / CLOCKS_PER_SEC);
	cout << "Total run time: " << totalMilliseconds << " ms" << endl;
	cout << "Average voxel time: " << totalMilliseconds / voxCount << " ms" << endl;
}

void printLoopTime(const clock_t &loopStart, const int voxCount) {
	clock_t loopEnd = clock();
	if (voxCount > 0) {
		cout << voxCount << " unmasked voxels, CPU time per voxel was "
		     << ((loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC)) << " s" << endl;
	}
}

} // End namespace QUIT
