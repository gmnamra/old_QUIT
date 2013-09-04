/*
 *  afi_main.cpp
 *
 *  Created by Tobias Wood on 03/08/2012.
 *  Copyright (c) 2012-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <string>
#include <time.h>
#include <getopt.h>

using namespace std;

#include "Nifti.h"
#ifdef AGILENT
	#include "procpar.h"
#endif
#include "DESPOT.h"

const string usage("Usage is: afi [options] input output \n\
\
Options:\n\
	--mask, -m file  : Mask input with specified file.\n");
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	static struct option long_options[] = {
		{"mask", required_argument, 0, 'm'},
		{0, 0, 0, 0}
	};
	
	int indexptr = 0, c;
	string procPath, outPrefix;
	double n, nomFlip;
	vector<double> tr1, tr2, flip, B1, mask;
	Nifti::File maskFile, inFile;
	while ((c = getopt_long(argc, argv, "m:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'm':
				cout << "Reading mask." << endl;
				maskFile.open(optarg, Nifti::Modes::Read);
				mask = maskFile.readVolume<double>(0);
				break;
		}
	}
	if ((argc - optind) != 2) {
		cout << usage << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Opening input file " << argv[optind] << endl;
	inFile.open(argv[optind], Nifti::Modes::Read);
	if (maskFile.isOpen() && !maskFile.matchesSpace(inFile)) {
		cerr << "Mask dimensions/transform do not match SPGR file." << endl;
		exit(EXIT_FAILURE);
	}
	#ifdef AGILENT
	Agilent::ProcPar pp;
	if (ReadPP(inFile, pp)) {
		// From Sam Hurley. The sequence is implemented by waiting afi_dummy
		// periods after the first afi_tr.
		n = pp.realValue("afi_dummy") + 1;
		nomFlip = pp.realValue("flip1");
		cout << "Read TR2/TR1 ratio of " << n << " and flip-angle " << nomFlip << " degrees from procpar." << endl;
	} else
	#endif
	{
		cout << "Enter TR2/TR1 (ratio) and flip-angle (degrees): ";
		cin >> n >> nomFlip;
	}
	nomFlip = nomFlip * M_PI / 180.;
	tr1 = inFile.readVolume<double>(0);
	tr2 = inFile.readVolume<double>(1);
	inFile.close();
	outPrefix = string(argv[++optind]);
	flip.resize(inFile.voxelsPerVolume());
	B1.resize(inFile.voxelsPerVolume());
	cout << "Allocated output memory." << endl;
	cout << "Processing..." << endl;
	for (int vox = 0; vox < inFile.voxelsPerVolume(); vox++) {
		if (!maskFile.isOpen() || mask[vox] > 0.) {
			double r = tr2[vox] / tr1[vox];
			double temp = (r*n - 1.) / (n - r);
			if (temp > 1.)
				temp = 1.;
			if (temp < -1.)
				temp = -1.;
			double alpha = acos(temp);
			flip[vox] = alpha * 180. / M_PI;
			B1[vox]   = alpha / nomFlip;
		}
		else
			B1[vox] = 1.; // So smoothing doesn't get messed up
	}
	Nifti::File outFile(inFile);
	outFile.setDim(4, 1);
	outFile.setDatatype(DT_FLOAT32);
	string outPath = outPrefix + "_flip.nii.gz";
	cout << "Writing actual flip angle to " << outPath << "..." << endl;
	outFile.open(outPath, Nifti::Modes::Write);
	outFile.writeVolume(0, flip);
	outFile.close();
	
	outPath = outPrefix + "_B1.nii.gz";
	cout << "Writing B1 ratio to " << outPath << "..." << endl;
	outFile.open(outPath, Nifti::Modes::Write);
	outFile.writeVolume(0, B1);
	outFile.close();
	
	cout << "Finished." << endl;
    return EXIT_SUCCESS;
}

