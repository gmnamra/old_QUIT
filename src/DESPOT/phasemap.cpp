/*
 *  phasemap_main.cpp
 *
 *  Created by Tobias Wood on 20/06/2012.
 *  Copyright (c) 2012-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <string>
#include <getopt.h>

#include "Nifti/Nifti.h"
#include "QUIT/QUIT.h"
#include "DESPOT.h"

using namespace std;
using namespace QUIT;

const string usage
{
"Usage is: phasemap input_1 input_2 outprefix\n\
\n\
Echo times will be read from procpar if present.\n\
Options:\n\
	--mask, -m mask_file : Mask input with specified file\n\
	--phasetime T        : Calculate the phase accumulated in time T\n\
	--smooth             : Smooth output with a gaussian.\n"
};

int main(int argc, char** argv) {
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	cout << version << endl << credit_me << endl;
	static int smooth = false;
	static struct option long_options[] =
	{
		{"mask", required_argument, 0, 'm'},
		{"phasetime", required_argument, 0, 'p'},
		{"smooth", no_argument, &smooth, true},
		{0, 0, 0, 0}
	};
	
	int indexptr = 0, c;
	double TE1, TE2, deltaTE, phasetime = 0.;
	vector<double> data1, data2, B0, mask;
	Nifti maskFile, inFile;
	while ((c = getopt_long(argc, argv, "m:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'm':
				cout << "Reading mask from " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				mask.resize(maskFile.dims().head(3).prod());
				maskFile.readVolumes(mask.begin(), mask.end(), 0, 1);
				break;
			case 'p':
				phasetime = atof(optarg);
				break;
		}
	}
	if ((argc - optind) == 2) {
		cout << "Opening input file " << argv[optind] << "." << endl;
		inFile.open(argv[optind], Nifti::Mode::Read);
		if (maskFile.isOpen() && !maskFile.matchesSpace(inFile)) {
			cerr << "Mask dimensions/transform do not match input file." << endl;
			exit(EXIT_FAILURE);
		}
		Agilent::ProcPar pp;
		if (ReadPP(inFile, pp)) {
			TE1 = pp.realValue("te", 0);
			TE2 = pp.realValue("te", 1);
		} else {
			cout << "Enter TE2 & TE2 (seconds): ";
			cin >> TE1 >> TE2;
		}
		data1.resize(inFile.dims().head(3).prod());
		data2.resize(inFile.dims().head(3).prod());
		inFile.readVolumes(data1.begin(), data1.end(), 0, 1);
		inFile.readVolumes(data2.begin(), data2.end(), 1, 1);
		inFile.close();
	} else if ((argc - optind) == 3) {
		cout << "Opening input file 1" << argv[optind] << "." << endl;
		inFile.open(argv[optind], Nifti::Mode::Read);
		if (maskFile.isOpen() && !maskFile.matchesSpace(inFile)) {
			cerr << "Mask dimensions/transform do not match input file." << endl;
			exit(EXIT_FAILURE);
		}
		Agilent::ProcPar pp;
		if (ReadPP(inFile, pp)) {
			TE1 = pp.realValue("te", 0);
		} else {
			cout << "Enter TE1 (seconds): ";
			cin >> TE1;
		}
		data1.resize(inFile.dims().head(3).prod());
		inFile.readVolumes(data1.begin(), data1.end(), 0, 1);
		inFile.close();
		cout << "Opening input file 2" << argv[++optind] << "." << endl;
		inFile.open(argv[optind], Nifti::Mode::Read);
		if (maskFile.isOpen() && !maskFile.matchesSpace(inFile)) {
			cerr << "Mask dimensions/transform do not match input file." << endl;
			exit(EXIT_FAILURE);
		}
		if (ReadPP(inFile, pp)){
			TE2 = pp.realValue("te", 0);
		} else {
			cout << "Enter TE2 (seconds): ";
			cin >> TE2;
		}
		data2.resize(inFile.dims().head(3).prod());
		inFile.readVolumes(data2.begin(), data2.end(), 0, 1);
		inFile.close();
	} else {
		cerr << usage << endl;
		exit(EXIT_FAILURE);
	}
	string outPrefix(argv[++optind]);
	if (TE2 < TE1) {	// Swap them
		fprintf(stdout, "TE2 < TE1, swapping.\n");
		data1.swap(data2);
		double tmpTE = TE2;
		TE2 = TE1;
		TE1 = tmpTE;
	}
	deltaTE = TE2 - TE1;
	cout << "Delta TE = " << deltaTE << endl;
	B0.resize(inFile.dims().head(3).prod());
	cout << "Processing..." << endl;
	for (size_t vox = 0; vox < inFile.dims().head(3).prod(); vox++) {
		if (!maskFile.isOpen() || mask[vox] > 0.) {
			double deltaPhase = data2[vox] - data1[vox];
			B0[vox] = deltaPhase / (2 * M_PI * deltaTE);
			if (phasetime > 0.) {
				double ph = fmod(B0[vox] * 2 * M_PI * phasetime, 2 * M_PI);
				if (ph > M_PI) ph -= (2 * M_PI);
				if (ph < -M_PI) ph += (2 * M_PI);
				B0[vox] = ph;
			}
		}
	}
	cout << "Writing off-resonance map (in Hz)." << endl;
	string outPath = outPrefix + "f0.nii.gz";
	
	Nifti outFile(inFile, 1);
	outFile.description = version;
	outFile.open(outPath, Nifti::Mode::Write);
	outFile.writeVolumes(B0.begin(), B0.end(), 0, 1);
	outFile.close();
	cout << "Finished." << endl;
    return EXIT_SUCCESS;
}

