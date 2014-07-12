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

#include "Nifti/Nifti.h"
#include "QUIT/QUIT.h"
#include "DESPOT.h"

using namespace std;
using namespace QUIT;

const string usage{
"Usage is: afi [options] input \n\
\
Options:\n\
	--mask, -m file  : Mask input with specified file.\n\
	--out, -o path   : Add a prefix to the output filenames.\n"
};
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	cout << version << endl << credit_me << endl;
	static struct option long_options[] = {
		{"mask", required_argument, 0, 'm'},
		{"out", required_argument, 0, 'o'},
		{0, 0, 0, 0}
	};
	
	try { // MacOS can't cope with uncaught exceptions
	
	int indexptr = 0, c;
	string procPath, outPrefix = "";
	double n, nomFlip;
	vector<double> mask;
	Nifti::File maskFile, inFile;
	while ((c = getopt_long(argc, argv, "m:o:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'm':
				cout << "Reading mask." << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				mask.resize(maskFile.dims().head(3).prod());
				maskFile.readVolumes(mask.begin(), mask.end(), 0, 1);
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case '?': // getopt will print an error message
			default:
				cout << usage << endl;
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != 1) {
		cout << usage << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Opening input file " << argv[optind] << endl;
	inFile.open(argv[optind], Nifti::Mode::Read);
	if (maskFile.isOpen() && !maskFile.header().matchesSpace(inFile.header())) {
		cerr << "Mask dimensions/transform do not match SPGR file." << endl;
		exit(EXIT_FAILURE);
	}
	Agilent::ProcPar pp;
	if (ReadPP(inFile, pp)) {
		if (pp.contains("afi_dummy")) {
			// From Sam Hurley. The sequence is implemented by waiting afi_dummy
			// periods after the first afi_tr.
			n = pp.realValue("afi_dummy") + 1;
		} else if (pp.contains("afiflag")) {
			// It's come from the ge3d_toby sequence
			n = 5;
		} else {
			cerr << "Could not find afi_dummy or afiflag parameters, does not appear to be an AFI sequence. Enter TR2/TR1: ";
			cin >> n;
		}
		nomFlip = pp.realValue("flip1");
		cout << "Read TR2/TR1 ratio of " << n << " and flip-angle " << nomFlip << " degrees from procpar." << endl;
	} else {
		cout << "Enter TR2/TR1 (ratio) and flip-angle (degrees): ";
		cin >> n >> nomFlip;
	}
	nomFlip = nomFlip * M_PI / 180.;
	size_t nVoxels = inFile.dims().head(3).prod();
	vector<double> tr1(nVoxels), tr2(nVoxels);
	inFile.readVolumes(tr1.begin(), tr1.end(), 0, 1);
	inFile.readVolumes(tr2.begin(), tr2.end(), 1, 1);
	inFile.close();
	vector<double> flip(nVoxels);
	vector<double> B1(nVoxels);
	cout << "Allocated output memory." << endl;
	cout << "Processing..." << endl;
	for (size_t vox = 0; vox < nVoxels; vox++) {
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
	Nifti::Header outHdr = inFile.header();
	outHdr.description = version;
	outHdr.setDim(4, 1);
	string outPath = outPrefix + "angle" + OutExt();
	cout << "Writing actual flip angle to " << outPath << "..." << endl;
	Nifti::File outAngle(outHdr, outPath);
	outAngle.writeVolumes(flip.begin(), flip.end(), 0, 1);
	outAngle.close();
	
	outPath = outPrefix + "B1" + OutExt();
	cout << "Writing B1 ratio to " << outPath << "..." << endl;
	Nifti::File outB1(outHdr, outPath);
	outB1.writeVolumes(B1.begin(), B1.end(), 0, 1);
	outB1.close();
	
	cout << "Finished." << endl;
	
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
    return EXIT_SUCCESS;
}

