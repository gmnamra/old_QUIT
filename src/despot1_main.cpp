/*
 *  despot1_main.cpp
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright (c) 2011-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <time.h>
#include <getopt.h>
#include <iostream>
#include <atomic>
#include <Eigen/Dense>

#include "Nifti/Nifti.h"
#include "Nifti/Volume.h"
#include "DESPOT.h"
#include "ThreadPool.h"

#ifdef AGILENT
#include "procpar.h"
using namespace Agilent;
#endif

using namespace std;
using namespace Eigen;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: despot1 [options] spgr_input \n\
\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print writeResiduals information\n\
	--mask, -m file   : Mask input with specified file\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--B1, -b file     : B1 Map file (ratio)\n\
	--algo, -a l      : LLS algorithm (default)\n\
	           w      : WLLS algorithm\n\
			   n      : NLLS algorithm\n"
};

enum class Algos { LLS, WLLS, NLLS };
static int verbose = false;
static string outPrefix;
static Algos algo;
static struct option long_options[] =
{
	{"B1", required_argument, 0, 'b'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"verbose", no_argument, 0, 'v'},
	{0, 0, 0, 0}
};
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	cout << version << endl << credit_shared << endl;
	double spgrTR = 0.;
	Nifti spgrFile, B1File, maskFile;
	Volume<float> spgrVol, B1Vol;
	Volume<bool> maskVol;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvm:o:b:a:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskVol.readFrom(maskFile);
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'b':
				cout << "Reading B1 file: " << optarg << endl;
				B1File.open(optarg, Nifti::Mode::Read);
				B1Vol.readFrom(B1File);
				break;
			case 'a':
				switch (*optarg) {
					case 'l': algo = Algos::LLS; break;
					case 'w': algo = Algos::WLLS; break;
					case 'n': algo = Algos::NLLS; break;
					default:
						cout << "Unknown algorithm type " << optarg << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 'h':
			case '?': // getopt will print an error message
				exit(EXIT_FAILURE);
		}
	}
	if ((argc - optind) != 1) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		exit(EXIT_FAILURE);
	}
	
	//**************************************************************************
	#pragma mark Gather SPGR data
	//**************************************************************************
	cout << "Opening SPGR file: " << argv[optind] << endl;
	spgrFile.open(argv[optind], Nifti::Mode::Read);
	if ((maskFile.isOpen() && !maskFile.matchesSpace(spgrFile)) ||
	    (B1File.isOpen() && !B1File.matchesSpace(spgrFile))) {
		cerr << "Mask or B1 dimensions/transform do not match SPGR file." << endl;
		exit(EXIT_FAILURE);
	}
	size_t nSPGR = spgrFile.dim(4);
	ArrayXd spgrAngles(nSPGR);
	
	#ifdef AGILENT
	ProcPar pp;
	if (ReadPP(spgrFile, pp)) {
		spgrTR = pp.realValue("tr");
		for (size_t i = 0; i < nSPGR; i++) spgrAngles[i] = pp.realValue("flip1", i);
	} else
	#endif
	{
		cout << "Enter SPGR TR (seconds):"; cin >> spgrTR;
		cout << "Enter SPGR " << nSPGR << " Flip Angles (degrees):";
		for (size_t i = 0; i < nSPGR; i++) cin >> spgrAngles[i];
	}
	spgrAngles *= M_PI / 180.;
	if (verbose)
	{
		cout << "SPGR TR=" << spgrTR
		          << " s. Flip-angles: " << spgrAngles.transpose() * 180. / M_PI << endl;
		cout << "Ouput prefix will be: " << outPrefix << endl;
	}
	//**************************************************************************
	// Allocate memory for slices
	//**************************************************************************	
	size_t voxelsPerSlice = spgrFile.dims().head(2).prod();
	cout << "Reading SPGR data..." << flush;
	spgrVol.readFrom(spgrFile);
	cout << "done." << endl;
	//**************************************************************************
	// Create results data storage
	//**************************************************************************
	Volume<float> T1Vol(spgrVol.dims().head(3)), PDVol(spgrVol.dims().head(3)),
	              SoSVol(spgrVol.dims().head(3));
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	ThreadPool pool;
	for (size_t slice = 0; slice < spgrFile.dim(3); slice++) {
		clock_t loopStart;
		// Read in data
		if (verbose)
			cout << "Starting slice " << slice << "..." << flush;
		loopStart = clock();
		atomic<int> voxCount{0};
		size_t sliceOffset = slice * voxelsPerSlice;
		
		function<void (const int&)> processVox = [&] (const int &vox) {
			if (!maskFile.isOpen() || (maskVol.at(sliceOffset + vox))) {
				voxCount++;
				ArrayXd localAngles(spgrAngles);
				if (B1File.isOpen()) { // Correct for B1
					localAngles *= B1Vol.at(sliceOffset + vox);
				}
				ArrayXd signal = spgrVol.series(sliceOffset + vox).cast<double>();
				VectorXd Y = signal / localAngles.sin();
				MatrixXd X(Y.rows(), 2);
				X.col(0) = signal / localAngles.tan();
				X.col(1).setOnes();
				VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
				T1Vol.at(sliceOffset + vox) = static_cast<float>(-spgrTR / log(b[0]));
				PDVol.at(sliceOffset + vox) = static_cast<float>(b[1] / (1. - b[0]));
				SoSVol.at(sliceOffset + vox) = static_cast<float>((Y - X*b).array().square().sum());
			}
		};
		pool.for_loop(processVox, voxelsPerSlice);
		
		if (verbose) {
			clock_t loopEnd = clock();
			if (voxCount > 0)
				cout << voxCount << " unmasked voxels, CPU time per voxel was "
				          << ((loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC)) << " s, ";
			cout << "finished." << endl;
		}
	}

	if (verbose)
		cout << "Writing results." << endl;
	Nifti outFile(spgrFile, 1);
	outFile.description = version;
	outFile.open(outPrefix + "D1_T1.nii.gz", Nifti::Mode::Write);
	T1Vol.writeTo(outFile);
	outFile.close();
	outFile.open(outPrefix + "D1_PD.nii.gz", Nifti::Mode::Write);
	PDVol.writeTo(outFile);
	outFile.close();
	outFile.open(outPrefix + "D1_SoS.nii.gz", Nifti::Mode::Write);
	SoSVol.writeTo(outFile);
	outFile.close();

	cout << "All done." << endl;
	exit(EXIT_SUCCESS);
}
