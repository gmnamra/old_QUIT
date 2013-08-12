/*
 *  despot2-fm_main.cpp
 *
 *  Created by Tobias Wood on 2013/08/12.
 *  Copyright (c) 2013 Tobias Wood.
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

#include "Nifti.h"
#include "DESPOT.h"
#include "DESPOT_Functors.h"
#include "RegionContraction.h"
#include "ThreadPool.h"

#ifdef HAVE_NRECON
	#include "procpar.h"
	using namespace Recon;
#endif

using namespace std;
using namespace Eigen;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string credit {
"despot2 - Written by tobias.wood@kcl.ac.uk, based on work by Sean Deoni. \n\
Acknowledgements greatfully received, grant discussions welcome."
};

const string usage {
"Usage is: despot2-fm [options] T1_map ssfp_files\n\
\
Options:\n\
	--help, -h        : Print this message.\n\
	--mask, -m file   : Mask input with specified file.\n\
	--out, -o path    : Add a prefix to the output filenames.\n\
	--B0 file         : B0 Map file.\n\
	--B1 file         : B1 Map file.\n\
	--verbose, -v     : Print slice processing times.\n\
	--start_slice N   : Start processing from slice N.\n\
	--end_slice   N   : Finish processing at slice N.\n\
	--tesla, -t 3     : Enables DESPOT-FM with boundaries suitable for 3T\n\
	            7     : Boundaries suitable for 7T (default)\n\
	            u     : User specified boundaries from stdin.\n"

};

// tesla == 0 means NO DESPOT-FM
static int tesla = 0, fitB0 = false, verbose = false, start_slice = -1, end_slice = -1;
static string outPrefix;
static struct option long_options[] =
{
	{"B0", required_argument, 0, '0'},
	{"B1", required_argument, 0, '1'},
	{"help", no_argument, 0, 'h'},
	{"mask", required_argument, 0, 'm'},
	{"tesla", required_argument, 0, 't'},
	{"verbose", no_argument, 0, 'v'},
	{"start_slice", required_argument, 0, 'S'},
	{"end_slice", required_argument, 0, 'E'},
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
	cout << credit << endl;
	Eigen::initParallel();
	Nifti::File maskFile, B0File, B1File, inFile, savedHeader;
	vector<double> maskData, B0Data, B1Data, T1Data;
	string procPath;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hm:o:vt:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Modes::Read);
				maskData = maskFile.readVolume<double>(0);
				break;
			case '0':
				cout << "Reading B0 file: " << optarg << endl;
				B0File.open(optarg, Nifti::Modes::Read);
				B0Data = B0File.readVolume<double>(0);
				break;
			case '1':
				cout << "Reading B1 file: " << optarg << endl;
				B1File.open(optarg, Nifti::Modes::Read);
				B1Data = B1File.readVolume<double>(0);
				break;
			case 't':
				switch (*optarg) {
					case '3': tesla = 3; break;
					case '7': tesla = 7; break;
					case 'u': tesla = -1; break;
					default:
						cout << "Unknown boundaries type " << optarg << endl;
						abort();
						break;
				}
				cout << "Using " << tesla << "T boundaries." << endl;
				break;
			case 'v':
				verbose = true;
				break;
			case 'S':
				start_slice = atoi(optarg);
				break;
			case 'E':
				end_slice = atoi(optarg);
				break;
			case 0:
				// Just a flag
				break;
			case '?': // getopt will print an error message
			case 'h':
				cout << usage << endl;				
				return EXIT_FAILURE;
		}
	}
	if ((tesla != 0) && !B0File.isOpen())
		fitB0 = true;
	if ((argc - optind) < 2) {
		cout << "Wrong number of arguments. Need at least a T1 map and 1 SSFP file." << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Reading T1 Map from: " << argv[optind] << endl;
	savedHeader.open(argv[optind++], Nifti::Modes::Read);
	T1Data = savedHeader.readVolume<double>(0);
	savedHeader.close();
	if ((maskFile.isOpen() && !savedHeader.matchesSpace(maskFile)) ||
	    (B0File.isOpen() && !savedHeader.matchesSpace(B0File)) ||
		(B1File.isOpen() && !savedHeader.matchesSpace(B1File))){
		cerr << "Dimensions/transforms do not match in input files." << endl;
		exit(EXIT_FAILURE);
	}
	//**************************************************************************
	// Gather SSFP Data
	//**************************************************************************
	int nFlip, nPhases, nResiduals = 0;
	nPhases = argc - optind;
	vector<DESPOTData> data(nPhases);
	int voxelsPerSlice, voxelsPerVolume;
	vector<vector<double>> ssfpData(nPhases);
	VectorXd inFlip;
	double inTR;
	for (size_t p = 0; p < nPhases; p++) {
		cout << "Reading SSFP header from " << argv[optind] << endl;
		inFile.open(argv[optind], Nifti::Modes::Read);
		if (!inFile.matchesSpace(savedHeader)) {
			cerr << "Input file dimensions and/or transforms do not match." << endl;
			exit(EXIT_FAILURE);
		}
		if (p == 0) { // Read nFlip, TR and flip angles from first file
			nFlip = inFile.dim(4);
			inFlip.resize(nFlip);
			voxelsPerSlice = inFile.voxelsPerSlice();
			voxelsPerVolume = inFile.voxelsPerVolume();
			#ifdef HAVE_NRECON
			ParameterList pars;
			if (ReadProcpar(inFile.basePath() + ".procpar", pars)) {
				inTR = RealValue(pars, "tr");
				for (int i = 0; i < nFlip; i++)
					inFlip[i] = RealValue(pars, "flip1", i);
			} else
			#endif
			{
				cout << "Enter SSFP TR (seconds): " << flush;
				cin >> inTR;
				cout << "Enter " << nFlip << " flip angles (degrees): " << flush;
				for (int i = 0; i < nFlip; i++)
					cin >> inFlip[i];
			}
			inFlip *= M_PI / 180.;
		}
		data[p].resize(nFlip);
		data[p].TR = inTR;
		data[p].setFlip(inFlip);
		#ifdef HAVE_NRECON
		ParameterList pars;
		if (ReadProcpar(inFile.basePath() + ".procpar", pars)) {
			data[p].phase = RealValue(pars, "rfphase") * M_PI / 180.;
		} else
		#endif
		{
			cout << "Enter phase-cycling (degrees): " << flush;
			cin >> data[p].phase; data[p].phase *= M_PI / 180.;
		}
		cout << "Reading SSFP data..." << endl;
		ssfpData[p] = inFile.readAllVolumes<double>();
		inFile.close();
		nResiduals += nFlip;
		optind++;
	}
	
	if (optind != argc) {
		cerr << "Unprocessed arguments supplied.\n" << usage;
		exit(EXIT_FAILURE);
	}
	
	// Set up boundaries for DESPOT-FM if needed
	const long nP = fitB0 ? 3 : 2;
	ArrayXd bounds(nP);
	if (tesla != 0) {
		if (tesla > 0) {
			bounds.block(0, 0, 2, 2) = DESPOT2FM::defaultBounds(tesla);
		} else if (tesla < 0) {
			cout << "Enter parameter pairs (low then high)" << endl;
			for (int i = 0; i < nP; i++) {
				cout << DESPOT2FM::names()[i] << ": " << flush;
				cin >> bounds(i, 0) >> bounds(i, 1);
			}
		}
		// If fitting, give a suitable range and allocate results memory
		if (fitB0) {
			bounds(2, 0) =  0.0 / data[0].TR;
			bounds(2, 1) =  0.5 / data[0].TR;
			B0Data.resize(voxelsPerVolume);
		}
	}
	
	if (verbose) {
		cout << "SSFP Angles (deg): " << data[0].flip().transpose() * 180 / M_PI << endl;
		if (tesla != 0)
			cout << "Low bounds: " << bounds.col(0).transpose() << endl
				 << "Hi bounds:  " << bounds.col(1).transpose() << endl;
	}
	//**************************************************************************
	// Set up results data
	//**************************************************************************
	vector<vector<double>> paramsData(nP);
	for (int p = 0; p < nP; p++)
		paramsData[p].resize(voxelsPerVolume);
	if (tesla == 0)
		nResiduals = 1;
	vector<vector<double>> residuals(nResiduals);
	for (int i = 0; i < nResiduals; i++)
		residuals[i].resize(voxelsPerVolume);
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	if ((start_slice < 0) || (start_slice >= inFile.dim(3)))
		start_slice = 0;
	if ((end_slice < 0) || (end_slice > inFile.dim(3)))
		end_slice = inFile.dim(3);
	ThreadPool pool;
    time_t procStart = time(NULL);
	char theTime[512];
	strftime(theTime, 512, "%H:%M:%S", localtime(&procStart));
	cout << "Started processing at " << theTime << endl;
	for (int slice = start_slice; slice < end_slice; slice++) {
		// Read in data
		if (verbose)
			cout << "Starting slice " << slice << "..." << flush;
		
		atomic<int> voxCount{0};
		const int sliceOffset = slice * voxelsPerSlice;
		clock_t loopStart = clock();
		function<void (const int&)> processVox = [&] (const int &vox) {
			// Set up parameters and constants
			double T1 = 0.;
			ArrayXd params(nP); params.setZero();
			ArrayXd resid(nResiduals); resid.setZero();
			if (!maskFile.isOpen() || ((maskData[sliceOffset + vox] > 0.) && (T1Data[sliceOffset + vox] > 0.)))
			{	// Zero T1 causes zero-pivot error.
				voxCount++;
				T1 = T1Data[sliceOffset + vox];
				// Gather signals.
				vector<DESPOTData> localData = data;
				for (int p = 0; p < nPhases; p++) {
					localData[p].f0_off = B0File.isOpen() ? B0Data[sliceOffset + vox] : 0.;
					localData[p].B1 = B1File.isOpen() ? B1Data[sliceOffset + vox] : 1.;
					VectorXd sig(nFlip);
					for (int i = 0; i < nFlip; i++)
						sig(i) = ssfpData[p][i*voxelsPerVolume + sliceOffset + vox];
					localData[p].setSignal(sig);
				}
				
				// DESPOT2-FM
				ArrayXd weights(nResiduals);
				weights.setConstant(1.0);
				DESPOT2FM tc(localData, T1, false, fitB0);
				RegionContraction<DESPOT2FM> rc(tc, bounds, weights);
				rc.optimise(params);
				resid = rc.residuals();
			}
			for (int p = 0; p < nP; p++) {
				paramsData[p][sliceOffset + vox] = params[p];
			}
			for (int i = 0; i < nResiduals; i++) {
				residuals[i][sliceOffset + vox] = resid[i];
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
    time_t procEnd = time(NULL);
    strftime(theTime, 512, "%H:%M:%S", localtime(&procEnd));
	cout << "Finished processing at " << theTime << ". Run-time was " 
	     << difftime(procEnd, procStart) << " s." << endl;
	
	if (tesla == 0) {
		const vector<string> classic_names { "D2_PD", "D2_T2" };
		for (int p = 0; p < 2; p++) {
			savedHeader.open(outPrefix + classic_names[p] + ".nii.gz", Nifti::Modes::Write);
			savedHeader.writeVolume(0, paramsData[p]);
			savedHeader.close();
		}
		savedHeader.setDim(4, nResiduals);
		savedHeader.open(outPrefix + "D2_Residual.nii.gz", Nifti::Modes::Write);
		for (int i = 0; i < nResiduals; i++)
			savedHeader.writeSubvolume(0, 0, 0, i, -1, -1, -1, i+1, residuals[i]);
		savedHeader.close();
	} else {
		for (int p = 0; p < nP; p++) {
			savedHeader.open(outPrefix + DESPOT2FM::names()[p] + ".nii.gz", Nifti::Modes::Write);
			savedHeader.writeVolume(0, paramsData[p]);
			savedHeader.close();
		}
		savedHeader.setDim(4, nResiduals);
		savedHeader.open(outPrefix + "FM_Residual.nii.gz", Nifti::Modes::Write);
		for (int i = 0; i < nResiduals; i++)
			savedHeader.writeSubvolume(0, 0, 0, i, -1, -1, -1, i+1, residuals[i]);
		savedHeader.close();
	}
	return EXIT_SUCCESS;
}
