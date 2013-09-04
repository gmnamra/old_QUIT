/*
 *  despot2_main.cpp
 *
 *  Created by Tobias Wood on 23/01/2012.
 *  Copyright (c) 2012-2013 Tobias Wood.
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

#ifdef AGILENT
	#include "procpar.h"
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
"Usage is: despot2 [options] T1_map ssfp_file\n\
\n\
Options:\n\
	--help, -h        : Print this message.\n\
	--mask, -m file   : Mask input with specified file.\n\
	--out, -o path    : Add a prefix to the output filenames.\n\
	--B0 file         : B0 Map file.\n\
	--B1 file         : B1 Map file.\n\
	--verbose, -v     : Print slice processing times.\n\
	--start_slice N   : Start processing from slice N.\n\
	--end_slice   N   : Finish processing at slice N.\n"
};

static int verbose = false, start_slice = -1, end_slice = -1;
static string outPrefix;
static struct option long_options[] =
{
	{"B0", required_argument, 0, '0'},
	{"B1", required_argument, 0, '1'},
	{"help", no_argument, 0, 'h'},
	{"mask", required_argument, 0, 'm'},
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
	while ((c = getopt_long(argc, argv, "hm:o:v:", long_options, &indexptr)) != -1) {
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
			case 'v': verbose = true; break;
			case 'S': start_slice = atoi(optarg); break;
			case 'E': end_slice = atoi(optarg); break;
			case 0:
				// Just a flag
				break;
			case '?': // getopt will print an error message
			case 'h':
				cout << usage << endl;				
				exit(EXIT_SUCCESS);
		}
	}
	if ((argc - optind) != 2) {
		cout << "Wrong number of arguments. Need a least a T1 map and 1 SSFP (180 degree phase cycling) file." << endl;
		cout << usage << endl;
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
	int nFlip, nResiduals = 0;
	int voxelsPerSlice, voxelsPerVolume;
	VectorXd inFlip;
	double inTR;
	cout << "Reading SSFP header from " << argv[optind] << endl;
	inFile.open(argv[optind], Nifti::Modes::Read);
	nFlip = inFile.dim(4);
	inFlip.resize(nFlip);
	voxelsPerSlice = inFile.voxelsPerSlice();
	voxelsPerVolume = inFile.voxelsPerVolume();
	#ifdef AGILENT
	Agilent::ProcPar pp;
	if (ReadPP(inFile, pp)) {
		inTR = pp.realValue("tr");
		for (int i = 0; i < nFlip; i++)
			inFlip[i] = pp.realValue("flip1", i);
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
	cout << "Reading SSFP data..." << endl;
	vector<double> ssfpData = inFile.readAllVolumes<double>();
	inFile.close();
	nResiduals += nFlip;
	optind++;
	if (verbose) {
		cout << "SSFP Angles (deg): " << inFlip.transpose() * 180 / M_PI << endl;
	}
	//**************************************************************************
	// Set up results data
	//**************************************************************************
	vector<double> PDData(voxelsPerVolume), T2Data(voxelsPerVolume), residualData(voxelsPerVolume);
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
			double PD = 0., T2 = 0., resid = 0.;
			if (!maskFile.isOpen() || ((maskData[sliceOffset + vox] > 0.) && (T1Data[sliceOffset + vox] > 0.)))
			{	// Zero T1 causes zero-pivot error.
				voxCount++;
				double T1 = T1Data[sliceOffset + vox];
				// Gather signals.
				double B1 = B1File.isOpen() ? B1Data[sliceOffset + vox] : 1.;
				VectorXd sig(nFlip);
				for (int i = 0; i < nFlip; i++) {
					sig(i) = ssfpData.at(i*voxelsPerVolume + sliceOffset + vox);
				}
				resid = classicDESPOT2(inFlip, sig, inTR, T1, B1, PD, T2);
			}
			PDData.at(sliceOffset + vox) = PD;
			T2Data.at(sliceOffset + vox) = T2;
			residualData.at(sliceOffset + vox) = resid;
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
	
	const vector<string> classic_names { "D2_PD", "D2_T2" };
	savedHeader.open(outPrefix + "D2_PD.nii.gz", Nifti::Modes::Write);
	savedHeader.writeVolume(0, PDData);
	savedHeader.close();
	savedHeader.open(outPrefix + "D2_T2.nii.gz", Nifti::Modes::Write);
	savedHeader.writeVolume(0, T2Data);
	savedHeader.close();
	savedHeader.open(outPrefix + "D2_Residual.nii.gz", Nifti::Modes::Write);
	savedHeader.writeVolume(0, residualData);
	savedHeader.close();
	cout << "Finished writing data." << endl;
	exit(EXIT_SUCCESS);
}
