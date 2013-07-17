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

#include "Nifti.h"
#include "DESPOT.h"

#ifdef HAVE_NRECON
#include "procpar.h"
using namespace Recon;
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
	-m, --mask file : Mask input with specified file.\n\
	--out, -o path    : Add a prefix to the output filenames.\n\
	--B1 file       : Correct flip angles with specified B1 ratio.\n\
	-v, --verbose   : Print out more messages.\n\
	-d, --drop      : Drop certain flip-angles (Read from stdin).\n"
};

const string credit {
"despot1 - Written by tobias.wood@kcl.ac.uk, based on work by Sean Deoni. \n\
Acknowledgements greatfully received, grant discussions welcome."
};

static int verbose = false, drop = false;
static string outPrefix;
static struct option long_options[] =
{
	{"B1", required_argument, 0, '1'},
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
	cout << credit << endl;
	double spgrTR = 0.;
	int nSPGR;
	vector<double> B1Data, maskData;
	Nifti::File spgrFile, B1File, maskFile;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "m:o:vd", long_options, &indexptr)) != -1) {
		switch (c) {
			case '1':
				cout << "Opening B1 file: " << optarg << endl;
				B1File.open(optarg, 'r');
				B1Data = B1File.readVolume<double>(0);
				B1File.close();
				break;
			case 'm':
				cout << "Opening mask file: " << optarg << endl;
				maskFile.open(optarg, 'r');
				maskData = maskFile.readVolume<double>(0);
				maskFile.close();
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'v': verbose = true; break;
			case 'd': drop = true; break;
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
	spgrFile.open(argv[optind], 'r');
	if ((maskFile.isOpen() && !maskFile.matchesSpace(spgrFile)) ||
	    (B1File.isOpen() && !B1File.matchesSpace(spgrFile))) {
		cerr << "Mask or B1 dimensions/transform do not match SPGR file." << endl;
		exit(EXIT_FAILURE);
	}
	nSPGR = spgrFile.dim(4);
	VectorXd spgrAngles(nSPGR);
	
	#ifdef HAVE_NRECON
	ParameterList pars;
	if (ReadProcpar(spgrFile.basePath() + ".procpar", pars)) {
		spgrTR = RealValue(pars, "tr");
		for (int i = 0; i < nSPGR; i++) spgrAngles[i] = RealValue(pars, "flip1", i);
	} else
	#endif
	{
		cout << "Enter SPGR TR (seconds):"; cin >> spgrTR;
		cout << "Enter SPGR " << nSPGR << " Flip Angles (degrees):";
		for (int i = 0; i < nSPGR; i++) cin >> spgrAngles[i];
	}
	spgrAngles *= M_PI / 180.;
	//**************************************************************************
	#pragma mark Select which angles to use in the analysis
	//**************************************************************************	
	VectorXi spgrKeep(nSPGR);
	spgrKeep.setOnes();
	if (drop) {
		cout << "Choose SPGR angles to use (1 to keep, 0 to drop, " << nSPGR << " values): ";
		for (int i = 0; i < nSPGR; i++) cin >> spgrKeep[i];
		VectorXd temp = spgrAngles;
		spgrAngles.resize(spgrKeep.sum());
		int angle = 0;
		for (int i = 0; i < nSPGR; i++)
			if (spgrKeep(i)) spgrAngles(angle++) = temp(i);
	}
	if (verbose)
	{
		cout << "SPGR TR=" << spgrTR
		          << " s. Flip-angles: " << spgrAngles.transpose() * 180. / M_PI << endl;
		cout << "Ouput prefix will be: " << outPrefix << endl;
	}
	//**************************************************************************
	// Allocate memory for slices
	//**************************************************************************	
	int voxelsPerSlice = spgrFile.voxelsPerSlice();
	int voxelsPerVolume = spgrFile.voxelsPerVolume();
	
	cout << "Reading SPGR data..." << flush;
	vector<double> spgrData = spgrFile.readAllVolumes<double>();
	spgrFile.close();
	cout << "done." << endl;
	//**************************************************************************
	// Create results data storage
	//**************************************************************************
	#define NR 3
	vector<vector<double>> resultsData(NR);
	for (int i = 0; i < NR; i++)
		resultsData[i].resize(voxelsPerVolume);
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	for (int slice = 0; slice < spgrFile.dim(3); slice++) {
		clock_t loopStart;
		// Read in data
		if (verbose)
			cout << "Starting slice " << slice << "..." << flush;
		loopStart = clock();
		atomic<int> voxCount{0};
		int sliceOffset = slice * voxelsPerSlice;
		
		function<void (const int&)> processVox = [&] (const int &vox) {
			double T1 = 0., M0 = 0., B1 = 1., res = 0.; // Place to restore per-voxel return values, assume B1 field is uniform for classic DESPOT
			if (!maskFile.isOpen() || (maskData[sliceOffset + vox] > 0.))
			{
				voxCount++;
				if (B1File.isOpen())
					B1 = B1Data[sliceOffset + vox];
				ArrayXd spgrs(nSPGR);
				int vol = 0;
				for (int img = 0; img < nSPGR; img++) {
					if (spgrKeep(img))
						spgrs[vol++] = spgrData[img * voxelsPerVolume + sliceOffset + vox];
				}
				res = classicDESPOT1(spgrAngles, spgrs, spgrTR, B1, M0, T1);
			}
			resultsData[0][sliceOffset + vox] = M0;
			resultsData[1][sliceOffset + vox] = T1;
			resultsData[2][sliceOffset + vox] = res;
		};
		apply_for(voxelsPerSlice, processVox);
		
		if (verbose) {
			clock_t loopEnd = clock();
			if (voxCount > 0)
				cout << voxCount << " unmasked voxels, CPU time per voxel was "
				          << ((loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC)) << " s, ";
			cout << "finished." << endl;
		}
	}
	const string names[NR] = { "D1_PD", "D1_T1", "D1_Residual" };
	Nifti::File outFile(spgrFile);
	outFile.setDim(4, 1);
	outFile.setDatatype(DT_FLOAT32);
	for (int r = 0; r < NR; r++) {
		string outName = outPrefix + names[r] + ".nii.gz";
		if (verbose)
			cout << "Writing result header: " << outName << endl;
		outFile.open(outName, 'w');
		outFile.writeVolume<double>(0, resultsData[r]);
		outFile.close();
	}
	cout << "All done." << endl;
	exit(EXIT_SUCCESS);
}
