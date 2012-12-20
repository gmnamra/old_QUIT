/*
 *  despot1_main.cpp
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright 2011 Tobias Wood. All rights reserved.
 *
 */

#include <time.h>
#include <getopt.h>
#include <iostream>
#include <atomic>
#include <Eigen/Dense>

#include "NiftiImage.h"
#include "DESPOT.h"

#define USE_PROCPAR
#ifdef USE_PROCPAR
#include "procpar.h"
using namespace Recon;
#endif

using namespace std;
using namespace Eigen;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: despot1 [options] spgr_input output_prefix \n\
\
Options:\n\
	-m, --mask file : Mask input with specified file.\n\
	--B1 file       : Correct flip angles with specified B1 ratio.\n\
	-v, --verbose   : Print out more messages.\n\
	-d, --drop      : Drop certain flip-angles (Read from stdin).\n"
};

const string credit {
"despot1 - Written by tobias.wood@kcl.ac.uk, based on work by Sean Deoni. \n\
Acknowledgements greatfully received, grant discussions welcome."
};

static int verbose = false, drop = false;
static struct option long_options[] =
{
	{"B1", required_argument, 0, '1'},
	{"mask", required_argument, 0, 'm'},
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
	double *B1Data = NULL, *maskData = NULL;
	NiftiImage spgrFile, inFile;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "m:vd", long_options, &indexptr)) != -1) {
		switch (c) {
			case '1':
				cout << "Opening B1 file: " << optarg << endl;
				inFile.open(optarg, 'r');
				B1Data = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case 'm':
				cout << "Opening mask file: " << optarg << endl;
				inFile.open(optarg, 'r');
				maskData = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case 'v': verbose = true; break;
			case 'd': drop = true; break;
			case '?': // getopt will print an error message
				exit(EXIT_FAILURE);
		}
	}
	if ((argc - optind) != 2) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		exit(EXIT_FAILURE);
	}
	
	//**************************************************************************
	#pragma mark Gather SPGR data
	//**************************************************************************
	cout << "Opening SPGR file: " << argv[optind] << endl;
	spgrFile.open(argv[optind], 'r');
	nSPGR = spgrFile.dim(4);
	VectorXd spgrAngles(nSPGR);
	
	#ifdef USE_PROCPAR
	ParameterList pars;
	if (ReadProcpar(spgrFile.basename() + ".procpar", pars)) {
		spgrTR = RealValue(pars, "tr");
		for (int i = 0; i < nSPGR; i++) spgrAngles[i] = RealValue(pars, "flip1", i);
	} else
	#endif
	{
		cout << "Enter SPGR TR (s):"; cin >> spgrTR;
		cout << "Enter SPGR Flip Angles (degrees):";
		for (int i = 0; i < nSPGR; i++) cin >> spgrAngles[i];
	}
	spgrAngles *= M_PI / 180.;
	const string outPrefix(argv[++optind]);
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
	int totalVoxels = spgrFile.voxelsPerVolume();
	
	cout << "Reading SPGR data..." << flush;
	double *SPGR = spgrFile.readAllVolumes<double>();
	spgrFile.close();
	cout << "done." << endl;
	//**************************************************************************
	// Create results data storage
	//**************************************************************************
	#define NR 3
	double **resultsData   = new double*[NR];
	for (int i = 0; i < NR; i++)
		resultsData[i] = new double[totalVoxels];
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	for (int slice = 0; slice < spgrFile.dim(3); slice++)
	{
		clock_t loopStart;
		// Read in data
		if (verbose)
			cout << "Starting slice " << slice << "..." << flush;
		loopStart = clock();
		atomic<int> voxCount{0};
		int sliceOffset = slice * voxelsPerSlice;
		
		function<void (const int&)> processVox = [&] (const int &vox) {
			double T1 = 0., M0 = 0., B1 = 1., res = 0.; // Place to restore per-voxel return values, assume B1 field is uniform for classic DESPOT
			if ((!maskData) || (maskData[sliceOffset + vox] > 0.))
			{
				voxCount++;
				if (B1Data)
					B1 = B1Data[sliceOffset + vox];
				ArrayXd spgrs(nSPGR);
				int vol = 0;
				for (int img = 0; img < nSPGR; img++) {
					if (spgrKeep(img))
						spgrs[vol++] = SPGR[img * totalVoxels + sliceOffset + vox];
				}
				res = classicDESPOT1(spgrAngles, spgrs, spgrTR, B1, M0, T1);
				
				// Sanity check
				M0 = clamp(M0, 0., 1.e7);
				T1 = clamp(T1, 0., 15.);
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
	const string names[NR] = { "_M0", "_T1", "_despot1_res" };
	NiftiImage outFile(spgrFile);
	outFile.setDatatype(DT_FLOAT32);
	outFile.setDim(4, 1);
	for (int r = 0; r < NR; r++)
	{
		string outName = outPrefix + names[r] + ".nii.gz";
		if (verbose)
			cout << "Writing result header: " << outName << endl;
		outFile.open(outName, 'w');
		outFile.writeVolume<double>(0, resultsData[r]);
		outFile.close();
		delete[] resultsData[r];
	}
	// Clean up memory
	delete[] resultsData;
	delete[] SPGR;
	delete[] B1Data;
	delete[] maskData;
	cout << "All done." << endl;
	exit(EXIT_SUCCESS);
}
