/*
 *  despot_hifi.cpp
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright 2011 Tobias Wood. All rights reserved.
 *
 *  Based on code by Sean Deoni
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
const string credit {
"despot-hifi - Written by tobias.wood@kcl.ac.uk, based on work by Sean Deoni. \n\
Acknowledgements greatfully received, grant discussions welcome."
};

const string usage {
"Usage is: despot-hifi [options] spgr_input ir-spgr_input\n\
\
Options:\n\
	-m, --mask file  : Mask input with specified file.\n\
	--out, -o path    : Add a prefix to the output filenames.\n\
	-v, --verbose    : Print out more messages.\n\
	-i, --inv 0-3    : Specify the scanner Inversion mode:\n\
	                   0 (Default) Use raw segment TR from input file\n\
			           1 = 1.5T scanner, readout pulses div 2 + 2\n\
	                   2 = 3T, scale TI by 0.9, readout pulses div 2 + 2\n\
	                   3 = 3T, scale TI by 0.84, readout pulses + 2)\n\
			           1,2 & 3 MUST be followed by the PE Step Count\n\
	-p, --pe N       : PE Step Count.\n"
};

static int verbose = false, inversionMode = 0, peReadout = 0;
static string outPrefix;
static double inversionEfficiency = 0.;
static struct option long_options[] =
{
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"verbose", no_argument, 0, 'v'},
	{"inv", required_argument, 0, 'i'},
	{"pe", required_argument, 0, 'p'},
	{0, 0, 0, 0}
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	cout << credit << endl;
	int nSPGR = 0, nIR = 0;
	NiftiImage inFile, spgrFile, irFile;
	double *maskData = NULL;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "i:m:o:vp:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'i':
				inversionMode = atoi(optarg);
				break;
			case 'p':
				peReadout = atoi(optarg);
				break;
			case 'm':
				cout << "Opening mask file: " << optarg << endl;
				inFile.open(optarg, 'r');
				maskData = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'v': verbose = true; break;
			case '?': // getopt will print an error message
				exit(EXIT_FAILURE);
		}
	}
	if ((argc - optind) != 2) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		exit(EXIT_FAILURE);
	}
	if (inversionMode && !peReadout) {
		cerr << "Inversion mode was specified but not readout pulse count." << endl;
		exit(EXIT_FAILURE);
	}
	
	double TIScale = 0.;
	switch (inversionMode) {
		case 1:
			TIScale = 1.0;
			peReadout = (peReadout / 2) + 2;
			inversionEfficiency = 0.97;
			break;
		case 2:
			TIScale = 0.9; // From Sean's code
			peReadout = (peReadout / 2) + 2;
			inversionEfficiency = 0.97;
			break;
		case 3:
			TIScale = 0.84; // From Sean's code
			peReadout = peReadout + 2;
			inversionEfficiency = 0.97;
			break;
		case 0:
			TIScale = 1.0;
			peReadout = 0;
			break;
		default:
			fprintf(stderr, "Inversion mode must be 0 - 3\n");
			exit(EXIT_FAILURE);
			break;
	}

	//**************************************************************************
	#pragma mark Gather SPGR data
	//**************************************************************************
	cout << "Opening SPGR file: " << argv[optind] << endl;
	spgrFile.open(argv[optind], 'r');
	nSPGR = spgrFile.dim(4);
	VectorXd spgrAngles(nSPGR);
	double spgrTR;
	
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
	
	//**************************************************************************
	#pragma mark Gather IR-SPGR data
	//**************************************************************************	
	cout << "Opening IR-SPGR file: " << argv[++optind] << endl;
	irFile.open(argv[optind], 'r');
	irFile.checkVoxelsCompatible(spgrFile);
	nIR = irFile.dim(4);
	VectorXd irTI(nIR);
	double irAngle, irTR;
	
	#ifdef USE_PROCPAR
	if (ReadProcpar(irFile.basename() + ".procpar", pars)) {
		irAngle = RealValue(pars, "flip1") * M_PI / 180.;
		for (int i = 0; i < nIR; i++) irTI[i] = RealValue(pars, "ti", i);
		irTR = RealValue(pars, "trseg") - irTI[0];
	} else
	#endif
	{
		cout << "Enter IR-SPGR Flip Angle (degrees):"; cin >> irAngle; irAngle *= M_PI / 180.;
		if (peReadout > 0) {
			cout << "Enter IR-SPGR TR (s):"; cin >> irTR;
			irTR = irTR * peReadout;
			cout << "Enter IR-SPGR TI times (s):";
			for (int i = 0; i < nIR; i++) {
				cin >> irTI[i];
				irTI[i] *= TIScale;
			}
		} else {
			cout << "Enter IR-SPGR TI times (s):";
			for (int i = 0; i < nIR; i++) cin >> irTI[i];
			fprintf(stdout, "Enter first scan Segment TR (s):"); cin >> irTR;
			irTR -= irTI[0]; // Subtract off TI to get 
		}
	}
	if (verbose) {
		cout << "Found " << nIR << " SPGR-IR images with flip angle " << irAngle * 180. / M_PI << " degrees." << endl;
		cout << "Segment TR is " << irTR << " seconds." << endl;
		cout << "Inversion time(s) are ";
		for (int i = 0; i < nIR; i++) cout << irTI[i] << " ";
		cout << "seconds." << endl;
	}
	//**************************************************************************
	// Allocate memory for slices
	//**************************************************************************	
	int voxelsPerSlice = spgrFile.voxelsPerSlice();
	int voxelsPerVolume = spgrFile.voxelsPerVolume();
	
	cout << "Reading image data..." << flush;
	double *SPGR = spgrFile.readAllVolumes<double>();
	double *IR   = irFile.readAllVolumes<double>();
	spgrFile.close();
	irFile.close();
	cout << "done." << endl;
	//**************************************************************************
	// Create results data storage
	//**************************************************************************
	#define NR 4
	double **resultsData   = new double*[NR];
	for (int i = 0; i < NR; i++)
		resultsData[i] = new double[voxelsPerVolume];
	const string names[NR] = { "HIFI_M0", "HIFI_T1", "HIFI_B1", "HIFI_residual" };
	
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
			double T1 = 0., M0 = 0., B1 = 1., res = 0.; // Assume B1 field is uniform for classic DESPOT
			if ((!maskData) || (maskData[sliceOffset + vox] > 0.))
			{
				voxCount++;
				ArrayXd spgrs(nSPGR), irs(nIR);
				int vol = 0;
				for (int img = 0; img < nSPGR; img++)
						spgrs[vol++] = SPGR[img * voxelsPerVolume + sliceOffset + vox];
				for (int img = 0; img < nIR; img++)
						irs[img] = IR[img * voxelsPerVolume + sliceOffset + vox];
				res = calcHIFI(spgrAngles, spgrs, spgrTR,
				               irTI, irs, irAngle, irTR, peReadout, inversionEfficiency,
							   M0, T1, B1);
				// Sanity check
				M0 = clamp(M0, 0., 1.e7);
				T1 = clamp(T1, 0., 15.);
				B1 = clamp(B1, 0., 2.);
			}
			resultsData[0][sliceOffset + vox] = M0;
			resultsData[1][sliceOffset + vox] = T1;
			resultsData[2][sliceOffset + vox] = B1;
			resultsData[3][sliceOffset + vox] = res;
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
	
	//**************************************************************************
	#pragma mark Write out data
	//**************************************************************************
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
	delete[] IR;
	delete[] maskData;
	cout << "All done." << endl;
	exit(EXIT_SUCCESS);
}
