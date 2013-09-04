/*
 *  despot1hifi.cpp
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright (c) 2011-2013 Tobias Wood.
 *
 *  Based on code by Sean Deoni
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
#include "ThreadPool.h"

#ifdef AGILENT
	#include "procpar.h"
#endif

using namespace std;
using namespace Eigen;

//******************************************************************************
#pragma mark HIFI Calculation functions
//******************************************************************************
double HIFIResidual(const ArrayXd &flipAngles, const ArrayXd &spgrVals, const double spgrTR,
				const ArrayXd &TI, const ArrayXd &irVals, const double irFlipAngle,
				const double irTR, const double nReadout, const double eff,
                double &M0, double &T1, double &B1);
double HIFIResidual(const ArrayXd &flipAngles, const ArrayXd &spgrVals, const double spgrTR,
				const ArrayXd &TI, const ArrayXd &irVals, const double irFlipAngle,
				const double irTR, const double nReadout, const double eff,
                double &M0, double &T1, double &B1) {
	ArrayXd st = SPGR(flipAngles, spgrTR, B1, M0, T1);
	ArrayXd it = IRSPGR(TI, irTR, B1, irFlipAngle, eff, M0, T1);
	double res = (spgrVals - st).square().sum() +
	              (irVals - it).square().sum();
	return res;
}

double calcHIFI(const ArrayXd &flipAngles, const ArrayXd &spgrVals, const double spgrTR,
				const ArrayXd &TI, const ArrayXd &irVals, const double irFlipAngle,
				const double irTR, const double nReadout, const double eff,
                double &M0, double &T1, double &B1);
double calcHIFI(const ArrayXd &flipAngles, const ArrayXd &spgrVals, const double spgrTR,
				const ArrayXd &TI, const ArrayXd &irVals, const double irFlipAngle,
				const double irTR, const double nReadout, const double eff,
                double &M0, double &T1, double &B1) {
	// Golden Section Search to find B1
	// From www.mae.wvu.edu/~smirnov/nr/c10-1.pdf
	double R = 0.61803399; // Golden ratio - 1
	double C = 1 - R;
	double precision = 0.001;
    
	// Set up initial bracket using some guesses
	double B1_0 = 0.3; double B1_3 = 1.8; double B1_1, B1_2;
	
	B1 = B1_0;
	classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
	double res1 = HIFIResidual(flipAngles, spgrVals, spgrTR, TI, irVals, irFlipAngle, irTR, nReadout, eff, M0, T1, B1);
	B1 = B1_3;
	classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
	double res2 = HIFIResidual(flipAngles, spgrVals, spgrTR, TI, irVals, irFlipAngle, irTR, nReadout, eff, M0, T1, B1);
	if (res1 < res2) {
		B1_1 = B1_0 + 0.2;
		B1_2 = B1_1 + C * (B1_3 - B1_1);
	} else {
		B1_2 = B1_3 - 0.2;
		B1_1 = B1_2 - C * (B1_2 - B1_0);
	}
	
	B1 = B1_1;
	classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
	res1 = HIFIResidual(flipAngles, spgrVals, spgrTR, TI, irVals, irFlipAngle, irTR, nReadout, eff, M0, T1, B1);
	B1 = B1_2;
	classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
	res2 = HIFIResidual(flipAngles, spgrVals, spgrTR, TI, irVals, irFlipAngle, irTR, nReadout, eff, M0, T1, B1);
	while ( fabs(B1_3 - B1_0) > precision * (fabs(B1_1) + fabs(B1_2))) {
		if (res2 < res1) {
			B1_0 = B1_1; B1_1 = B1_2;
			B1_2 = R * B1_1 + C * B1_3;
			res1 = res2;
			B1 = B1_2;
			classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
			res2 = HIFIResidual(flipAngles, spgrVals, spgrTR, TI, irVals, irFlipAngle, irTR, nReadout, eff, M0, T1, B1);
		} else {
			B1_3 = B1_2; B1_2 = B1_1;
			B1_1 = R * B1_2 + C * B1_0;
			res2 = res1;
			B1 = B1_1;
			classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
			res1 = HIFIResidual(flipAngles, spgrVals, spgrTR, TI, irVals, irFlipAngle, irTR, nReadout, eff, M0, T1, B1);
		}
	}
	// Best value for B1
	if (res1 < res2) {
		B1 = B1_1;
		return res1;
	} else {
		B1 = B1_2;
		return res2;
	}
	std::cout << "Finished Golden Ratio Search" << std::endl<< std::endl<< std::endl<< std::endl;
}

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
	                   0 = Use raw segment TR from input file\n\
			           1 = 1.5T scanner, readout pulses div 2 + 2\n\
	                   2 (Default) = 3T, scale TI by 0.9, readout pulses div 2 + 2\n\
	                   3 = 3T, scale TI by 0.84, readout pulses + 2)\n\
"
};

static int verbose = false, inversionMode = 2, peReadout = 0;
static string outPrefix;
static double inversionEfficiency = 0.;
static struct option long_options[] =
{
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"verbose", no_argument, 0, 'v'},
	{"inv", required_argument, 0, 'i'},
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
	Nifti::File maskFile, spgrFile, irFile;
	vector<double> maskData;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "i:m:o:vp:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'i':
				inversionMode = atoi(optarg);
				if ((inversionMode < 0) || (inversionMode > 3)) {
					cout << "Bad inversion mode (" << inversionMode << "). Must be 0-3" << endl;
					exit(EXIT_FAILURE);
				}
				break;
			case 'm':
				cout << "Opening mask file: " << optarg << endl;
				maskFile.open(optarg, Nifti::Modes::Read);
				maskData = maskFile.readVolume<double>(0);
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
	
	//**************************************************************************
	#pragma mark Gather SPGR data
	//**************************************************************************
	cout << "Opening SPGR file: " << argv[optind] << endl;
	spgrFile.open(argv[optind], Nifti::Modes::Read);
	if (maskFile.isOpen() && !maskFile.matchesSpace(spgrFile)) {
		cerr << "SPGR file dimensions or transform do not match mask." << endl;
		exit(EXIT_FAILURE);
	}
	nSPGR = spgrFile.dim(4);
	VectorXd spgrAngles(nSPGR);
	double spgrTR;
	
	#ifdef AGILENT
	Agilent::ProcPar pp;
	if (ReadPP(spgrFile, pp)) {
		spgrTR = pp.realValue("tr");
		for (int i = 0; i < nSPGR; i++) spgrAngles[i] = pp.realValue("flip1", i);
	} else
	#endif
	{
		cout << "Enter SPGR TR (seconds):"; cin >> spgrTR;
		cout << "Enter " << nSPGR << " Flip Angles (degrees):";
		for (int i = 0; i < nSPGR; i++) cin >> spgrAngles[i];
	}
	spgrAngles *= M_PI / 180.;
	
	//**************************************************************************
	#pragma mark Gather IR-SPGR data
	//**************************************************************************	
	cout << "Opening IR-SPGR file: " << argv[++optind] << endl;
	irFile.open(argv[optind], Nifti::Modes::Read);
	if (!irFile.matchesSpace(spgrFile)) {
		cerr << "Header of " << spgrFile.imagePath() << " does not match " << irFile.imagePath() << endl;
		exit(EXIT_FAILURE);
	}
	nIR = irFile.dim(4);
	VectorXd irTI(nIR);
	double irAngle, irTR;
	
	#ifdef AGILENT
	if (ReadPP(irFile, pp)) {
		irAngle = pp.realValue("flip1") * M_PI / 180.;
		for (int i = 0; i < nIR; i++) irTI[i] = pp.realValue("ti", i);
		irTR = pp.realValue("trseg") - irTI[0];
	} else
	#endif
	{
		cout << "Enter IR-SPGR Flip Angle (degrees):"; cin >> irAngle; irAngle *= M_PI / 180.;
		if (inversionMode > 0) {
			cout << "Enter IR-SPGR TR (seconds):"; cin >> irTR;
			cout << "Enter original number of slices (PE2):"; cin >> peReadout;
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
			}
			irTR = irTR * peReadout;
			cout << "Enter " << nIR << " IR-SPGR TI times (seconds):";
			for (int i = 0; i < nIR; i++) {
				cin >> irTI[i];
				irTI[i] *= TIScale;
			}
		} else {
			cout << "Enter " << nIR << " IR-SPGR TI times (seconds):";
			for (int i = 0; i < nIR; i++) cin >> irTI[i];
			fprintf(stdout, "Enter first scan Segment TR (seconds):"); cin >> irTR;
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
	vector<double> SPGR = spgrFile.readAllVolumes<double>();
	vector<double> IR   = irFile.readAllVolumes<double>();
	spgrFile.close();
	irFile.close();
	cout << "done." << endl;
	//**************************************************************************
	// Create results data storage
	//**************************************************************************
	#define NR 4
	vector<vector<double>> resultsData(NR);
	for (int i = 0; i < NR; i++)
		resultsData[i].resize(voxelsPerVolume);
	const string names[NR] = { "HIFI_M0", "HIFI_T1", "HIFI_B1", "HIFI_residual" };
	
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	ThreadPool pool;
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
			if (!maskFile.isOpen() || (maskData[sliceOffset + vox] > 0.))
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
		pool.for_loop(processVox, voxelsPerSlice);
		
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
	Nifti::File outFile(spgrFile);
	outFile.setDim(4, 1);
	outFile.setDatatype(DT_FLOAT32);
	for (int r = 0; r < NR; r++) {
		string outName = outPrefix + names[r] + ".nii.gz";
		if (verbose)
			cout << "Writing result header: " << outName << endl;
		outFile.open(outName, Nifti::Modes::Write);
		outFile.writeVolume<double>(0, resultsData[r]);
		outFile.close();
	}
	cout << "All done." << endl;
	exit(EXIT_SUCCESS);
}
