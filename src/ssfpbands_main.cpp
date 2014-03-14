/*
 *  ssfpbands_main.cpp
 *
 *  Created by Tobias Wood on 14/03/2014.
 *  Copyright (c) 2014 Tobias Wood.
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
#include "ThreadPool.h"
#include "DESPOT.h"

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
"Usage is: ssfpbands [options] input1 [input2] \n\
\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--mask, -m file   : Mask input with specified file\n\
	--type, -t p      : Input is magnitude and phase (default)\n\
	           i      : Input is real/imaginary\n\
			   c      : Input is complex\n"
};

enum class Type { Phase, Imag, Complex };
static bool verbose = false;
static string outPrefix;
static Type inputType = Type::Phase;
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"type", required_argument, 0, 't'},
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
	cout << version << endl << credit_me << endl;
	Nifti maskFile;
	Volume<bool> maskVol;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvo:m:t:", long_options, &indexptr)) != -1) {
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
			case 'i':
				switch (*optarg) {
					case 'p': inputType = Type::Phase;  cout << "Input is magnitude and phase." << endl; break;
					case 'i': inputType = Type::Imag; cout << "Input is real and imaginary." << endl; break;
					case 'c': inputType = Type::Complex; cout << "Input is complex." << endl; break;
					default:
						cout << "Unknown input type " << optarg << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 'h':
			case '?': // getopt will print an error message
				exit(EXIT_FAILURE);
		}
	}
	//**************************************************************************
	#pragma mark Gather data
	//**************************************************************************
	cout << "Opening input file: " << argv[optind] << endl;
	Nifti inputFile;
	inputFile.open(argv[optind++], Nifti::Mode::Read);
	if (maskFile.isOpen() && !maskFile.matchesSpace(inputFile)) {
		cerr << "Mask does not match input file." << endl;
		exit(EXIT_FAILURE);
	}
	Nifti templateFile(inputFile, 1);
	if ((inputFile.rank() < 4) || ((inputFile.dim(4) % 4) != 0)) {
		cout << "Input must contain 4 phase-cycles (0, 90, 180, 270)." << endl;
		exit(EXIT_FAILURE);
	}
	size_t nFlip = inputFile.dim(4) / 4;
	
	VolumeSeries<float> input1, input2;
	VolumeSeries<complex<float>> inputC;
	if (inputType != Type::Complex) {
		input1.readFrom(inputFile);
		inputFile.close();
		inputFile.open(argv[optind++], Nifti::Mode::Read);
		if (!inputFile.matchesSpace(templateFile)) {
			cerr << "Input files do not match." << endl;
			exit(EXIT_FAILURE);
		}
		input2.readFrom(inputFile);
		inputFile.close();
	} else {
		inputC.readFrom(inputFile);
		inputFile.close();
	}
	if ((argc - optind) != 0) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		exit(EXIT_FAILURE);
	}
	// Results storage
	VolumeSeries<float> noBands(input1.dims().head(3), nFlip);
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	ThreadPool pool;
	for (size_t k = 0; k < templateFile.dim(3); k++) {
		clock_t loopStart;
		// Read in data
		if (verbose)
			cout << "Starting slice " << k << "..." << flush;
		loopStart = clock();
		atomic<int> voxCount{0};
		
		for (size_t j = 0; j < templateFile.dim(2); j++) {
			function<void (const size_t)> processVox = [&] (const size_t i) {
				const typename Volume<float>::IndexArray vox{i, j, k};
				if (!maskFile.isOpen() || (maskVol[vox])) {
					voxCount++;

					ArrayXcd I1(nFlip), I2(nFlip), I3(nFlip), I4(nFlip);
					switch (inputType) {
						case (Type::Phase): {
							ArrayXd mag = input1.series(vox).cast<double>();
							ArrayXd ph  = input2.series(vox).cast<double>();
							
							I1.real() = mag.head(nFlip) * ph.head(nFlip).cos();
							I1.imag() = mag.head(nFlip) * ph.head(nFlip).sin();
							I1.real() = mag.segment(nFlip, nFlip) * ph.segment(nFlip, nFlip).cos();
							I1.imag() = mag.segment(nFlip, nFlip) * ph.segment(nFlip, nFlip).sin();
							I1.real() = mag.segment(2*nFlip, nFlip) * ph.segment(2*nFlip, nFlip).cos();
							I1.imag() = mag.segment(2*nFlip, nFlip) * ph.segment(2*nFlip, nFlip).sin();
							I1.real() = mag.tail(nFlip) * ph.tail(nFlip).cos();
							I1.imag() = mag.tail(nFlip) * ph.tail(nFlip).sin();
						}	break;
						case (Type::Imag): {
							ArrayXd re = input1.series(vox).cast<double>();
							ArrayXd im = input2.series(vox).cast<double>();
							
							I1.real() = re.head(nFlip);
							I1.imag() = im.head(nFlip);
							I1.real() = re.segment(nFlip, nFlip);
							I1.imag() = im.segment(nFlip, nFlip);
							I1.real() = re.segment(2*nFlip, nFlip);
							I1.imag() = im.segment(2*nFlip, nFlip);
							I1.real() = re.tail(nFlip);
							I1.imag() = im.tail(nFlip);
						 }	break;
						case (Type::Complex): {
							ArrayXcd input = inputC.series(vox).cast<complex<double>>();
							I1 = input.head(nFlip);
							I2 = input.segment(nFlip, nFlip);
							I3 = input.segment(2*nFlip, nFlip);
							I4 = input.tail(nFlip);
						}	break;
					}
					
					ArrayXcd l1 = I3 - I1;
					ArrayXcd l2 = I4 - I2;
					
					ArrayXcd crossPoint = (I2 - I1) / (l2 - l1);
					
					noBands.series(vox) = crossPoint.abs().cast<float>();
				}
			};
			
			pool.for_loop(processVox, templateFile.dim(1));
		}
		
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
	templateFile.setDim(4, nFlip);
	templateFile.open("no_bands.nii.gz", Nifti::Mode::Write);
	noBands.writeTo(templateFile);
	templateFile.close();
	cout << "All done." << endl;
	exit(EXIT_SUCCESS);
}
