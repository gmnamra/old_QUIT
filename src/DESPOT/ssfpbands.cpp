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
#include "Eigen/Dense"

#include "Nifti/Nifti.h"
#include "QUIT/QUIT.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: ssfpbands [options] input \n\
\n\
Input must be a single complex image with 0, 90, 180, 360 phase-cycles in order\n\
along the 4th dimension.\n\
\n\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--out, -o path    : Specify an output filename (default image base)\n\
	--mask, -m file   : Mask input with specified file\n\
	--save, -s r      : Save the line-regularised image (default)\n\
	           m      : Save the magnitude-regularised image\n\
	           e      : Save the ellipse cross-point\n\
	           c      : Save the complex sum/average\n"
};

enum class SaveMode { LineReg, MagReg, CrossPoint, ComplexSum, Lambda, Mu };
static bool verbose = false;
static string outname;
static SaveMode save = SaveMode::LineReg;
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"save", required_argument, 0, 's'},
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
	MultiArray<int8_t, 3> maskData;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvo:m:s:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskData.resize(maskFile.dims().head(3));
				maskFile.readVolumes(maskData.begin(), maskData.end(), 0, 1);
				break;
			case 'o':
				outname = optarg;
				cout << "Output prefix will be: " << outname << endl;
				break;
			case 's':
				switch (*optarg) {
					case 'r': save = SaveMode::LineReg; break;
					case 'm': save = SaveMode::MagReg; break;
					case 'e': save = SaveMode::CrossPoint; break;
					case 'c': save = SaveMode::ComplexSum; break;
					case 'l': save = SaveMode::Lambda; break;
					case 'u': save = SaveMode::Mu; break;
					default:
						cout << "Unknown regularisation mode '" << *optarg << "'" << endl;
						exit(EXIT_FAILURE);
				}
				break;
			case 'h':
			case '?': // getopt will print an error message
				exit(EXIT_FAILURE);
		}
	}
	//**************************************************************************
	#pragma mark Gather data
	//**************************************************************************
	if ((argc - optind) != 1) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		exit(EXIT_FAILURE);
	}
	if (verbose) cout << "Opening input file: " << argv[optind] << endl;
	Nifti inputFile;
	inputFile.open(argv[optind++], Nifti::Mode::Read);
	if (maskFile.isOpen() && !maskFile.matchesSpace(inputFile)) {
		cerr << "Mask does not match input file." << endl;
		exit(EXIT_FAILURE);
	}
	if ((inputFile.rank() < 4) || ((inputFile.dim(4) % 4) != 0)) {
		cout << "Input must contain 4 phase-cycles (0, 90, 180, 270)." << endl;
		exit(EXIT_FAILURE);
	}
	MultiArray<complex<float>, 4> input(inputFile.dims());
	inputFile.readVolumes(input.begin(), input.end());
	if (outname == "") {
		outname = inputFile.basePath();
	}
	inputFile.close();

	// Results storage
	auto d = input.dims().head(3);
	size_t nFlip = input.dims()[3] / 4;
	MultiArray<complex<float>, 5>::Index nd; nd << d, 4, nFlip;
	auto CData = input.reshape<5>(nd);
	MultiArray<complex<float>, 4> output(d, nFlip);
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	ThreadPool pool;
	for (size_t vol = 0; vol < nFlip; vol++) {
		if (verbose) cout << "Processing volume " << vol << "..." << endl;
		function<void (const size_t)> processVox = [&] (const size_t k) {
			for (size_t j = 0; j < d[1]; j++) {
				decltype(CData)::Index idx; idx << 0,j,k,0,vol;
				decltype(CData)::Index sz; sz << -1,0,0,0,0;
				auto C1 = CData.slice<1>(idx, sz).asArray(); idx[3] = 1;
				auto C2 = CData.slice<1>(idx, sz).asArray(); idx[3] = 2;
				auto C3 = CData.slice<1>(idx, sz).asArray(); idx[3] = 3;
				auto C4 = CData.slice<1>(idx, sz).asArray();
				ArrayXXf I1(d[0], 2), I2(d[0], 2), I3(d[0], 2), I4(d[0], 2),
				         d1(d[0], 2), d2(d[0], 2), n1(d[0], 2), n2(d[0], 2);
				I1.col(0) = C1.real(); I1.col(1) = C1.imag();
				I2.col(0) = C2.real(); I2.col(1) = C2.imag();
				I3.col(0) = C3.real(); I3.col(1) = C3.imag();
				I4.col(0) = C4.real(); I4.col(1) = C4.imag();
				d1 = (I3 - I1);
				d2 = (I4 - I2);
				n1.col(0) = d1.col(1); n1.col(1) = -d1.col(0);
				n2.col(0) = d2.col(1); n2.col(1) = -d2.col(0);
				auto lm = (n2 * (I2 - I1)).rowwise().sum() / (n2 * d1).rowwise().sum();
				auto mu = (n1 * (I1 - I2)).rowwise().sum() / (n1 * d2).rowwise().sum();
				auto cs = (C1 + C2 + C3 + C4) / 4;
				auto cp2d = I1 + d1.colwise()*lm;
				ArrayXcf cp(d[0]); cp.real() = cp2d.col(0); cp.imag() = cp2d.col(1);
				auto regLine = ((lm > 0) && (lm < 1) && (mu > 0) && (mu < 1));
				auto regMag  = (cp.abs() < C1.abs()) && (cp.abs() < C2.abs()) &&
				               (cp.abs() < C3.abs()) && (cp.abs() < C4.abs());
				auto out = output.slice<1>({0,j,k,vol}, {-1,0,0,0}).asArray();
				switch (save) {
					case SaveMode::LineReg:    out = regLine.select(cp, cs); break;
					case SaveMode::MagReg:     out = regMag.select(cp, cs); break;
					case SaveMode::CrossPoint: out = cp; break;
					case SaveMode::ComplexSum: out = cs; break;
					case SaveMode::Lambda:     out = lm.cast<complex<float>>(); break;
					case SaveMode::Mu:         out = mu.cast<complex<float>>(); break;
				}
				if (maskFile.isOpen()) {
					auto m = maskData.slice<1>({0,j,k},{-1,0,0}).asArray();
					out = m.select(out, 0);
				}
			}
		};
		pool.for_loop(processVox, d[2]);
	}
	inputFile.setDim(4, nFlip);
	inputFile.setDatatype(Nifti::DataType::COMPLEX128);
	switch (save) {
		case SaveMode::LineReg:    outname += "_lreg.nii.gz"; break;
		case SaveMode::MagReg:     outname += "_mreg.nii.gz"; break;
		case SaveMode::CrossPoint: outname += "_cross.nii.gz"; break;
		case SaveMode::ComplexSum: outname += "_sum.nii.gz"; break;
		case SaveMode::Lambda:     outname += "_lambda.nii.gz"; inputFile.setDatatype(Nifti::DataType::FLOAT64); break;
		case SaveMode::Mu:         outname += "_mu.nii.gz";     inputFile.setDatatype(Nifti::DataType::FLOAT64); break;
	}
	if (verbose) cout << "Writing output file: " << outname << endl;
	inputFile.open(outname, Nifti::Mode::Write);
	inputFile.writeVolumes(output.begin(), output.end());
	inputFile.close();
	cout << "Finished." << endl;
	exit(EXIT_SUCCESS);
}
