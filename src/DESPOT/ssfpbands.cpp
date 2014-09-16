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
using namespace QUIT;

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
	--flip, -f        : Data order is flip-angle, then phase (default opposite)\n\
	--save, -s r      : Save the line-regularised image (default)\n\
	           m      : Save the magnitude-regularised image\n\
	           e      : Save the ellipse cross-point\n\
	           c      : Save the complex sum/average\n"
};

enum class SaveMode { LineReg, MagReg, CrossPoint, ComplexSum, Lambda, Mu };
static bool verbose = false;
static size_t phase_dim = 3, flip_dim = 4;
static string prefix;
static SaveMode save = SaveMode::LineReg;
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"flip", required_argument, 0, 'f'},
	{"save", required_argument, 0, 's'},
	{0, 0, 0, 0}
};
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	Nifti::File maskFile;
	MultiArray<int8_t, 3> maskData;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvo:m:fs:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskData.resize(maskFile.matrix());
				maskFile.readVolumes(maskData.begin(), maskData.end(), 0, 1);
				break;
			case 'o':
				prefix = optarg;
				cout << "Output prefix will be: " << prefix << endl;
				break;
			case 'f':
				phase_dim = 4; flip_dim = 3; break;
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
						return EXIT_FAILURE;
				}
				break;
			case 'h':
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
		}
	}
	if (verbose) cout << version << endl << credit_me << endl;
	if ((argc - optind) != 1) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		return EXIT_FAILURE;
	}
	if (verbose) cout << "Opening input file: " << argv[optind] << endl;
	string fname(argv[optind++]);
	Nifti::File inputFile(fname);
	Nifti::Header inHdr = inputFile.header();
	if (maskFile && !maskFile.header().matchesSpace(inHdr)) {
		cerr << "Mask does not match input file." << endl;
		return EXIT_FAILURE;
	}
	if ((inputFile.rank() < 4) || ((inputFile.dim(4) % 4) != 0)) {
		cout << "Input must contain 4 phase-cycles (0, 90, 180, 270)." << endl;
		return EXIT_FAILURE;
	}
	MultiArray<complex<float>, 4> input(inputFile.dims().head(4));
	inputFile.readVolumes(input.begin(), input.end());
	if (prefix == "") {
		prefix = inputFile.basePath();
	}
	inputFile.close();

	// Results storage
	const auto d = inputFile.matrix();
	size_t nFlip = inputFile.dim(4) / 4;
	MultiArray<complex<float>, 5>::Index nd; nd << d, 0, 0;
	nd[phase_dim] = 4;
	nd[flip_dim] = nFlip;
	MultiArray<complex<float>, 5> CData = input.reshape<5>(nd);
	MultiArray<float, 4> lm_out(d, nFlip), mu_out(d, nFlip), mabs_out(d, nFlip), leeway_out(d, nFlip);
	MultiArray<complex<float>, 4> cs_out(d, nFlip), gs_out(d, nFlip), reg_out(d, nFlip);
	MultiArray<int, 4> lreg_out(d, nFlip);
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	clock_t startClock = clock();
	ThreadPool pool;
	for (size_t vol = 0; vol < nFlip; vol++) {
		if (verbose) cout << "Processing volume " << vol << "..." << endl;
		function<void (const size_t)> processVox = [&] (const size_t k) {
			for (size_t j = 0; j < d[1]; j++) {
				for (size_t i = 0; i < d[0]; i++) {
					MultiArray<complex<float>, 5>::Index idx;
					idx << i,j,k,0,0;
					idx[flip_dim] = vol;

					complex<float> c_0   = CData[idx]; idx[phase_dim] = 1;
					complex<float> c_90  = CData[idx]; idx[phase_dim] = 2;
					complex<float> c_180 = CData[idx]; idx[phase_dim] = 3;
					complex<float> c_270 = CData[idx];

					Vector2f v_0{c_0.real(), c_0.imag()};
					Vector2f v_90{c_90.real(), c_90.imag()};
					Vector2f v_180{c_180.real(), c_180.imag()};
					Vector2f v_270{c_270.real(), c_270.imag()};

					Vector2f d_0 = v_180 - v_0;
					Vector2f d_90 = v_270 - v_90;
					Vector2f n_0(d_0[1], -d_0[0]);
					Vector2f n_90(d_90[1], -d_90[0]);

					complex<float> cs = (c_0 + c_90 + c_180 + c_270) / complex<float>(4.0,0.0);
					cs_out[{i,j,k,vol}] = cs;

					float lm_0 = ((v_90 - v_0).dot(n_90)) / (d_0.dot(n_90));
					float lm_90 = ((v_0 - v_90).dot(n_0)) / (d_90.dot(n_0));

					Vector2f gs_0 = (v_0 + lm_0 * d_0);
					Vector2f gs_90 = (v_90 + lm_90 * d_90);

					gs_out[{i,j,k,vol}] = {gs_0[0], gs_0[1]};
					lm_out[{i,j,k,vol}] = lm_0;
					mu_out[{i,j,k,vol}] = lm_90;

					float max_abs = max(max(max(abs(c_0),abs(c_90)),abs(c_180)),abs(c_270));
					bool mreg = gs_0.norm() > max_abs;
					bool lreg = false;

					float leeway = 1. - fabs(d_0.dot(d_90) / (d_0.norm() * d_90.norm()));


					if (lm_0 < -leeway) { lreg = true; lm_0 = 0; }
					if (lm_0 > (1. + leeway)) { lreg = true; lm_0 = 1; }
					if (lm_90 < -leeway) { lreg = true; lm_90 = 0; }
					if (lm_90 > (1. + leeway)) { lreg = true; lm_90 = 1; }
					gs_0 = (v_0 + lm_0 * d_0);
					gs_90 = (v_90 + lm_90 * d_90);

					mabs_out[{i,j,k,vol}] = max_abs;
					leeway_out[{i,j,k,vol}] = leeway;
					lreg_out[{i,j,k,vol}] = lreg;

					reg_out[{i,j,k,vol}] = lreg ? cs : gs_out[{i,j,k,vol}];
				}
			}
		};
		pool.for_loop(processVox, d[2]);
	}
	printElapsedClock(startClock, d.prod());
	inHdr.setDim(4, nFlip);
	inHdr.setDatatype(Nifti::DataType::COMPLEX64);

	string outname = prefix + "_lreg" + OutExt();
	Nifti::File outFile(inHdr, outname);
	outFile.writeVolumes(lreg_out.begin(), lreg_out.end());
	outFile.close();
	outname = prefix + "_lw" + OutExt();
	outFile.open(outname, Nifti::Mode::Write);
	outFile.writeVolumes(leeway_out.begin(), leeway_out.end());
	outFile.close();
	outname = prefix + "_mabs" + OutExt();
	outFile.open(outname, Nifti::Mode::Write);
	outFile.writeVolumes(mabs_out.begin(), mabs_out.end());
	outFile.close();
	outname = prefix + "_cs" + OutExt();
	outFile.open(outname, Nifti::Mode::Write);
	outFile.writeVolumes(cs_out.begin(), cs_out.end());
	outFile.close();
	outname = prefix + "_gs" + OutExt();
	outFile.open(outname, Nifti::Mode::Write);
	outFile.writeVolumes(gs_out.begin(), gs_out.end());
	outFile.close();
	outname = prefix + "_reg" + OutExt();
	outFile.open(outname, Nifti::Mode::Write);
	outFile.writeVolumes(reg_out.begin(), reg_out.end());
	outFile.close();
	inHdr.setDatatype(Nifti::DataType::FLOAT64);
	outname = prefix + "_lm" + OutExt();
	outFile.open(outname, Nifti::Mode::Write);
	outFile.writeVolumes(lm_out.begin(), lm_out.end());
	outFile.close();
	outname = prefix + "_mu" + OutExt();
	outFile.open(outname, Nifti::Mode::Write);
	outFile.writeVolumes(mu_out.begin(), mu_out.end());
	outFile.close();
	if (verbose) cout << "Finished." << endl;
	return EXIT_SUCCESS;
}
