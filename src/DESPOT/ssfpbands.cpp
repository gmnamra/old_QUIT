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
	--help, -h        : Print this message.\n\
	--verbose, -v     : Print more information.\n\
	--out, -o path    : Specify an output filename (default image base).\n\
	--mask, -m file   : Mask input with specified file.\n\
	--flip, -f        : Data order is flip-angle, then phase (default opposite).\n\
	--phases, -p N    : Number of phase-cycling patterns used (default is 4).\n\
	--threads, -T N   : Use N threads (default=hardware limit).\n\
	--save, -sR       : Save the robustly regularised GS (default)\n\
	          M       : Save the magnitude regularised GS\n\
	          G       : Save the unregularised GS\n\
	          C       : Save the CS\n"
};

enum class Save { RR, MR, GS, CS };
static Save mode = Save::RR;
static bool verbose = false;
static size_t phase_dim = 3, flip_dim = 4, nPhases = 4;
static string prefix;
static struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"flip", required_argument, 0, 'f'},
	{"phases", required_argument, 0, 'p'},
	{"threads", required_argument, 0, 'T'},
	{"save", required_argument, 0, 's'},
	{0, 0, 0, 0}
};

// From Knuth, surprised this isn't in STL
unsigned long long choose(unsigned long long n, unsigned long long k) {
	if (k > n)
		return 0;

	unsigned long long r = 1;
	for (unsigned long long d = 1; d <= k; ++d) {
		r *= n--;
		r /= d;
	}
	return r;
}

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	Nifti::File maskFile;
	MultiArray<int8_t, 3> maskData;
	ThreadPool threads;

	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvo:m:fs:p:T:", long_options, &indexptr)) != -1) {
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
			case 'p':
				nPhases = atoi(optarg);
				if ((nPhases % 2) != 0) {
					cerr << "Number of phase-cycling patterns must be divisible by 2." << endl;
					return EXIT_FAILURE;
				}
				if (nPhases < 4) {
					cerr << "Must have a minimum of 4 phase-cycling patterns." << endl;
					return EXIT_FAILURE;
				}
				break;
			case 's':
				switch(*optarg) {
					case 'R': mode = Save::RR; break;
					case 'M': mode = Save::MR; break;
					case 'G': mode = Save::GS; break;
					case 'C': mode = Save::CS; break;
					default:
						cerr << "Unknown desired save image: " << *optarg << endl;
						return EXIT_FAILURE;
				}
				break;
			case 'T':
				threads.resize(atoi(optarg));
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
	const auto d = inputFile.matrix();
	size_t nFlip = inputFile.dim(4) / nPhases;
	if (nFlip < 1) {
		cerr << "The specified number of phase-cycling patterns is inconsistent with the number of volumes in the input file." << endl;
		return EXIT_FAILURE;
	}
	size_t nLines = nPhases / 2;
	size_t nCrossings = choose(nLines, 2);

	if (verbose) cout << "Reading data..." << endl;
	MultiArray<complex<float>, 4> input(inputFile.dims().head(4));
	inputFile.readVolumes(input.begin(), input.end());
	if (prefix == "") {
		prefix = inputFile.basePath() + "_";
	}
	inputFile.close();

	if (verbose) cout << "Prepping data..." << endl;
	typedef typename MultiArray<complex<float>, 5>::Index idx_t;
	idx_t reshape_dims; reshape_dims << d, 0, 0;
	reshape_dims[phase_dim] = 4;
	reshape_dims[flip_dim] = nFlip;
	MultiArray<complex<float>, 5> reshaped = input.reshape<5>(reshape_dims);
	idx_t split_start = idx_t::Zero();
	reshape_dims[phase_dim] = 2;
	MultiArray<complex<float>, 5> aData = reshaped.slice<5>(split_start, reshape_dims);
	split_start[phase_dim] = 2;
	MultiArray<complex<float>, 5> bData = reshaped.slice<5>(split_start, reshape_dims);
	// For results
	MultiArray<complex<float>, 4> output(d, nFlip);
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	clock_t startClock = clock();
	for (size_t vol = nFlip; vol-- > 0;) { // Reverse iterate over volumes
		if (verbose) cout << "Processing volume " << vol << "..." << endl;
		for (size_t vk = 0; vk < d[2]; vk++) {
			function<void (const size_t, const size_t)> processVox = [&] (const size_t vi, const size_t vj) {
				MultiArray<complex<float>, 5>::Index idx;
				idx << vi,vj,vk,0,0;
				idx[flip_dim] = vol;

				MatrixXf sols(2, nCrossings); sols.setZero();
				size_t si = 0;
				for (size_t li = 0; li < nLines; li++) {
					idx_t idx_i; idx_i << vi,vj,vk,0,0;
					idx_i[flip_dim] = vol;
					idx_i[phase_dim] = li;
					for (size_t lj = li + 1; lj < nLines; lj++) {
						idx_t idx_j = idx_i;
						idx_j[phase_dim] = lj;

						Vector2f a_i{aData[idx_i].real(), aData[idx_i].imag()};
						Vector2f a_j{aData[idx_j].real(), aData[idx_j].imag()};
						Vector2f b_i{bData[idx_i].real(), bData[idx_i].imag()};
						Vector2f b_j{bData[idx_j].real(), bData[idx_j].imag()};

						Vector2f d_i = (b_i - a_i);
						Vector2f d_j = (b_j - a_j);
						Vector2f n_i{d_i[1], -d_i[0]};
						Vector2f n_j{d_j[1], -d_j[0]};

						float mu = (a_j - a_i).dot(n_j) / d_i.dot(n_j);
						float nu = (a_i - a_j).dot(n_i) / d_j.dot(n_i);
						float xi = 1.0 - pow(d_i.dot(d_j) / (d_i.norm() * d_j.norm()),2.0);

						Vector2f cs = (a_i + a_j + b_i + b_j) / 4.0;
						Vector2f gs = a_i + mu * d_i;

						Vector2f rs = cs;
						if (vol < (nFlip - 1)) { // Use the phase of the last flip-angle for regularisation
							float phase = arg(output[{vi,vj,vk,nFlip-1}]);
							Vector2f d_p{cos(phase),sin(phase)};
							float lm_i = (a_i).dot(n_i) / d_p.dot(n_i);
							float lm_j = (a_j).dot(n_j) / d_p.dot(n_j);
							Vector2f p_i = lm_i * d_p;
							Vector2f p_j = lm_j * d_p;
							rs = (p_i + p_j) / 2.0;
						}

						bool rob_reg = true;
						// Do the logic this way round so NaN does not propagate
						if ((mu > -xi) && (mu < 1 + xi) && (nu > -xi) && (nu < 1 + xi))
							rob_reg = false;

						float gs_norm = gs.norm();
						bool mag_reg = true;
						if ((gs_norm < a_i.norm()) &&
						    (gs_norm < a_j.norm()) &&
						    (gs_norm < b_i.norm()) &&
						    (gs_norm < b_j.norm())) {
						    mag_reg = false;
						}
						switch (mode) {
							case Save::RR: sols.col(si++) = rob_reg ? rs : gs; break;
							case Save::MR: sols.col(si++) = mag_reg ? cs : gs; break;
							case Save::GS: sols.col(si++) = gs; break;
							case Save::CS: sols.col(si++) = cs; break;
						}
					}
				}
				Vector2f mean_sol = sols.rowwise().mean();
				output[{vi,vj,vk,vol}] = {mean_sol[0], mean_sol[1]};
			};
			threads.for_loop2(processVox, d[0], d[1]);
		}
	}
	printElapsedClock(startClock, d.prod());
	inHdr.setDim(4, nFlip);
	inHdr.setDatatype(Nifti::DataType::COMPLEX64);

	string outname = prefix;
	switch (mode) {
		case Save::RR: outname.append("reg" + OutExt()); break;
		case Save::MR: outname.append("magreg" + OutExt()); break;
		case Save::GS: outname.append("gs" + OutExt()); break;
		case Save::CS: outname.append("cs" + OutExt()); break;
	}
	if (verbose) cout << "Output filename: " << outname << endl;
	Nifti::File outFile(inHdr, outname);
	outFile.writeVolumes(output.begin(), output.end());
	outFile.close();
	if (verbose) cout << "Finished." << endl;
	return EXIT_SUCCESS;
}
