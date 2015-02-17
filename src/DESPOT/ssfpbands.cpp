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
	--flip, -F        : Data order is phase, then flip-angle (default opposite).\n\
	--phases, -p N    : Number of phase-cycling patterns used (default is 4).\n\
	--alternate, -a   : Opposing phase-cycles alternate (default is two blocks).\n\
	--threads, -T N   : Use N threads (default=hardware limit).\n\
	--save, -sL       : Save the line regularised GS (default)\n\
	          M       : Save the magnitude regularised GS\n\
	          G       : Save the unregularised GS\n\
	          C       : Save the CS\n\
	          P       : Save the PS\n\
	--secondpass, -2  : Perform a 2nd pass as per Xiang and Hoff\n"
};

enum class Save { LR, MR, GS, CS, PS };
static Save mode = Save::LR;
static bool verbose = false, pass2 = false, pass22d = false, alternate = false;
static size_t phase_dim = 4, flip_dim = 3, nPhases = 4;
static string prefix;
const struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"flip", required_argument, 0, 'F'},
	{"phases", required_argument, 0, 'p'},
	{"alternate", no_argument, 0, 'a'},
	{"threads", required_argument, 0, 'T'},
	{"save", required_argument, 0, 's'},
	{"secondpass", optional_argument, 0, '2'},
	{0, 0, 0, 0}
};
const char *short_options = "hvo:m:Fs:p:aT:2::";
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
	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
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
			case 'F':
				phase_dim = 3; flip_dim = 4; break;
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
			case 'a': alternate = true; break;
			case 's':
				switch(*optarg) {
					case 'L': mode = Save::LR; break;
					case 'M': mode = Save::MR; break;
					case 'G': mode = Save::GS; break;
					case 'C': mode = Save::CS; break;
					case 'P': mode = Save::PS; break;
					default:
						cerr << "Unknown desired save image: " << *optarg << endl;
						return EXIT_FAILURE;
				}
				break;
			case '2': pass2 = true; if (optarg) pass22d = true; break;
			case 'T':
				threads.resize(atoi(optarg));
				break;
			case 'h':
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
		}
	}
	if (pass22d) cout << "2D second pass selected." << endl;
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
	reshape_dims[phase_dim] = nPhases;
	reshape_dims[flip_dim] = nFlip;
	MultiArray<complex<float>, 5> reshaped = input.reshape<5>(reshape_dims);
	idx_t split_start = idx_t::Zero();
	reshape_dims[phase_dim] = 2;
	idx_t split_strides = idx_t::Ones();
	if (alternate)
		split_strides[phase_dim] = 2;
	MultiArray<complex<float>, 5> aData = reshaped.slice<5>(split_start, reshape_dims, split_strides);
	if (alternate)
		split_start[phase_dim] = 1;
	else
		split_start[phase_dim] = nPhases / 2;
	MultiArray<complex<float>, 5> bData = reshaped.slice<5>(split_start, reshape_dims, split_strides);
	// For results
	MultiArray<complex<float>, 4> output(d, nFlip), second_pass(d, nFlip);
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	clock_t startClock = clock();
	for (size_t vol = nFlip; vol-- > 0;) { // Reverse iterate over volumes
		if (verbose) cout << "Processing volume " << vol << "..." << endl;
		for (size_t vk = 0; vk < d[2]; vk++) {
			function<void (const size_t, const size_t)> processVox = [&] (const size_t vi, const size_t vj) {
				MatrixXf sols(2, nCrossings); sols.setZero();
				if (!maskFile || (maskData[{vi,vj,vk}])) {
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

							Vector2f ps;
							if (vol < (nFlip - 1)) { // Use the phase of the last flip-angle for regularisation
								float phase = arg(output[{vi,vj,vk,nFlip-1}]);
								Vector2f d_p{cos(phase),sin(phase)};
								float lm_i = (a_i).dot(n_i) / d_p.dot(n_i);
								float lm_j = (a_j).dot(n_j) / d_p.dot(n_j);
								Vector2f p_i = lm_i * d_p;
								Vector2f p_j = lm_j * d_p;
								ps = (p_i + p_j) / 2.0;
							}

							bool line_reg = true;
							// Do the logic this way round so NaN does not propagate
							if ((mu > -xi) && (mu < 1 + xi) && (nu > -xi) && (nu < 1 + xi))
								line_reg = false;

							bool mag_reg = true;
							float maxnorm = max(max(max(a_i.norm(), a_j.norm()), b_i.norm()), b_j.norm());
							if (gs.norm() < maxnorm) {
								mag_reg = false;
							}
							if (ps.norm() > maxnorm) {
								ps = cs;
							}
							switch (mode) {
								case Save::LR: sols.col(si) = line_reg ? ps : gs; break;
								case Save::MR: sols.col(si) = mag_reg ? ps : gs; break;
								case Save::GS: sols.col(si) = gs; break;
								case Save::PS: sols.col(si) = ps; break;
								case Save::CS: sols.col(si) = cs; break;
							}
							si++;
						}
					}
				}
				Vector2f mean_sol = sols.rowwise().mean();
				output[{vi,vj,vk,vol}] = {mean_sol[0], mean_sol[1]};
			};
			threads.for_loop2(processVox, d[0], d[1]);
			if (threads.interrupted())
				break;
		}

		if (pass2) {
			for (size_t vk = 0; vk < d[2]; vk++) {
			//size_t vk = 27;
				function<void (const size_t, const size_t)> processVox = [&] (const size_t vi, const size_t vj) {
					complex<float> sp(0.,0.);
					if (!maskFile || (maskData[{vi,vj,vk}])) {
						for (size_t li = 0; li < nLines; li++) {
							float num = 0, den = 0;
							for (int k = ((!pass22d && (vk > 0)) ? -1 : 0); k < ((!pass22d && (vk < d[2] - 1)) ? 2 : 1); k++) {
								for (int j = (vj > 0 ? -1 : 0); j < (vj < d[1] - 1 ? 2 : 1); j++) {
									for (int i = (vi > 0 ? -1 : 0); i < (vi < d[0] - 1 ? 2 : 1); i++) {
										idx_t idx; idx << vi + i,vj + j,vk + k,0,0;
										idx[flip_dim] = vol;
										idx[phase_dim] = li;
										complex<float> a_i = aData[idx];
										complex<float> b_i = bData[idx];
										complex<float> s_i = output[{vi+i,vj+j,vk+k,vol}];
										num += real(conj(b_i - s_i)*(a_i - b_i) + conj(a_i - b_i)*(b_i - s_i));
										den += real(conj(a_i - b_i)*(a_i - b_i));
									}
								}
							}
							float w = -num / (2. * den);
							if (isfinite(w)) {
								idx_t sidx; sidx << vi,vj,vk,0,0;
								sidx[flip_dim] = vol;
								sidx[phase_dim] = li;
								complex<float> s = (aData[sidx]*w + (1.f-w)*bData[sidx]);
								sp += (s / static_cast<float>(nLines));
							}
						}
					}
					second_pass[{vi,vj,vk,vol}] = sp;
				};
				threads.for_loop2(processVox, 0, d[0], 1, 0, d[1], 1);
				if (threads.interrupted())
					break;
			}
		}
	}
	if (verbose) printElapsedClock(startClock, d.prod());
	inHdr.setDim(4, nFlip);
	inHdr.setDatatype(Nifti::DataType::COMPLEX64);

	string outname = prefix;
	switch (mode) {
		case Save::LR: outname.append("lreg"); break;
		case Save::MR: outname.append("mreg"); break;
		case Save::GS: outname.append("gs"); break;
		case Save::CS: outname.append("cs"); break;
		case Save::PS: outname.append("ps"); break;
	}
	if (pass2) outname.append("_2p");
	outname.append(OutExt());
	if (verbose) cout << "Output filename: " << outname << endl;
	Nifti::File outFile(inHdr, outname);
	if (pass2 == true) {
		outFile.writeVolumes(second_pass.begin(), second_pass.end());
	} else {
		outFile.writeVolumes(output.begin(), output.end());
	}
	outFile.close();
	if (verbose) cout << "Finished." << endl;
	return EXIT_SUCCESS;
}
