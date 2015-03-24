/*
 *  despot1_main.cpp
 *
 *  Created by Tobias Wood on 27/01/2015.
 *  Copyright (c) 2015 Tobias Wood.
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
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "Nifti/Nifti.h"
#include "QUIT/QUIT.h"
#include "SignalEquations.h"

using namespace std;
using namespace Eigen;
using namespace QUIT;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: t2star [options] input_file \n\
\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--no-prompt, -n   : Suppress input prompts\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--mask, -m file   : Mask input with specified file\n\
	--echos, -e N     : Number of echos (timeseries input)\n\
	--thresh, -t n    : Threshold maps at PD < n\n\
	--clamp, -c n     : Clamp T2* between 0 and n\n\
	--algo, -a L      : LLS algorithm (default)\n\
	           A      : ARLO algorithm\n\
	--resids, -r      : Write out per flip-angle residuals\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};

enum class Algo { LogLin, ARLO };
static Algo algo = Algo::LogLin;
static int NE = 0;
static bool verbose = false, prompt = true, all_residuals = false;
static string outPrefix;
static double thresh = -numeric_limits<double>::infinity();
static double clamp_lo = -numeric_limits<double>::infinity(), clamp_hi = numeric_limits<double>::infinity();
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"echos", required_argument, 0, 'e'},
	{"thresh", required_argument, 0, 't'},
	{"clamp", required_argument, 0, 'c'},
	{"algo", required_argument, 0, 'a'},
	{"threads", required_argument, 0, 'T'},
	{"resids", no_argument, 0, 'r'},
	{0, 0, 0, 0}
};
static const char *short_opts = "hvnm:e:o:b:t:c:a:T:r";
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	try { // To fix uncaught exceptions on Mac
	cout << version << endl << credit_shared << endl;
	Eigen::initParallel();
	Nifti::File inputFile, maskFile;
	MultiArray<int8_t, 3> maskVol;
	ThreadPool threads;
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskVol.resize(maskFile.matrix());
				maskFile.readVolumes(maskVol.begin(), maskVol.end(), 0, 1);
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'e': NE = atoi(optarg); break;
			case 't': thresh = atof(optarg); break;
			case 'c':
				clamp_lo = 0;
				clamp_hi = atof(optarg);
				break;
			case 'a':
				switch (*optarg) {
					case 'L': algo = Algo::LogLin; if (verbose) cout << "LogLin algorithm selected." << endl; break;
					case 'A': algo = Algo::ARLO; if (verbose) cout << "ARLO algorithm selected." << endl; break;
					default:
						cout << "Unknown algorithm type " << optarg << endl;
						return EXIT_FAILURE;
						break;
				} break;
			case 'T':
				threads.resize(atoi(optarg));
				break;
			case 'r': all_residuals = true; break;
			case 'h':
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != 1) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		return EXIT_FAILURE;
	}
	// Gather input data
	cout << "Opening input file: " << argv[optind] << endl;
	inputFile.open(argv[optind], Nifti::Mode::Read);
	checkHeaders(inputFile.header(), {maskFile});
	Agilent::ProcPar pp; ReadPP(inputFile, pp);
	double TE1, ESP;
	if (pp) {
		TE1 = pp.realValue("te");
		ESP = pp.realValue("te2");
	} else {
		if (prompt) cout << "Enter first echo-time: " << flush;
		QUIT::Read<double>::FromLine(cin, TE1);
		if (prompt) cout << "Enter echo spacing: " << flush;
		QUIT::Read<double>::FromLine(cin, ESP);
	}
	if (NE < 1) {
		NE = inputFile.dim(4);
	} else {
		// Check that NE makes sense
		if ((inputFile.dim(4) % NE) != 0) {
			throw(runtime_error("Number of volumes is not a multiple of NE."));
		}
	}
	int NVols = inputFile.dim(4) / NE;

	// Set up echo times array
	MatrixXd X(NE, 2);
	X(0, 0) = TE1;
	for (int i = 1; i < X.rows(); i++) {
		X(i, 0) = X(i-1, 0) + ESP;
	}
	X.col(1).setOnes();

	if (verbose) {
		cout << "Ouput prefix will be: " << outPrefix << endl;
		cout << "Echo times are: " << X.col(0).transpose() << endl;
		cout << "Clamp: " << clamp_lo << " " << clamp_hi << endl;
		cout << "Thresh: " << thresh << endl;
	}
	cout << "Reading input data..." << flush;
	MultiArray<complex<float>, 4> inputVols(inputFile.dims().head(4));
	inputFile.readVolumes(inputVols.begin(), inputVols.end());
	cout << "done." << endl;
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	const auto dims = inputFile.matrix();
	MultiArray<float, 4> T2starVol(dims, NVols), PDVol(dims, NVols), ResVol(dims, NVols);
	MultiArray<float, 4> ResidsVols;
	if (all_residuals) {
		ResidsVols = MultiArray<float, 4>(dims, inputFile.dim(4));
	}
	for (size_t k = 0; k < dims[2]; k++) {
		clock_t loopStart;
		if (verbose) cout << "Starting slice " << k << "..." << flush;
		loopStart = clock();
		atomic<int> voxCount{0};
		function<void (const size_t, const size_t)> process = [&] (const size_t i, const size_t j) {
			if (!maskFile || (maskVol[{i,j,k}])) {
				voxCount++;
				for (size_t outVol = 0; outVol < NVols; outVol++) {
					double T2star = 0, PD = 0;
					ArrayXd signal = inputVols.slice<1>({i,j,k,outVol*NE},{0,0,0,NE}).asArray().abs().cast<double>();
					switch (algo) {
					case Algo::LogLin: {
						VectorXd Y = signal.log();
						VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
						T2star = -1 / b[0];
						PD = exp(b[1]);
					} break;
					case Algo::ARLO: {
						double si2sum = 0, sidisum = 0;
						for (int i = 0; i < NE - 2; i++) {
							double si = (ESP / 3) * (signal(i) + 4*signal(i+1) + signal(i+2));
							double di = signal(i) - signal(i+2);
							si2sum += si*si;
							sidisum = si*di;
						}
						T2star = (si2sum + (ESP/3)*sidisum) / ((ESP/3)*si2sum + sidisum);
						PD = (signal / (-X.col(0).array() / T2star).exp()).mean();
					} break;
					}

					if (PD < thresh) {
						PD = 0.;
						T2star = 0.;
					}
					T2star = clamp(T2star, clamp_lo, clamp_hi);
					ArrayXd theory = PD * (-X.col(0).array() / T2star).exp();
					ArrayXd resids = (signal - theory);
					if (all_residuals) {
						ResidsVols.slice<1>({i,j,k,outVol*NE},{0,0,0,NE}).asArray() = resids.cast<float>();
					}
					const MultiArray<float, 4>::Index idx{i,j,k,outVol};
					T2starVol[idx]  = static_cast<float>(T2star);
					PDVol[idx]  = static_cast<float>(PD);
					ResVol[idx] = static_cast<float>(sqrt(resids.square().sum() / resids.rows()) / PD);
				}
			}
		};

		threads.for_loop2(process, dims[0], dims[1]);

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
	outPrefix = outPrefix + "ME_";
	Nifti::Header outHdr = inputFile.header();
	outHdr.description = version;
	outHdr.setDim(4, NVols);
	outHdr.setDatatype(Nifti::DataType::FLOAT32);
	outHdr.intent = Nifti::Intent::Estimate;
	outHdr.intent_name = "T2* (seconds)";
	Nifti::File outFile(outHdr, outPrefix + "T2star" + OutExt());
	outFile.writeVolumes(T2starVol.begin(), T2starVol.end());
	outFile.close();
	outHdr.intent_name = "PD (au)";
	outFile.setHeader(outHdr);
	outFile.open(outPrefix + "PD" + OutExt(), Nifti::Mode::Write);
	outFile.writeVolumes(PDVol.begin(), PDVol.end());
	outFile.close();
	outHdr.intent_name = "Fractional Residual";
	outFile.setHeader(outHdr);
	outFile.open(outPrefix + "residual" + OutExt(), Nifti::Mode::Write);
	outFile.writeVolumes(ResVol.begin(), ResVol.end());
	outFile.close();
	if (all_residuals) {
		outHdr.intent_name = "Residuals";
		outHdr.setDim(4, inputFile.dim(4));
		outFile.setHeader(outHdr);
		outFile.open(outPrefix + "residuals" + OutExt(), Nifti::Mode::Write);
		outFile.writeVolumes(ResidsVols.begin(), ResidsVols.end(), 0, inputFile.dim(4));
		outFile.close();
	}
	cout << "All done." << endl;
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

