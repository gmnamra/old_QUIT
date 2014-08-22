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

#include "Nifti/Nifti.h"
#include "QUIT/QUIT.h"
#include "Model.h"
#include "DESPOT_Functors.h"
#include "unsupported/Eigen/NonLinearOptimization"
#include "unsupported/Eigen/NumericalDiff"

using namespace std;
using namespace Eigen;
using namespace QUIT;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: despot1 [options] spgr_input \n\
\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--mask, -m file   : Mask input with specified file\n\
	--B1, -b file     : B1 Map file (ratio)\n\
	--algo, -a l      : LLS algorithm (default)\n\
	           w      : WLLS algorithm\n\
	           n      : NLLS (Levenberg-Marquardt)\n\
	--its, -i N      : Max iterations for WLLS (default 4)\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};

enum class Algos { LLS, WLLS, NLLS };
static bool verbose = false;
static size_t nIterations = 4;
static string outPrefix;
static Algos algo;
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"B1", required_argument, 0, 'b'},
	{"algo", required_argument, 0, 'a'},
	{"its", required_argument, 0, 'i'},
	{"threads", required_argument, 0, 'T'},
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
	cout << version << endl << credit_shared << endl;
	Eigen::initParallel();

	try { // To fix uncaught exceptions on Mac

	Nifti::File spgrFile, B1File, maskFile;
	MultiArray<float, 3> B1Vol;
	MultiArray<int8_t, 3> maskVol;
	ThreadPool threads;

	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvm:o:b:a:i:T:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
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
			case 'b':
				cout << "Reading B1 file: " << optarg << endl;
				B1File.open(optarg, Nifti::Mode::Read);
				B1Vol.resize(B1File.matrix());
				B1File.readVolumes(B1Vol.begin(), B1Vol.end());
				break;
			case 'a':
				switch (*optarg) {
					case 'l': algo = Algos::LLS;  cout << "LLS algorithm selected." << endl; break;
					case 'w': algo = Algos::WLLS; cout << "WLLS algorithm selected." << endl; break;
					case 'n': algo = Algos::NLLS; cout << "NLLS algorithm selected." << endl; break;
					default:
						cout << "Unknown algorithm type " << optarg << endl;
						return EXIT_FAILURE;
						break;
				} break;
			case 'i':
				nIterations = atoi(optarg);
				break;
			case 'T':
				threads.resize(atoi(optarg));
				break;
			case 'h':
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != 1) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		return EXIT_FAILURE;
	}
	
	//**************************************************************************
	#pragma mark Gather SPGR data
	//**************************************************************************
	cout << "Opening SPGR file: " << argv[optind] << endl;
	spgrFile.open(argv[optind], Nifti::Mode::Read);
	checkHeaders(spgrFile.header(), {maskFile, B1File});
	Agilent::ProcPar pp; ReadPP(spgrFile, pp);
	SPGRSimple spgrSequence(true, pp);
	if (verbose) {
		cout << spgrSequence;
		cout << "Ouput prefix will be: " << outPrefix << endl;
	}
	if (spgrSequence.size() != spgrFile.header().dim(4)) {
		throw(std::runtime_error("Specified number of flip-angles does not match number of volumes in file: " + spgrFile.imagePath()));
	}
	double TR = spgrSequence.m_TR;
	cout << "Reading SPGR data..." << flush;
	MultiArray<complex<float>, 4> spgrVols(spgrFile.dims().head(4));
	spgrFile.readVolumes(spgrVols.begin(), spgrVols.end());
	cout << "done." << endl;
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	const auto dims = spgrFile.matrix();
	MultiArray<float, 3> T1Vol(dims), PDVol(dims), SoSVol(dims);
	for (size_t k = 0; k < spgrFile.dim(3); k++) {
		clock_t loopStart;
		if (verbose)
			cout << "Starting slice " << k << "..." << flush;
		loopStart = clock();
		atomic<int> voxCount{0};
		function<void (const size_t)> process = [&] (const size_t j) {
			for (size_t i = 0; i < spgrFile.dim(1); i++) {
				if (!maskFile || (maskVol[{i,j,k}])) {
					voxCount++;
					double B1 = B1File ? B1Vol[{i,j,k}] : 1.;
					ArrayXd localAngles(spgrSequence.B1flip(B1));
					double T1, PD, SoS;
					ArrayXd signal = spgrVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray().abs().cast<double>();
					VectorXd Y = signal / localAngles.sin();
					MatrixXd X(Y.rows(), 2);
					X.col(0) = signal / localAngles.tan();
					X.col(1).setOnes();
					VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
					T1 = -TR / log(b[0]);
					PD = b[1] / (1. - b[0]);
					if (algo == Algos::WLLS) {
						VectorXd W(spgrSequence.size());
						for (size_t n = 0; n < nIterations; n++) {
							W = (localAngles.sin() / (1. - (exp(-TR/T1)*localAngles.cos()))).square();
							b = (X.transpose() * W.asDiagonal() * X).partialPivLu().solve(X.transpose() * W.asDiagonal() * Y);
							T1 = -TR / log(b[0]);
							PD = b[1] / (1. - b[0]);
						}
					} else if (algo == Algos::NLLS) {
						DESPOTFunctor f(spgrSequence, Pools::One, signal.cast<complex<double>>(), B1, false, false);
						NumericalDiff<DESPOTFunctor> nDiff(f);
						LevenbergMarquardt<NumericalDiff<DESPOTFunctor>> lm(nDiff);
						lm.parameters.maxfev = nIterations;
						VectorXd p(4);
						p << PD, T1, 0., 0.; // Don't need T2 of f0 for this (yet)
						lm.lmder1(p);
						PD = p(0); T1 = p(1);
					}
					ArrayXd theory = spgrSequence.signal(Pools::One, Vector4d(PD, T1, 0., 0.), B1).abs();
					SoS = (signal.abs() - theory.abs()).square().sum();
					T1Vol[{i,j,k}]  = static_cast<float>(T1);
					PDVol[{i,j,k}]  = static_cast<float>(PD);
					SoSVol[{i,j,k}] = static_cast<float>(SoS);
				}
			}
		};
			
		threads.for_loop(process, spgrFile.dim(2));
		
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
	Nifti::Header outHdr = spgrFile.header();
	outHdr.description = version;
	outHdr.setDim(4, 1);
	outHdr.setDatatype(Nifti::DataType::FLOAT32);
	outHdr.intent = Nifti::Intent::Estimate;
	outHdr.intent_name = "T1 (seconds)";
	Nifti::File outFile(outHdr, outPrefix + "D1_T1" + OutExt());
	outFile.writeVolumes(T1Vol.begin(), T1Vol.end());
	outFile.close();
	outHdr.intent_name = "PD (au)";
	outFile.open(outPrefix + "D1_PD" + OutExt(), Nifti::Mode::Write);
	outFile.writeVolumes(PDVol.begin(), PDVol.end());
	outFile.close();
	outHdr.intent_name = "Sum of Squared Residuals";
	outFile.open(outPrefix + "D1_SoS" + OutExt(), Nifti::Mode::Write);
	outFile.writeVolumes(SoSVol.begin(), SoSVol.end());
	outFile.close();

	cout << "All done." << endl;
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
