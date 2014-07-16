/*
 *  despot2_main.cpp
 *
 *  Created by Tobias Wood on 23/01/2012.
 *  Copyright (c) 2012-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <getopt.h>
#include <iostream>
#include <atomic>
#include <Eigen/Dense>
#include "unsupported/Eigen/NonLinearOptimization"
#include "unsupported/Eigen/NumericalDiff"

#include "Nifti/Nifti.h"
#include "QUIT/QUIT.h"
#include "DESPOT.h"
#include "DESPOT_Functors.h"

using namespace std;
using namespace Eigen;
using namespace QUIT;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: despot2 [options] T1_map ssfp_file\n\
\n\
Options:\n\
	--help, -h        : Print this message.\n\
	--mask, -m file   : Mask input with specified file.\n\
	--out, -o path    : Add a prefix to the output filenames.\n\
	--B1 file         : B1 Map file.\n\
	--elliptical, -e  : Input is band-free elliptical data.\n\
	--verbose, -v     : Print slice processing times.\n\
	--no-prompt, -n   : Suppress input prompts.\n\
	--algo, -a l      : LLS algorithm (default)\n\
	           w      : WLLS algorithm\n\
	           n      : NLLS (Levenberg-Marquardt)\n\
	--its, -i N      : Max iterations for WLLS (default 4)\n"
};

enum class Algos { LLS, WLLS, NLLS };
static int verbose = false, prompt = true, elliptical = false;
static size_t nIterations = 4;
static Algos algo;
static string outPrefix;
static struct option long_options[] =
{
	{"B1", required_argument, 0, '1'},
	{"elliptical", no_argument, 0, 'e'},
	{"help", no_argument, 0, 'h'},
	{"mask", required_argument, 0, 'm'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"algo", required_argument, 0, 'a'},
	{"its", required_argument, 0, 'i'},
	{0, 0, 0, 0}
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	try { // To fix uncaught exceptions on Mac

	Nifti::File maskFile, B0File, B1File;
	MultiArray<double, 3> maskVol, B1Vol;
	string procPath;

	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hm:o:b:vna:i:e", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'o':
				outPrefix = optarg;
				if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskVol.resize(maskFile.dims().head(3));
				maskFile.readVolumes(maskVol.begin(), maskVol.end(), 0, 1);
				break;
			case 'b':
				cout << "Reading B1 file: " << optarg << endl;
				B1File.open(optarg, Nifti::Mode::Read);
				B1Vol.resize(B1File.dims().head(3));
				B1File.readVolumes(B1Vol.begin(), B1Vol.end(), 0, 1);
				break;
			case 'a':
				switch (*optarg) {
					case 'l': algo = Algos::LLS;  if (verbose) cout << "LLS algorithm selected." << endl; break;
					case 'w': algo = Algos::WLLS; if (verbose) cout << "WLLS algorithm selected." << endl; break;
					case 'n': algo = Algos::NLLS; if (verbose) cout << "NLLS algorithm selected." << endl; break;
					default:
						cout << "Unknown algorithm type " << optarg << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 'i':
				nIterations = atoi(optarg);
				break;
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'e': elliptical = true; break;
			case 0:
				// Just a flag
				break;
			case '?': // getopt will print an error message
			case 'h':
				cout << usage << endl;				
				exit(EXIT_SUCCESS);
		}
	}
	if (verbose) cout << version << credit_shared << endl;
	Eigen::initParallel();
	if ((argc - optind) != 2) {
		cout << "Wrong number of arguments. Need a T1 map and 1 SSFP file." << endl;
		cout << usage << endl;
		exit(EXIT_FAILURE);
	}
	if (verbose) cout << "Reading T1 Map from: " << argv[optind] << endl;
	Nifti::File inFile(argv[optind++]);
	if ((maskFile.isOpen() && !inFile.header().matchesSpace(maskFile.header())) ||
		(B1File.isOpen() && !inFile.header().matchesSpace(B1File.header()))){
		cerr << "Dimensions/transforms do not match in input files." << endl;
		exit(EXIT_FAILURE);
	}
	MultiArray<double, 3> T1Vol(inFile.header().fulldims().head(3));
	inFile.readVolumes(T1Vol.begin(), T1Vol.end(), 0, 1);
	inFile.close();
	Nifti::Header outHdr = inFile.header(); // Save the header data to write out files
	//**************************************************************************
	// Gather SSFP Data
	//**************************************************************************
	Model ssfpMdl(Signal::Components::One, Model::Scaling::None);
	if (verbose) cout << "Reading SSFP data from: " << argv[optind] << endl;
	inFile.open(argv[optind], Nifti::Mode::Read);
	if (!inFile.header().matchesSpace(outHdr)) {
		cerr << "Dimensions/transforms do not match in input files." << endl;
		exit(EXIT_FAILURE);
	}
	MultiArray<complex<double>, 4> ssfpVols(inFile.header().fulldims().head(4));
	inFile.readVolumes(ssfpVols.begin(), ssfpVols.end());
	Agilent::ProcPar pp; ReadPP(inFile, pp);
	if (elliptical) {
		ssfpMdl.addSignal(SignalType::SSFP_Ellipse, prompt, pp);
	} else {
		ssfpMdl.addSignal(SignalType::SSFP, prompt, pp);
	}
	if (verbose) {
		cout << ssfpMdl;
		cout << "Ouput prefix will be: " << outPrefix << endl;
	}
	double TR = ssfpMdl.m_signals.at(0)->m_TR;
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	const auto dims = ssfpVols.dims().head(3);
	MultiArray<float, 3> T2Vol(dims), PDVol(dims), offResVol(dims), SoSVol(dims);
	time_t startTime;
	if (verbose)
		startTime = printStartTime();
	clock_t startClock = clock();
	int voxCount = 0;
	ThreadPool pool;
	for (size_t k = 0; k < dims(2); k++) {
		if (verbose)
			cout << "Starting slice " << k << "..." << flush;
		clock_t loopStart = clock();
		atomic<int> sliceCount{0};
		function<void (const int&)> process = [&] (const int &j) {
			for (size_t i = 0; i < dims(0); i++) {
				if (!maskFile.isOpen() || (maskVol[{i,j,k}])) {
					sliceCount++;
					double B1, T1, T2, E1, E2, PD, offRes, SoS;
					B1 = B1File.isOpen() ? B1Vol[{i,j,k}] : 1.;
					T1 = T1Vol[{i,j,k}];
					E1 = exp(-TR / T1);
					ArrayXd localAngles(ssfpMdl.m_signals.at(0)->B1flip(B1));
					auto data = ssfpVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray().cast<complex<double>>();
					auto s = data.abs();
					auto p = data.imag().binaryExpr(data.real(), ptr_fun<double,double,double>(atan2));
					offRes = p.mean() / (TR * M_PI);
					VectorXd Y = s / localAngles.sin();
					MatrixXd X(Y.rows(), 2);
					X.col(0) = s / localAngles.tan();
					X.col(1).setOnes();
					VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
					if (elliptical) {
						T2 = 2. * TR / log((b[0]*E1 - 1.) / (b[0] - E1));
						E2 = exp(-TR / T2);
						PD = b[1] * (1. - E1*E2*E2) / (sqrt(E2) * (1. - E1));
					} else {
						T2 = TR / log((b[0]*E1 - 1.)/(b[0] - E1));
						E2 = exp(-TR / T2);
						PD = b[1] * (1. - E1*E2) / (1. - E1);
					}
					if (algo == Algos::WLLS) {
						VectorXd W(ssfpMdl.size());
						for (size_t n = 0; n < nIterations; n++) {
							if (elliptical) {
								W = ((1. - E1*E2) * localAngles.sin() / (1. - E1*E2*E2 - (E1 - E2*E2)*localAngles.cos())).square();
							} else {
								W = ((1. - E1*E2) * localAngles.sin() / (1. - E1*E2 - (E1 - E2)*localAngles.cos())).square();
							}
							b = (X.transpose() * W.asDiagonal() * X).partialPivLu().solve(X.transpose() * W.asDiagonal() * Y);
							if (elliptical) {
								T2 = 2. * TR / log((b[0]*E1 - 1.) / (b[0] - E1));
								E2 = exp(-TR / T2);
								PD = b[1] * (1. - E1*E2*E2) / (sqrt(E2) * (1. - E1));
							} else {
								T2 = TR / log((b[0]*E1 - 1.)/(b[0] - E1));
								E2 = exp(-TR / T2);
								PD = b[1] * (1. - E1*E2) / (1. - E1);
							}
						}
					} else if (algo == Algos::NLLS) {
						DESPOTFunctor f(ssfpMdl, s.cast<complex<double>>(), B1, false, false);
						NumericalDiff<DESPOTFunctor> nDiff(f);
						LevenbergMarquardt<NumericalDiff<DESPOTFunctor>> lm(nDiff);
						lm.parameters.maxfev = nIterations;
						VectorXd p(4);
						p << PD, T1, T2, offRes;
						lm.lmder1(p);
						PD = p(0); T1 = p(1); T2 = p(2); offRes = p(3);
					}
					ArrayXcd theory = ssfpMdl.signal(Vector4d(PD, T1, T2, offRes), B1);
					SoS = (data - theory).abs2().sum();
					T2Vol[{i,j,k}]  = static_cast<float>(T2);
					PDVol[{i,j,k}]  = static_cast<float>(PD);
					offResVol[{i,j,k}] = static_cast<float>(offRes);
					SoSVol[{i,j,k}] = static_cast<float>(SoS);
				}
			}
		};
		pool.for_loop(process, dims(1));
		if (verbose) printLoopTime(loopStart, sliceCount);
		voxCount += sliceCount;
	}
	if (verbose) {
		printElapsedTime(startTime);
		cout << "Writing results." << endl;
	}
	printElapsedClock(startClock, voxCount);
	outHdr.description = version;
	outHdr.setDim(4, 1);
	outHdr.setDatatype(Nifti::DataType::FLOAT32);
	outHdr.intent = Nifti::Intent::Estimate;
	outHdr.intent_name = "T2 (s)";
	Nifti::File outFile(outHdr, outPrefix + "D2_T2" + OutExt());
	outFile.writeVolumes(T2Vol.begin(), T2Vol.end());
	outFile.close();
	outHdr.intent_name = "PD (au)";
	outFile.open(outPrefix + "D2_PD" + OutExt(), Nifti::Mode::Write);
	outFile.writeVolumes(PDVol.begin(), PDVol.end());
	outFile.close();
	outHdr.intent_name = "Sum of Squared Residuals";
	outFile.open(outPrefix + "D2_SoS" + OutExt(), Nifti::Mode::Write);
	outFile.writeVolumes(SoSVol.begin(), SoSVol.end());
	outFile.close();
	outHdr.intent_name = "Off-resonance (Hz)";
	outFile.open(outPrefix + "D2_f0" + OutExt(), Nifti::Mode::Write);
	outFile.writeVolumes(offResVol.begin(), offResVol.end());
	outFile.close();

	if (verbose) cout << "All done." << endl;
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	exit(EXIT_SUCCESS);
}
