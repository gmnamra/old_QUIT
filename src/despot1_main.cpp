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
#include "Nifti/Volume.h"
#include "ThreadPool.h"
#include "Model.h"
#include "DESPOT_Functors.h"
#include "unsupported/Eigen/NonLinearOptimization"
#include "unsupported/Eigen/NumericalDiff"

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
	--its, -i N      : Max iterations for WLLS (default 4)\n"
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
	Nifti spgrFile, B1File, maskFile;
	Volume<float> B1Vol;
	Volume<bool> maskVol;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvm:o:b:a:i:", long_options, &indexptr)) != -1) {
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
			case 'b':
				cout << "Reading B1 file: " << optarg << endl;
				B1File.open(optarg, Nifti::Mode::Read);
				B1Vol.readFrom(B1File);
				break;
			case 'a':
				switch (*optarg) {
					case 'l': algo = Algos::LLS;  cout << "LLS algorithm selected." << endl; break;
					case 'w': algo = Algos::WLLS; cout << "WLLS algorithm selected." << endl; break;
					case 'n': algo = Algos::NLLS; cout << "NLLS algorithm selected." << endl; break;
					default:
						cout << "Unknown algorithm type " << optarg << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 'n':
				nIterations = atoi(optarg);
				break;
			case 'h':
			case '?': // getopt will print an error message
				exit(EXIT_FAILURE);
		}
	}
	if ((argc - optind) != 1) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		exit(EXIT_FAILURE);
	}
	
	//**************************************************************************
	#pragma mark Gather SPGR data
	//**************************************************************************
	SimpleModel spgrMdl(Signal::Components::One, Model::Scaling::None);
	cout << "Opening SPGR file: " << argv[optind] << endl;
	spgrFile.open(argv[optind], Nifti::Mode::Read);
	if ((maskFile.isOpen() && !maskFile.matchesSpace(spgrFile)) ||
	    (B1File.isOpen() && !B1File.matchesSpace(spgrFile))) {
		cerr << "Mask or B1 dimensions/transform do not match SPGR file." << endl;
		exit(EXIT_FAILURE);
	}
	#ifdef AGILENT
	ProcPar pp;
	if (ReadPP(spgrFile, pp)) {
		spgrMdl.procparseSPGR(pp);
	} else
	#endif
	{
		spgrMdl.parseSPGR(spgrFile.dim(4), true);
	}
	if (verbose) {
		cout << spgrMdl;
		cout << "Ouput prefix will be: " << outPrefix << endl;
	}
	double TR = spgrMdl.m_signals.at(0)->m_TR;
	//**************************************************************************
	// Allocate memory for slices
	//**************************************************************************
	cout << "Reading SPGR data..." << flush;
	VolumeSeries<complex<float>> spgrVol(spgrFile);
	cout << "done." << endl;
	//**************************************************************************
	// Create results data storage
	//**************************************************************************
	Volume<float> T1Vol(spgrVol.dims().head(3)), PDVol(spgrVol.dims().head(3)),
	              SoSVol(spgrVol.dims().head(3));
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	ThreadPool pool;
	for (size_t k = 0; k < spgrFile.dim(3); k++) {
		clock_t loopStart;
		// Read in data
		if (verbose)
			cout << "Starting slice " << k << "..." << flush;
		loopStart = clock();
		atomic<int> voxCount{0};
		
		for (size_t j = 0; j < spgrFile.dim(2); j++) {
			function<void (const size_t)> processVox = [&] (const size_t i) {
				const typename Volume<float>::IndexArray vox{i, j, k};
				if (!maskFile.isOpen() || (maskVol[vox])) {
					voxCount++;
					//cout << spgrMdl.m_signals[0]->m_flip.transpose() << endl;
					double B1 = B1File.isOpen() ? B1Vol[vox] : 1.;
					ArrayXd localAngles(spgrMdl.m_signals.at(0)->B1flip(B1));
					//cout << localAngles.transpose() << endl;
					double T1, PD, SoS;
					ArrayXd signal = spgrVol.series(vox).abs().cast<double>();
					VectorXd Y = signal / localAngles.sin();
					MatrixXd X(Y.rows(), 2);
					X.col(0) = signal / localAngles.tan();
					X.col(1).setOnes();
					VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
					T1 = -TR / log(b[0]);
					PD = b[1] / (1. - b[0]);
					if (algo == Algos::WLLS) {
						VectorXd W(spgrMdl.size());
						for (size_t n = 0; n < nIterations; n++) {
							W = (localAngles.sin() / (1. - (exp(-TR/T1)*localAngles.cos()))).square();
							b = (X.transpose() * W.asDiagonal() * X).partialPivLu().solve(X.transpose() * W.asDiagonal() * Y);
							T1 = -TR / log(b[0]);
							PD = b[1] / (1. - b[0]);
						}
					} else if (algo == Algos::NLLS) {
						DESPOTFunctor f(make_shared<SimpleModel>(spgrMdl), signal.cast<complex<double>>(), B1, false, false);
						NumericalDiff<DESPOTFunctor> nDiff(f);
						LevenbergMarquardt<NumericalDiff<DESPOTFunctor>> lm(nDiff);
						lm.parameters.maxfev = nIterations;
						VectorXd p(4);
						p << PD, T1, 0., 0.; // Don't need T2 of f0 for this (yet)
						lm.lmder1(p);
						PD = p(0); T1 = p(1);
					}
					ArrayXd theory = spgrMdl.signal(Vector4d(PD, T1, 0., 0.), B1).abs();
					SoS = (signal - theory).square().sum();
					T1Vol[vox]  = static_cast<float>(T1);
					PDVol[vox]  = static_cast<float>(PD);
					SoSVol[vox] = static_cast<float>(SoS);
				}
			};
			
			pool.for_loop(processVox, spgrFile.dim(1));
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
	Nifti outFile(spgrFile, 1);
	outFile.description = version;
	outFile.open(outPrefix + "D1_T1.nii.gz", Nifti::Mode::Write);
	T1Vol.writeTo(outFile);
	outFile.close();
	outFile.open(outPrefix + "D1_PD.nii.gz", Nifti::Mode::Write);
	PDVol.writeTo(outFile);
	outFile.close();
	outFile.open(outPrefix + "D1_SoS.nii.gz", Nifti::Mode::Write);
	SoSVol.writeTo(outFile);
	outFile.close();

	cout << "All done." << endl;
	exit(EXIT_SUCCESS);
}
