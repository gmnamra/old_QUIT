/*
 *  mcDESPOT_main.cpp
 *
 *  Created by Tobias Wood on 2013/08/12.
 *  Copyright (c) 2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <time.h>
#include <getopt.h>
#include <signal.h>
#include <iostream>
#include <atomic>
#include <Eigen/Dense>

#include "Nifti/Nifti.h"
#include "QUIT/QUIT.h"
#include "DESPOT.h"
#include "DESPOT_Functors.h"
#include "RegionContraction.h"

using namespace std;
using namespace Eigen;
using namespace QUIT;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: despot2-fm [options] T1_map ssfp_files\n\
\
Options:\n\
	--help, -h       : Print this message\n\
	--verbose, -v    : Print slice processing times\n\
	--no-prompt, -n  : Suppress input prompts\n\
	--mask, -m file  : Mask input with specified file\n\
	--out, -o path   : Add a prefix to the output filenames\n\
	--f0, -f SYM     : Fit symmetric f0 map (default)\n\
	         ASYM    : Fit asymmetric f0 map\n\
	         file    : Use f0 Map file (in Hertz)\n\
	--B1, -b file    : B1 Map file (ratio)\n\
	--start, -s N    : Start processing from slice N\n\
	--stop, -p  N    : Stop processing at slice N\n\
	--scale, -S 0    : Normalise signals to mean (default)\n\
	            1    : Fit a scaling factor/proton density\n\
	--tesla, -t 3    : Use boundaries suitable for 3T (default)\n\
	            7    : Boundaries suitable for 7T\n\
	            u    : User specified boundaries from stdin\n\
	--model, -M s    : Use simple model (default)\n\
	            f    : Use finite pulse length correction\n\
	--complex, -x    : Fit to complex data\n\
	--contract, -c n : Read contraction settings from stdin (Will prompt)\n"
};

static auto scale = Model::Scaling::NormToMean;
static auto tesla = Model::FieldStrength::Three;
static auto f0fit = OffRes::FitSym;
static size_t start_slice = 0, stop_slice = numeric_limits<size_t>::max();
static int verbose = false, prompt = true, writeResiduals = false,
           fitFinite = false, fitComplex = false,
           samples = 2000, retain = 20, contract = 10,
           voxI = 0, voxJ = 0;
static double expand = 0.;
static string outPrefix;
static struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"f0", required_argument, 0, 'f'},
	{"B1", required_argument, 0, 'b'},
	{"start", required_argument, 0, 's'},
	{"stop", required_argument, 0, 'p'},
	{"scale", required_argument, 0, 'S'},
	{"tesla", required_argument, 0, 't'},
	{"model", no_argument, 0, 'M'},
	{"complex", no_argument, 0, 'x'},
	{"contract", no_argument, 0, 'c'},
	{"resid", no_argument, 0, 'r'},
	{0, 0, 0, 0}
};

//******************************************************************************
#pragma mark SIGTERM interrupt handler and Threads
//******************************************************************************
ThreadPool threads;
bool interrupt_received = false;
void int_handler(int sig);
void int_handler(int) {
	cout << endl << "Stopping processing early." << endl;
	threads.stop();
	interrupt_received = true;
}

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	try { // To fix uncaught exceptions on Mac
	
	Nifti::Nifti1 maskFile, f0File, B1File;
	MultiArray<int8_t, 3> maskVol;
	MultiArray<float, 3> f0Vol, B1Vol;
	string procPath;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvnm:o:f:b:s:p:S:t:M:xcri:j:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'm':
				if (verbose) cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskVol.resize(maskFile.dims());
				maskFile.readVolumes(maskVol.begin(), maskVol.end(), 0, 1);
				break;
			case 'o':
				outPrefix = optarg;
				if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'f':
				if (string(optarg) == "SYM") {
					f0fit = OffRes::FitSym;
				} else if (string(optarg) == "ASYM") {
					f0fit = OffRes::Fit;
				} else {
					if (verbose) cout << "Reading f0 file: " << optarg << endl;
					f0File.open(optarg, Nifti::Mode::Read);
					f0Vol.resize(f0File.dims());
					f0File.readVolumes(f0Vol.begin(), f0Vol.end(), 0, 1);
					f0fit = OffRes::Map;
				}
				break;
			case 'b':
				if (verbose) cout << "Reading B1 file: " << optarg << endl;
				B1File.open(optarg, Nifti::Mode::Read);
				B1Vol.resize(B1File.dims());
				B1File.readVolumes(B1Vol.begin(), B1Vol.end(), 0, 1);
				break;
			case 's': start_slice = atoi(optarg); break;
			case 'p': stop_slice = atoi(optarg); break;
			case 'S':
				switch (atoi(optarg)) {
					case 0 : scale = Model::Scaling::NormToMean; break;
					case 1 : scale = Model::Scaling::None; break;
					default:
						cout << "Invalid scaling mode: " + to_string(atoi(optarg)) << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 't':
				switch (*optarg) {
					case '3': tesla = Model::FieldStrength::Three; break;
					case '7': tesla = Model::FieldStrength::Seven; break;
					case 'u': tesla = Model::FieldStrength::User; break;
					default:
						cout << "Unknown boundaries type " << optarg << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 'M':
				switch (*optarg) {
					case 's': fitFinite = false; cout << "Simple model selected." << endl; break;
					case 'f': fitFinite = true; cout << "Finite pulse correction selected." << endl; break;
					default:
						cout << "Unknown model type " << *optarg << endl;
						exit(EXIT_FAILURE);
						break;
				}
				break;
			case 'x':
				fitComplex = true;
				break;
			case 'c':
				cout << "Enter max number of contractions: " << flush; cin >> contract;
				cout << "Enter number of samples per contraction: " << flush; cin >> samples;
				cout << "Enter number of samples to retain: " << flush; cin >> retain;
				cout << "Enter fraction to expand region by: " << flush; cin >> expand;
				break;
			case 'r': writeResiduals = true; break;
			case 'i': voxI = atoi(optarg); break;
			case 'j': voxJ = atoi(optarg); break;
			case '?': // getopt will print an error message
			case 'h':
			default:
				cout << usage << endl;
				exit(EXIT_SUCCESS);
				break;
		}
	}
	if (verbose) cout << version << endl << credit_me << endl;
	Eigen::initParallel();
	if ((argc - optind) < 2) {
		cout << "Wrong number of arguments. Need at least a T1 map and 1 SSFP file." << endl;
		exit(EXIT_FAILURE);
	}
	if (verbose) cout << "Reading T1 Map from: " << argv[optind] << endl;
	Nifti::Nifti1 inFile(argv[optind++], Nifti::Mode::Read);
	MultiArray<float, 3> T1Vol{inFile.dims()};
	inFile.readVolumes(T1Vol.begin(), T1Vol.end(), 0, 1);
	inFile.close();
	if ((maskFile.isOpen() && !inFile.matchesSpace(maskFile)) ||
	    (f0File.isOpen() && !inFile.matchesSpace(f0File)) ||
		(B1File.isOpen() && !inFile.matchesSpace(B1File))){
		cerr << "Dimensions/transforms do not match in input files." << endl;
		exit(EXIT_FAILURE);
	}
	Nifti::Nifti1 templateFile(inFile, 1); // Save header data to write out results
	//**************************************************************************
	// Gather SSFP Data
	//**************************************************************************
	size_t nFiles = argc - optind;
	vector<MultiArray<complex<float>, 4>> ssfpData(nFiles);
	Model model(Signal::Components::One, scale);
	VectorXd inFlip;
	for (size_t p = 0; p < nFiles; p++) {
		if (verbose) cout << "Reading SSFP header from " << argv[optind] << endl;
		inFile.open(argv[optind], Nifti::Mode::Read);
		if (p == 0)
			templateFile = Nifti::Nifti1(inFile, 1);
		if (!inFile.matchesSpace(templateFile)) {
			cerr << "Input file dimensions and/or transforms do not match." << endl;
			exit(EXIT_FAILURE);
		}
		Agilent::ProcPar pp; ReadPP(inFile, pp);
		if (fitFinite) {
			model.addSignal(SignalType::SSFP_Finite, prompt, pp);
		} else {
			model.addSignal(SignalType::SSFP, prompt, pp);
		}
		if (verbose) cout << "Reading data." << endl;
		ssfpData.at(p).resize(inFile.dims());
		inFile.readVolumes(ssfpData.at(p).begin(), ssfpData.at(p).end());
		inFile.close();
		optind++;
	}
	if (optind != argc) {
		cerr << "Unprocessed arguments supplied.\n" << usage;
		exit(EXIT_FAILURE);
	}
	
	ArrayXd thresh(model.nParameters()); thresh.setConstant(0.05);
	ArrayXXd bounds = model.bounds(tesla);
	if (tesla == Model::FieldStrength::User) {
		cout << "Enter parameter pairs (low then high)" << endl;
		for (size_t i = 0; i < model.nParameters() - 1; i++) {
			cout << model.names()[i] << ": " << flush;
			cin >> bounds(i, 0) >> bounds(i, 1);
		}
	}
	if (f0fit == OffRes::FitSym) {
		bounds(model.nParameters() - 1, 0) = 0.;
	}
	ArrayXd weights(model.size()); weights.setOnes();
	
	if (verbose) {
		cout << model;
		cout << "Bounds:" << endl <<  bounds.transpose() << endl;
	}
	//**************************************************************************
	// Set up results data
	//**************************************************************************
	size_t nParams = 2;
	if (scale == Model::Scaling::None)
		nParams = 3;
	MultiArray<float, 4> paramsVols(templateFile.dims().head(3), nParams);
	MultiArray<float, 4> residualVols(templateFile.dims().head(3), model.size());
	MultiArray<float, 3> SoSVol(T1Vol.dims());
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	if (stop_slice > templateFile.dim(3))
		stop_slice = templateFile.dim(3);
	time_t startTime;
	if (verbose) startTime = printStartTime();
	clock_t startClock = clock();
	int voxCount = 0;
	signal(SIGINT, int_handler);	// If we've got here there's actually allocated data to save
	for (size_t k = start_slice; k < stop_slice; k++) {
		if (verbose) cout << "Starting slice " << k << "..." << flush;
		atomic<int> sliceCount{0};
		clock_t loopStart = clock();

		function<void (const size_t&)> processVox = [&] (const size_t &j) {
			for (size_t i = 0; i < T1Vol.dims()[0]; i++) {
				if (!maskFile.isOpen() || (maskVol[{i,j,k}] && T1Vol[{i,j,k}] > 0.)) {
					// -ve T1 is nonsensical, no point fitting
					sliceCount++;
					ArrayXcd signal = model.loadSignals(ssfpData, i, j, k);
					ArrayXXd localBounds = bounds;
					if (scale == Model::Scaling::None) {
						localBounds(0, 0) = 0.;
						localBounds(0, 1) = signal.abs().maxCoeff() * 25;
					}
					localBounds.row(1).setConstant(T1Vol[{i,j,k}]);
					if (f0fit == OffRes::Map) {
						localBounds.row(3).setConstant(f0Vol[{i,j,k}]);
					}
					double B1 = B1File.isOpen() ? B1Vol[{i,j,k}] : 1.;
					DESPOTFunctor func(model, signal, B1, fitComplex, false);
					RegionContraction<DESPOTFunctor> rc(func, localBounds, weights, thresh,
														samples, retain, contract, expand, (voxI > 0));
					ArrayXd params(model.nParameters()); params.setZero();
					rc.optimise(params, time(NULL) + i); // Add the voxel number to the time to get a decent random seed

					if (scale == Model::Scaling::None) {
						// Skip T1
						paramsVols[{i,j,k,0}] = params(0);
						paramsVols[{i,j,k,1}] = params(2);
						paramsVols[{i,j,k,2}] = params(3);
					} else {
						paramsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = params.tail(2).cast<float>(); // Skip PD & T1
					}
					SoSVol[{i,j,k}] = static_cast<float>(rc.SoS());
					if (writeResiduals) {
						residualVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = rc.residuals().cast<float>();
					}
				}
			}
		};
		if (voxI == 0)
			threads.for_loop(processVox, T1Vol.dims()[1]);
		else {
			processVox(voxI);
			exit(0);
		}
		
		if (verbose) printLoopTime(loopStart, sliceCount);
		voxCount += sliceCount;
		if (interrupt_received)
			break;
	}
    if (verbose) printElapsedTime(startTime);
    printElapsedClock(startClock, voxCount);
	
	outPrefix = outPrefix + "FM_";
	templateFile.description = version;
	if (scale == Model::Scaling::None) {
		templateFile.open(outPrefix + model.names().at(0) + ".nii.gz", Nifti::Mode::Write);
		auto p = paramsVols.slice<3>({0,0,0,0},{-1,-1,-1,0});
		templateFile.writeVolumes(p.begin(), p.end());
		templateFile.close();
		templateFile.open(outPrefix + model.names().at(2) + ".nii.gz", Nifti::Mode::Write);
		p = paramsVols.slice<3>({0,0,0,1},{-1,-1,-1,0});
		templateFile.writeVolumes(p.begin(), p.end());
		templateFile.close();
		templateFile.open(outPrefix + model.names().at(3) + ".nii.gz", Nifti::Mode::Write);
		p = paramsVols.slice<3>({0,0,0,2},{-1,-1,-1,0});
		templateFile.writeVolumes(p.begin(), p.end());
		templateFile.close();
	} else {
		templateFile.open(outPrefix + model.names().at(2) + ".nii.gz", Nifti::Mode::Write);
		auto p = paramsVols.slice<3>({0,0,0,0},{-1,-1,-1,0});
		templateFile.writeVolumes(p.begin(), p.end());
		templateFile.close();
		templateFile.open(outPrefix + model.names().at(3) + ".nii.gz", Nifti::Mode::Write);
		p = paramsVols.slice<3>({0,0,0,1},{-1,-1,-1,0});
		templateFile.writeVolumes(p.begin(), p.end());
		templateFile.close();
	}
	templateFile.open(outPrefix + "SoS.nii.gz", Nifti::Mode::Write);
	templateFile.writeVolumes(SoSVol.begin(), SoSVol.end());
	templateFile.close();
	if (writeResiduals) {
		templateFile.setDim(4, static_cast<int>(model.size()));
		templateFile.open(outPrefix + "residuals.nii.gz", Nifti::Mode::Write);
		templateFile.writeVolumes(residualVols.begin(), residualVols.end());
		templateFile.close();
	}
	
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
