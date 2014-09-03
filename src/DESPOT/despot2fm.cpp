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
	--help, -h        : Print this message\n\
	--verbose, -v     : Print slice processing times\n\
	--no-prompt, -n   : Suppress input prompts\n\
	--mask, -m file   : Mask input with specified file\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--f0, -f SYM      : Fit symmetric f0 map (default)\n\
	         ASYM     : Fit asymmetric f0 map\n\
	         file     : Use f0 Map file (in Hertz)\n\
	--B1, -b file     : B1 Map file (ratio)\n\
	--start, -s N     : Start processing from slice N\n\
	--stop, -p  N     : Stop processing at slice N\n\
	--scale, -S 0     : Normalise signals to mean (default)\n\
	            1     : Fit a scaling factor/proton density\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n\
	--sequences, -M s : Use simple sequences (default)\n\
	            f     : Use finite pulse length correction\n\
	--complex, -x     : Fit to complex data\n\
	--contract, -c n  : Read contraction settings from stdin (Will prompt)\n"
};

static auto scale = Scale::NormToMean;
static auto tesla = FieldStrength::Three;
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
	{"threads", required_argument, 0, 'T'},
	{"sequences", no_argument, 0, 'M'},
	{"complex", no_argument, 0, 'x'},
	{"contract", no_argument, 0, 'c'},
	{"resid", no_argument, 0, 'r'},
	{0, 0, 0, 0}
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	try { // To fix uncaught exceptions on Mac
	cout << version << endl << credit_me << endl;
	Eigen::initParallel();
	Nifti::File maskFile, f0File, B1File;
	MultiArray<int8_t, 3> maskVol;
	MultiArray<float, 3> f0Vol, B1Vol;
	string procPath;
	ThreadPool threads;
	//ThreadPool::EnableDebug = true;

	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvnm:o:f:b:s:p:S:T:M:xcri:j:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'm':
				if (verbose) cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskVol.resize(maskFile.matrix());
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
				B1Vol.resize(B1File.matrix());
				B1File.readVolumes(B1Vol.begin(), B1Vol.end(), 0, 1);
				break;
			case 's': start_slice = atoi(optarg); break;
			case 'p': stop_slice = atoi(optarg); break;
			case 'S':
				switch (atoi(optarg)) {
					case 0 : scale = Scale::NormToMean; break;
					case 1 : scale = Scale::None; break;
					default:
						cout << "Invalid scaling mode: " + to_string(atoi(optarg)) << endl;
						return EXIT_FAILURE;
						break;
				} break;
			case 'T':
				threads.resize(atoi(optarg));
				break;
			case 'M':
				switch (*optarg) {
					case 's': fitFinite = false; cout << "Simple sequences selected." << endl; break;
					case 'f': fitFinite = true; cout << "Finite pulse correction selected." << endl; break;
					default:
						cout << "Unknown sequences type " << *optarg << endl;
						return EXIT_FAILURE;
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
				return EXIT_SUCCESS;
				break;
		}
	}
	if ((argc - optind) < 2) {
		cout << "Wrong number of arguments. Need at least a T1 map and 1 SSFP file." << endl;
		return EXIT_FAILURE;
	}
	if (verbose) cout << "Reading T1 Map from: " << argv[optind] << endl;
	Nifti::File T1File(argv[optind++]);
	const auto dims = T1File.matrix();
	MultiArray<float, 3> T1Vol(dims);
	T1File.readVolumes(T1Vol.begin(), T1Vol.end(), 0, 1);
	checkHeaders(T1File.header(), {maskFile, f0File, B1File});
	//**************************************************************************
	// Gather SSFP Data
	//**************************************************************************
	size_t nFiles = argc - optind;
	vector<MultiArray<complex<float>, 4>> ssfpData(nFiles);
	Sequences sequences(scale);
	VectorXd inFlip;
	for (size_t p = 0; p < nFiles; p++) {
		if (verbose) cout << "Reading SSFP header from " << argv[optind] << endl;
		Nifti::File inFile(argv[optind]);
		checkHeaders(inFile.header(), {T1File});
		Agilent::ProcPar pp; ReadPP(inFile, pp);
		if (fitFinite) {
			sequences.addSequence(SequenceType::SSFP_Finite, prompt, pp);
		} else {
			sequences.addSequence(SequenceType::SSFP, prompt, pp);
		}
		if (verbose) cout << "Reading data." << endl;
		ssfpData.at(p).resize(inFile.dims().head(4));
		inFile.readVolumes(ssfpData.at(p).begin(), ssfpData.at(p).end());
		inFile.close();
		optind++;
	}
	if (optind != argc) {
		cerr << "Unprocessed arguments supplied.\n" << usage;
		return EXIT_FAILURE;
	}

	ArrayXd thresh(PoolInfo::nParameters(Pools::One)); thresh.setConstant(0.05);
	ArrayXd weights(sequences.size()); weights.setOnes();
	Array2d f0Bounds(-0.5/sequences.minTR(),0.5/sequences.minTR());
	if (f0fit == OffRes::FitSym) {
		f0Bounds(0) = 0.;
	}
	
	if (verbose) {
		cout << sequences;
	}
	//**************************************************************************
	// Set up results data
	//**************************************************************************
	size_t nParams = 2;
	if (scale == Scale::None)
		nParams = 3;
	MultiArray<float, 4> paramsVols(dims, nParams);
	MultiArray<float, 4> residualVols(dims, sequences.size());
	MultiArray<float, 3> SoSVol(dims);
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	if (stop_slice > dims[2])
		stop_slice = dims[2];
	time_t startTime;
	if (verbose) startTime = printStartTime();
	clock_t startClock = clock();
	int voxCount = 0;
	for (size_t k = start_slice; k < stop_slice; k++) {
		if (verbose) cout << "Starting slice " << k << "..." << flush;
		atomic<int> sliceCount{0};
		clock_t loopStart = clock();

		function<void (const size_t, const size_t)> processVox = [&] (const size_t i, const size_t j) {
			const MultiArray<float, 3>::Index idx{i,j,k};
			if (!maskFile || (maskVol[idx] && T1Vol[idx] > 0.)) {
				// -ve T1 is nonsensical, no point fitting
				sliceCount++;
				ArrayXcd signal = sequences.loadSignals(ssfpData, i, j, k);
				ArrayXXd bounds(PoolInfo::nParameters(Pools::One), 2);
				bounds.setZero();
				if (scale == Scale::None) {
					bounds(0, 0) = 0.;
					bounds(0, 1) = signal.abs().maxCoeff() * 25;
				} else {
					bounds.row(0).setConstant(1.);
				}
				bounds.row(1).setConstant(T1Vol[idx]);
				bounds(2,0) = 0.001;
				bounds(2,1) = T1Vol[idx];
				if (f0fit == OffRes::Map) {
					bounds.row(3).setConstant(f0Vol[idx]);
				} else {
					bounds.row(3) = f0Bounds;
				}
				double B1 = B1File ? B1Vol[{i,j,k}] : 1.;
				DESPOTFunctor func(sequences, Pools::One, signal, B1, fitComplex, false);
				RegionContraction<DESPOTFunctor> rc(func, bounds, weights, thresh,
													samples, retain, contract, expand, (voxI > 0));
				ArrayXd params(PoolInfo::nParameters(Pools::One)); params.setZero();
				rc.optimise(params, time(NULL) + i); // Add the voxel number to the time to get a decent random seed
				if (scale == Scale::None) {
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
		};
		if (voxI == 0) {
			threads.for_loop2(processVox, dims[0], dims[1]);
			if (threads.interrupted())
				break;
		} else {
			processVox(voxI, voxJ);
			voxCount = 1;
			break;
		}
		
		if (verbose) printLoopTime(loopStart, sliceCount);
		voxCount += sliceCount;
	}
	if (verbose) printElapsedTime(startTime);
	printElapsedClock(startClock, voxCount);
	if (voxI != 0)
		return EXIT_SUCCESS;

	outPrefix = outPrefix + "FM_";
	Nifti::Header hdr = T1File.header();
	hdr.setDim(4, 1);
	hdr.setDatatype(Nifti::DataType::FLOAT32);
	hdr.description = version;
	hdr.intent = Nifti::Intent::Estimate;
	if (scale == Scale::None) {
		hdr.intent_name = PoolInfo::Names(Pools::One).at(0);
		Nifti::File out(hdr, outPrefix + PoolInfo::Names(Pools::One).at(0) + OutExt());
		auto p = paramsVols.slice<3>({0,0,0,0},{-1,-1,-1,0});
		out.writeVolumes(p.begin(), p.end());
		out.close();
		hdr.intent_name = PoolInfo::Names(Pools::One).at(2);
		out.open(outPrefix + PoolInfo::Names(Pools::One).at(2) + OutExt(), Nifti::Mode::Write);
		p = paramsVols.slice<3>({0,0,0,1},{-1,-1,-1,0});
		out.writeVolumes(p.begin(), p.end());
		out.close();
		hdr.intent_name = PoolInfo::Names(Pools::One).at(3);
		out.open(outPrefix + PoolInfo::Names(Pools::One).at(3) + OutExt(), Nifti::Mode::Write);
		p = paramsVols.slice<3>({0,0,0,2},{-1,-1,-1,0});
		out.writeVolumes(p.begin(), p.end());
		out.close();
	} else {
		hdr.intent_name = PoolInfo::Names(Pools::One).at(2);
		Nifti::File out(hdr, outPrefix + PoolInfo::Names(Pools::One).at(2) + OutExt());
		auto p = paramsVols.slice<3>({0,0,0,0},{-1,-1,-1,0});
		out.writeVolumes(p.begin(), p.end());
		out.close();
		hdr.intent_name = PoolInfo::Names(Pools::One).at(3);
		out.open(outPrefix + PoolInfo::Names(Pools::One).at(3) + OutExt(), Nifti::Mode::Write);
		p = paramsVols.slice<3>({0,0,0,1},{-1,-1,-1,0});
		out.writeVolumes(p.begin(), p.end());
		out.close();
	}
	hdr.intent_name = "Sum of Squared Residuals";
	Nifti::File SoS(hdr, outPrefix + "SoS" + OutExt());
	SoS.writeVolumes(SoSVol.begin(), SoSVol.end());
	SoS.close();
	if (writeResiduals) {
		hdr.setDim(4, static_cast<int>(sequences.size()));
		Nifti::File res(hdr, outPrefix + "residuals" + OutExt());
		res.writeVolumes(residualVols.begin(), residualVols.end());
		res.close();
	}
	
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
