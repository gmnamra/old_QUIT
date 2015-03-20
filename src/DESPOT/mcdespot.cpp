/*
 *  mcdespot_main.cpp
 *
 *  Created by Tobias Wood on 14/02/2012.
 *  Copyright (c) 2012-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <atomic>
#include <getopt.h>
#include <time.h>
#include <fstream>
#include <Eigen/Dense>

#include "Nifti/Nifti.h"
#include "DESPOT_Functors.h"
#include "QUIT/QUIT.h"
#include "RegionContraction.h"

using namespace std;
using namespace Eigen;
using namespace QUIT;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: mcdespot [options]\n\
\n\
The program will prompt for input (unless --no-prompt specified)\n\
\n\
All times (TR) are in SECONDS. All angles are in degrees.\n\
\n\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--no-prompt, -n   : Don't print prompts for input\n\
	--mask, -m file   : Mask input with specified file\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--1, --2, --3     : Use 1, 2 or 3 component sequences (default 3)\n\
	--f0, -f SYM      : Fit symmetric f0 map (default)\n\
	         ASYM     : Fit asymmetric f0 map\n\
	         file     : Use f0 Map file (in Hertz)\n\
	--B1, -b file     : B1 Map file (ratio)\n\
	--start, -s n     : Only start processing at slice n.\n\
	--stop, -p n      : Finish at slice n-1\n\
	--scale, -S 0     : Normalise signals to mean (default)\n\
	            1     : Fit a scaling factor/proton density\n\
	--flip, -F        : Data order is phase, then flip-angle (default opposite)\n\
	--tesla, -t 3     : Boundaries suitable for 3T (default)\n\
	            7     : Boundaries suitable for 7T \n\
	            u     : User specified boundaries from stdin\n\
	--sequences, -M s : Use simple sequences (default)\n\
	            f     : Use Finite Pulse Length correction\n\
	--complex, -x     : Fit to complex data\n\
	--contract, -c n  : Read contraction settings from stdin (Will prompt)\n\
	--resids, -r      : Write out per flip-angle residuals\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};

static auto pools = Pools::Three;
static auto scale = Scale::NormToMean;
static auto tesla = FieldStrength::Three;
static auto f0fit = OffRes::FitSym;
static size_t start_slice = 0, stop_slice = numeric_limits<size_t>::max();
static int verbose = false, prompt = true, all_residuals = false,
           fitFinite = false, fitComplex = false, flipData = false,
           samples = 5000, retain = 50, contract = 10,
           voxI = 0, voxJ = 0;
static double expand = 0.;
static string outPrefix;
static const struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"f0", required_argument, 0, 'f'},
	{"B1", required_argument, 0, 'b'},
	{"start", required_argument, 0, 's'},
	{"stop", required_argument, 0, 'p'},
	{"scale", required_argument, 0, 'S'},
	{"flip", required_argument, 0, 'F'},
	{"tesla", required_argument, 0, 't'},
	{"sequences", no_argument, 0, 'M'},
	{"complex", no_argument, 0, 'x'},
	{"contract", no_argument, 0, 'c'},
	{"resids", no_argument, 0, 'r'},
	{"threads", required_argument, 0, 'T'},
	{"no-prompt", no_argument, 0, 'n'},
	{"1", no_argument, 0, '1'},
	{"2", no_argument, 0, '2'},
	{"3", no_argument, 0, '3'},
	{0, 0, 0, 0}
};
static const char* short_options = "hvm:o:f:b:s:p:S:t:FT:M:xcrn123i:j:";

//******************************************************************************
#pragma mark Read in all required files and data from cin
//******************************************************************************
Nifti::Header parseInput(Sequences &seq, vector<MultiArray<complex<float>, 4>> &signalVols);
Nifti::Header parseInput(Sequences &seq, vector<MultiArray<complex<float>, 4>> &signalVols)
{
	Nifti::Header hdr;
	string type, path;
	if (prompt) cout << "Specify next image type (SPGR/SSFP): " << flush;
	while (getline(cin, type) && (type != "END") && (type != "")) {
		if (type != "SPGR" && type != "SSFP") {
			throw(std::runtime_error("Unknown signal type: " + type));
		}
		if (prompt) cout << "Enter image path: " << flush;
		getline(cin, path);
		Nifti::File inFile(path);
		if (signalVols.size() == 0) {
			hdr = inFile.header(); // Save header info for later
		} else {
			checkHeaders(hdr, {inFile});
		}
		if (verbose) cout << "Opened: " << inFile.imagePath() << endl;
		Agilent::ProcPar pp; ReadPP(inFile, pp);
		if ((type == "SPGR") && !fitFinite) {
			seq.addSequence(SequenceType::SPGR, prompt, pp);
		} else if ((type == "SPGR" && fitFinite)) {
			seq.addSequence(SequenceType::SPGR_Finite, prompt, pp);
		} else if ((type == "SSFP" && !fitFinite)) {
			seq.addSequence(SequenceType::SSFP, prompt, pp);
		} else if ((type == "SSFP" && fitFinite)) {
			seq.addSequence(SequenceType::SSFP_Finite, prompt, pp);
		}
		if (seq.sequence(seq.count() - 1)->size() != inFile.dim(4)) {
			throw(std::runtime_error("Number of volumes in file " + inFile.imagePath() + " does not match input."));
		}
		MultiArray<complex<float>, 4> inData(inFile.dims().head(4));
		if (verbose) cout << "Reading data..." << flush;
		inFile.readVolumes(inData.begin(), inData.end());
		signalVols.push_back(inData);
		inFile.close();
		if (verbose) cout << "done." << endl;
		// Print message ready for next loop
		if (prompt) cout << "Specify next image type (SPGR/SSFP, END to finish input): " << flush;
	}
	return hdr;
}
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	try { // To fix uncaught exceptions on Mac
	cout << version << endl << credit_me << endl;
	Eigen::initParallel();
	Nifti::File maskFile, f0File, B1File;
	MultiArray<int8_t, 3> maskVol;
	MultiArray<float, 3> f0Vol, B1Vol;
	ThreadPool threads;

	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_options, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case '1': pools = Pools::One; break;
			case '2': pools = Pools::Two; break;
			case '3': pools = Pools::Three; break;
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
			case 'f':
				if (string(optarg) == "SYM") {
					f0fit = OffRes::FitSym;
				} else if (string(optarg) == "ASYM") {
					f0fit = OffRes::Fit;
				} else {
					cout << "Reading f0 file: " << optarg << endl;
					f0File.open(optarg, Nifti::Mode::Read);
					f0Vol.resize(f0File.matrix());
					f0File.readVolumes(f0Vol.begin(), f0Vol.end(), 0, 1);
					f0fit = OffRes::Map;
				}
				break;
			case 'b':
				cout << "Reading B1 file: " << optarg << endl;
				B1File.open(optarg, Nifti::Mode::Read);
				B1Vol.resize(B1File.matrix());
				B1File.readVolumes(B1Vol.begin(), B1Vol.end(), 0, 1);
				break;
			case 's': start_slice = atoi(optarg); break;
			case 'p': stop_slice = atoi(optarg); break;
			case 'S':
				switch (atoi(optarg)) {
					case 1 : scale = Scale::None; break;
					case 2 : scale = Scale::NormToMean; break;
					default:
						cout << "Invalid scaling mode: " + to_string(atoi(optarg)) << endl;
						return EXIT_FAILURE;
						break;
				} break;
			case 'F':
				flipData = true;
				break;
			case 'T':
				threads.resize(atoi(optarg));
				break;
			case 't':
				switch (*optarg) {
					case '3': tesla = FieldStrength::Three; break;
					case '7': tesla = FieldStrength::Seven; break;
					case 'u': tesla = FieldStrength::User; break;
					default:
						cout << "Unknown boundaries type " << *optarg << endl;
						return EXIT_FAILURE;
						break;
				} break;
			case 'M':
				switch (*optarg) {
					case 's': fitFinite = false; if (verbose) cout << "Simple sequences selected." << endl; break;
					case 'f': fitFinite = true;  if (verbose) cout << "Finite pulse correction selected." << endl; break;
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
				{ string dummy; getline(cin, dummy); } // Eat newlines
				break;
			case 'r': all_residuals = true; break;
			case 'i': voxI = atoi(optarg); break;
			case 'j': voxJ = atoi(optarg); break;
			case 'h':
			case '?': // getopt will print an error message
			default:
				cout << usage << endl;
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != 0) {
		cerr << usage << endl << "Incorrect number of arguments." << endl;
		return EXIT_FAILURE;
	}

	//**************************************************************************
	#pragma mark  Read input and set up corresponding SPGR & SSFP lists
	//**************************************************************************
	Sequences sequences(scale);
	// Build a Functor here so we can query number of parameters etc.
	cout << "Using " << to_string(pools) << " component sequences." << endl;
	vector<MultiArray<complex<float>, 4>> signalVols;
	Nifti::Header hdr = parseInput(sequences, signalVols);
	checkHeaders(hdr, {maskFile, f0File, B1File});
	//**************************************************************************
	#pragma mark Allocate memory and set up boundaries.
	//**************************************************************************
	MultiArray<float, 4> paramsVols(hdr.matrix(), PoolInfo::nParameters(pools));
	MultiArray<float, 4> ResidsVols(hdr.matrix(), sequences.size());;
	MultiArray<float, 3> ResVol(hdr.matrix());
	
	ArrayXd threshes(PoolInfo::nParameters(pools)); threshes.setConstant(0.05);
	ArrayXXd bounds = PoolInfo::Bounds(pools, tesla, sequences.minTR());
	if (tesla == FieldStrength::User) {
		if (prompt) cout << "Enter parameter pairs (low then high)" << endl;
		for (size_t i = 0; i < PoolInfo::nParameters(pools) - 1; i++) {
			if (prompt) cout << PoolInfo::Names(pools)[i] << ": " << flush;
			cin >> bounds(i, 0) >> bounds(i, 1);
		}
	}
	if (f0fit == OffRes::FitSym) {
		bounds(PoolInfo::nParameters(pools) - 1, 0) = 0.;
	}
	ArrayXd weights(sequences.size()); weights.setOnes();
	if (verbose) {
		cout << sequences;
		cout << "Bounds:" << endl <<  bounds.transpose() << endl;
		ofstream boundsFile(outPrefix + "bounds.txt");
		for (size_t p = 0; p < PoolInfo::nParameters(pools); p++) {
			boundsFile << PoolInfo::Names(pools)[p] << "\t" << bounds.row(p) << endl;
		}
		boundsFile.close();
	}
	
	//**************************************************************************
	#pragma mark Do the fitting
	//**************************************************************************
	if (stop_slice > hdr.dim(3))
		stop_slice = hdr.dim(3);
    time_t procStart = time(NULL);
	char theTime[512];
	strftime(theTime, 512, "%H:%M:%S", localtime(&procStart));
	cout << "Started processing at " << theTime << endl;
	for (size_t k = start_slice; k < stop_slice; k++) {
		if (verbose) cout << "Processing slice " << k << "..." << flush;
		atomic<int> voxCount{0};
		clock_t loopStart = clock();

		function<void (const size_t, const size_t)>
		processVox = [&] (const size_t i, const size_t j) {
			if (!maskFile || maskVol[{i,j,k}]) {
				ArrayXcd signal = sequences.loadSignals(signalVols, i, j, k, flipData);
				ArrayXXd localBounds = bounds;
				if (f0fit == OffRes::Map) {
					localBounds.row(PoolInfo::nParameters(pools) - 1).setConstant(f0Vol[{i,j,k}]);
				}
				if (scale == Scale::None) {
					localBounds(0, 0) = 0.;
					localBounds(0, 1) = signal.abs().maxCoeff() * 25;
				}
				double B1 = B1File ? B1Vol[{i,j,k}] : 1.;
				DESPOTFunctor func(sequences, pools, signal, B1, fitComplex, false);
				RegionContraction<DESPOTFunctor> rc(func, localBounds, weights, threshes,
													samples, retain, contract, expand, (voxI != 0));
				ArrayXd params(PoolInfo::nParameters(pools));
				rc.optimise(params); // Add the voxel number to the time to get a decent random seed
				paramsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = params.cast<float>();
				ResidsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = rc.residuals().cast<float>();
				ResVol[{i,j,k}] = static_cast<float>(rc.SoS());
				if ((rc.status() == RCStatus::Converged) || (rc.status() == RCStatus::IterationLimit)) {
					voxCount++;
				}
			}
		};
		if (voxI == 0) {
			threads.for_loop2(processVox, hdr.dim(1), hdr.dim(2));
			if (threads.interrupted())
				break;
		} else {
			processVox(voxI, voxJ);
			return EXIT_SUCCESS;
		}
		if (verbose) {
			clock_t loopEnd = clock();
			if (voxCount > 0)
				cout << voxCount << " unmasked voxels, CPU time per voxel was "
				          << ((loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC)) << " s, ";
			cout << "finished." << endl;
		}
		if (threads.interrupted())
			break;
	}
	time_t procEnd = time(NULL);
	strftime(theTime, 512, "%H:%M:%S", localtime(&procEnd));
	cout << "Finished processing at " << theTime << ". Run-time was " 
		 << difftime(procEnd, procStart) << " s." << endl;
	// Residuals can only be written here if we want them to go in a 4D gzipped file
	outPrefix = outPrefix + to_string(pools) + "C_";
	hdr.setDim(4, 1);
	hdr.setDatatype(Nifti::DataType::FLOAT32);
	hdr.description = version;
	hdr.intent = Nifti::Intent::Estimate;
	size_t start = (scale == Scale::None) ? 0 : 1;
	for (size_t p = start; p < PoolInfo::nParameters(pools); p++) {
		hdr.intent_name = PoolInfo::Names(pools).at(p);
		Nifti::File file(hdr, outPrefix + PoolInfo::Names(pools).at(p) + "" + OutExt());
		auto param = paramsVols.slice<3>({0,0,0,p},{-1,-1,-1,0});
		file.writeVolumes(param.begin(), param.end());
		file.close();
	}
	hdr.intent_name = "Fractional Residual";
	Nifti::File SoS(hdr, outPrefix + "residual" + OutExt());
	SoS.writeVolumes(ResVol.begin(), ResVol.end());
	SoS.close();
	if (all_residuals) {
		hdr.setDim(4, static_cast<int>(sequences.size()));
		hdr.intent_name = "Residuals";
		Nifti::File res(hdr, outPrefix + "residuals" + OutExt());
		res.writeVolumes(ResidsVols.begin(), ResidsVols.end());
		res.close();
	}
	cout << "Finished writing data." << endl;
	
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

