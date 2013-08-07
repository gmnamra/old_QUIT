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

#include <getopt.h>
#include <signal.h>
#include <iostream>
#include <atomic>
#include <chrono>
#include <Eigen/Dense>

#include "Nifti.h"
#include "DESPOT_Functors.h"
#include "ThreadPool.h"
#include "RegionContraction.h"

#ifdef HAVE_NRECON
	#include "procpar.h"
	using namespace Recon;
#endif

#ifdef USE_MCFINITE
	#define mcType mcFinite
#else
	#define mcType mcDESPOT
#endif

using namespace std;
using namespace Eigen;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string credit {
"mcdespot - written by tobias.wood@kcl.ac.uk. \n\
Acknowledgements greatfully received, grant discussions welcome."
};

const string usage {
"Usage is: mcdespot [options]\n\
\n\
The program will prompt for input (unless --no-prompt specified)\n\
\n\
All times (TR) are in SECONDS. All angles are in degrees.\n\
\n\
Options:\n\
	--help, -h		  : Print this message.\n\
	--verbose, -v     : Print extra information.\n\
	--no-prompt, -p   : Don't print prompts for input.\n\
	--mask, -m file   : Mask input with specified file.\n\
	--out, -o path    : Add a prefix to the output filenames.\n\
	--1, --2, --3     : Use 1, 2 or 3 component model (default 2).\n\
	--start_slice n   : Only start processing at slice n.\n\
	--end_slice n     : Finish at slice n-1.\n\
	--normalise, -n   : Normalise signals to maximum (Ignore M0).\n\
	--samples, -s n   : Use n samples for region contraction (Default 5000).\n\
	--retain, -r  n   : Retain n samples for new boundary (Default 50).\n\
	--contract, -c n  : Contract a maximum of n times (Default 10).\n\
	--expand, -e n    : Re-expand boundary by percentage n (Default 0).\n\
	--B0, -b 0        : Read B0 values from map files.\n\
	         1        : (Default) Fit one B0 value to all scans.\n\
	         2        : Fit a B0 value to each scan individually.\n\
	         3        : Fit one bounded B0 value to all scans.\n\
	         4        : Fit a bounded B0 value to each scan.\n\
	--tesla, -t 3     : Boundaries suitable for 3T (default)\n\
	            7     : Boundaries suitable for 7T \n\
	            u     : User specified boundaries from stdin.\n"
};

static auto B0fit = mcDESPOT::OffResMode::Single;
static auto components = mcDESPOT::Components::Two;
static auto tesla = mcDESPOT::FieldStrength::Three;
static auto PD = mcDESPOT::PDMode::Normalise;
static int verbose = false, prompt = true, finiteRF = false,
           start_slice = -1, end_slice = -1, slice = 0,
		   samples = 5000, retain = 50, contract = 10,
           voxI = -1, voxJ = -1;
static double expand = 0., weighting = 1.0;
static string outPrefix;
static struct option long_options[] =
{
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'p'},
	{"start_slice", required_argument, 0, 'S'},
	{"end_slice", required_argument, 0, 'E'},
	{"PD", required_argument, 0, 'P'},
	{"finite", no_argument, &finiteRF, true},
	{"tesla", required_argument, 0, 't'},
	{"samples", required_argument, 0, 's'},
	{"retain", required_argument, 0, 'r'},
	{"contract", required_argument, 0, 'c'},
	{"expand", required_argument, 0, 'e'},
	{"help", no_argument, 0, 'h'},
	{"1", no_argument, 0, '1'},
	{"2", no_argument, 0, '2'},
	{"3", no_argument, 0, '3'},
	{"B0", no_argument, 0, 'b'},
	{0, 0, 0, 0}
};
//******************************************************************************
#pragma mark SIGTERM interrupt handler and Threads
//******************************************************************************
ThreadPool threads;
bool interrupt_received = false;
void int_handler(int sig);
void int_handler(int sig)
{
	cout << endl << "Stopping processing early." << endl;
	threads.stop();
	interrupt_received = true;
}

//******************************************************************************
#pragma mark Read in all required files and data from cin
//******************************************************************************
//Utility function
Nifti::File openAndCheck(const string &path, const Nifti::File &saved, const string &type) {
	Nifti::File in(path, Nifti::Modes::Read);
	if (!(in.matchesSpace(saved))) {
		cerr << "Header for " << in.imagePath() << " does not match " << saved.imagePath() << endl;
		exit(EXIT_FAILURE);
	}
	if (verbose) cout << "Opened " << type << " image: " << in.imagePath() << endl;
	return in;
}

Nifti::File parseInput(vector<DESPOTData> &data,
				       vector<Nifti::File > &signalFiles,
				       vector<Nifti::File > &B1_files,
				       vector<Nifti::File > &B0_loFiles,
					   vector<Nifti::File > &B0_hiFiles,
					   const mcDESPOT::OffResMode &B0fit, const bool &finiteRF);
Nifti::File parseInput(vector<DESPOTData> &data,
				       vector<Nifti::File > &signalFiles,
				       vector<Nifti::File > &B1_files,
				       vector<Nifti::File > &B0_loFiles,
					   vector<Nifti::File > &B0_hiFiles,
					   const mcDESPOT::OffResMode &B0fit, const bool &finiteRF) {
	Nifti::File savedHeader;
	string type, path;
	if (prompt) cout << "Specify next image type (SPGR/SSFP): " << flush;
	while (getline(cin, type) && (type != "END") && (type != "")) {
		if (type != "SPGR" && type != "SSFP") {
			cerr << "Unknown signal type: " << type << endl;
			exit(EXIT_FAILURE);
		}
		bool spoil = (type == "SPGR");
		if (prompt) cout << "Enter image path: " << flush;
		getline(cin, path);
		if (signalFiles.size() == 0) {
			savedHeader.open(path, Nifti::Modes::ReadHeader);
		}
		signalFiles.push_back(openAndCheck(path, savedHeader, type));
		double inTR = 0., inTrf = 0., inPhase = 0., inTE = 0.;
		VectorXd inAngles(signalFiles.back().dim(4));
		#ifdef HAVE_NRECON
		ParameterList pars;
		if (ReadProcpar(signalFiles.back().basePath() + ".procpar", pars)) {
			inTR = RealValue(pars, "tr");
			for (int i = 0; i < inAngles.size(); i++)
				inAngles[i] = RealValue(pars, "flip1", i);
			if (!spoil)
				inPhase = RealValue(pars, "rfphase");
			#ifdef USE_MCFINITE
				if (spoil)
					inTE = RealValue(pars, "te");
				inTrf = RealValue(pars, "p1") / 1e6; // p1 is in microseconds
			#endif
		} else
		#endif
		{
			if (prompt) cout << "Enter TR (seconds): " << flush;
			cin >> inTR;
			if (prompt) cout << "Enter " << inAngles.size() << " Flip-angles (degrees): " << flush;
			for (int i = 0; i < inAngles.size(); i++)
				cin >> inAngles[i];
			getline(cin, path); // Just to eat the newline
			if (!spoil) {
				if (prompt) cout << "Enter SSFP Phase-Cycling (degrees): " << flush;
				cin >> inPhase;
				getline(cin, path); // Just to eat the newline
			}
			#ifdef USE_MCFINITE
				if (spoil) {
					if (prompt) cout << "Enter TE (seconds): " << flush;
					cin >> inTE;
				}
				if (prompt) cout << "Enter RF Pulse Length (seconds): " << flush;
				cin >> inTrf;
				getline(cin, path); // Just to eat the newline
			#endif
		}
		data.emplace_back(inAngles.size(), spoil, inTR, inTrf, inTE, inPhase * M_PI / 180.);
		data.back().setFlip(inAngles * M_PI / 180.);
		
		if (prompt) cout << "Enter B1 Map Path (Or NONE): " << flush;
		getline(cin, path);
		if ((path != "NONE") && (path != "")) {
			B1_files.push_back(openAndCheck(path, savedHeader, "B1"));
		} else {
			B1_files.push_back(Nifti::File());
		}
		
		if ((!spoil) && ((B0fit == mcDESPOT::OffResMode::Map) || (B0fit == mcDESPOT::OffResMode::Bounded) || (B0fit == mcDESPOT::OffResMode::MultiBounded))) {
			if (prompt && (B0fit == mcDESPOT::OffResMode::Map))
				cout << "Enter path to B0 map: " << flush;
			else if (prompt)
				cout << "Enter path to low B0 bound map: " << flush;
			getline(cin, path);
			B0_loFiles.push_back(openAndCheck(path, savedHeader, "B0"));
		} else {
			B0_loFiles.push_back(Nifti::File());
		}
		
		if ((!spoil) && ((B0fit == mcDESPOT::OffResMode::Bounded) || (B0fit == mcDESPOT::OffResMode::MultiBounded))) {
			if (prompt) cout << "Enter path to high B0 bound map: " << flush;
			getline(cin, path);
			B0_hiFiles.push_back(openAndCheck(path, savedHeader, "B0"));
		} else {
			B0_hiFiles.push_back(Nifti::File());
		}
		// Print message ready for next loop
		if (prompt) cout << "Specify next image type (SPGR/SSFP, END to finish input): " << flush;
	}
	if (data.size() == 0) {
		cerr << "No input images specified." << endl;
		exit(EXIT_FAILURE);
	}
	return savedHeader;
}
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	#pragma mark Argument Processing
	//**************************************************************************
	cout << credit << endl;
	Eigen::initParallel();
	Nifti::File maskFile, savedHeader;
	vector<double> maskData(0);
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "123hvpt:b:m:o:P:fw:s:r:c:e:i:j:", long_options, &indexptr)) != -1) {
		switch (c) {
			case '1': components = mcType::Components::One; break;
			case '2': components = mcType::Components::Two; break;
			case '3': components = mcType::Components::Three; break;
			case 'i': voxI = atoi(optarg); break;
			case 'j': voxJ = atoi(optarg); break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Modes::Read);
				maskData = maskFile.readVolume<double>(0);
				maskFile.close();
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'v': verbose = true; break;
			case 'p': prompt = false; break;
			case 'S': start_slice = atoi(optarg); break;
			case 'E': end_slice = atoi(optarg); break;
			case 'P':
				switch (*optarg) {
					case 'n' : PD = mcType::PDMode::Normalise; break;
					case 'i' : PD = mcType::PDMode::Individual; break;
					case 'g' : PD = mcType::PDMode::Global; break;
					default:
						cout << "Invalid PD fitting mode." << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 'f': finiteRF = true; break;
			case 'w': weighting = atof(optarg); break;
			case 's': samples  = atoi(optarg); break;
			case 'r': retain   = atoi(optarg); break;
			case 'c': contract = atoi(optarg); break;
			case 'e': expand   = atof(optarg); break;
			case 'b':
				switch (*optarg) {
					case '0' : B0fit = mcType::OffResMode::Map; break;
					case '1' : B0fit = mcType::OffResMode::Single; break;
					case '2' : B0fit = mcType::OffResMode::Multi; break;
					case '3' : B0fit = mcType::OffResMode::Bounded; break;
					case '4' : B0fit = mcType::OffResMode::MultiBounded; break;
					default:
						cout << "Invalid B0 Mode." << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 't':
				switch (*optarg) {
					case '3': tesla = mcType::FieldStrength::Three; break;
					case '7': tesla = mcType::FieldStrength::Seven; break;
					case 'u': tesla = mcType::FieldStrength::Unknown; break;
					default:
						cout << "Unknown boundaries type " << optarg << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 0:
				// Just a flag
				break;
			case 'h':
			case '?': // getopt will print an error message
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
	vector<DESPOTData> data;
	vector<Nifti::File > signalFiles, B1_files, B0_loFiles, B0_hiFiles;
	savedHeader = parseInput(data, signalFiles, B1_files, B0_loFiles, B0_hiFiles, B0fit, finiteRF);
	if ((maskData.size() > 0) && !(maskFile.matchesSpace(savedHeader))) {
		cerr << "Mask file has different dimensions/transform to input data." << endl;
		exit(EXIT_FAILURE);
	}
	//**************************************************************************
	#pragma mark Allocate memory and set up boundaries.
	// Use if files are open to indicate default values should be used -
	// 0 for B0, 1 for B1
	//**************************************************************************
	int voxelsPerSlice = savedHeader.voxelsPerSlice();
	int voxelsPerVolume = savedHeader.voxelsPerVolume();
	
	vector<vector<double>> signalVolumes(signalFiles.size()),
	                 B1Volumes(signalFiles.size()),
				     B0LoVolumes(signalFiles.size()),
					 B0HiVolumes(signalFiles.size());
	for (int i = 0; i < signalFiles.size(); i++) {
		signalVolumes[i].resize(voxelsPerSlice * data[i].size());
		if (B1_files[i].isOpen()) B1Volumes[i].resize(voxelsPerSlice);
		if (B0_loFiles[i].isOpen()) B0LoVolumes[i].resize(voxelsPerSlice);
		if (B0_hiFiles[i].isOpen()) B0HiVolumes[i].resize(voxelsPerSlice);
	}
	
	cout << "Using " << mcType::to_string(components) << " component model." << endl;
	size_t nP = mcType::nP(components);
	size_t nB0 = mcType::nOffRes(B0fit, signalFiles.size());
	size_t nPD = mcType::nPD(PD, signalFiles.size());
	const vector<string> names = mcType::names(components);
	
	int totalSignals = 0;
	for (int i = 0; i < data.size(); i++)
		totalSignals += data[i].size();
	ArrayXd weights(totalSignals);
	size_t index = 0;
	for (int i = 0; i < data.size(); i++) {
		if (data[i].spoil)
			weights.segment(index, data[i].size()).setConstant(weighting);
		else
			weights.segment(index, data[i].size()).setConstant(1.0);
		index += data[i].size();
	}
	savedHeader.setDim(4, 1);
	savedHeader.setDatatype(DT_FLOAT32);
	vector<Nifti::File> paramsHdrs(nP + nB0 + nPD, savedHeader);
	vector<vector<double>> paramsData(nP + nB0 + nPD);
	
	vector<vector<double>> residualData(totalSignals);
	for (int i = 0; i < residualData.size(); i ++)
		residualData[i].resize(voxelsPerVolume);
	Nifti::File residualHdr(savedHeader);
	residualHdr.setDim(4, totalSignals);
	residualHdr.open(outPrefix + "MCD_" + mcType::to_string(components) + "c_" + "Residual.nii.gz", Nifti::Modes::Write);
	
	ArrayXXd bounds(nP + nB0 + nPD, 2);
	bounds.block(0, 0, nP, 2) = mcType::defaultBounds(components, tesla);
	if (prompt && tesla == mcType::FieldStrength::Unknown) {
		cout << "Enter parameter pairs (low then high)" << endl;
	}
	for (int i = 0; i < nP; i++) {
		if (tesla == mcType::FieldStrength::Unknown) {
			if (prompt) cout << names[i] << ": " << flush;
			cin >> bounds(i, 0) >> bounds(i, 1);
		}
		paramsData.at(i).resize(voxelsPerSlice);
		paramsHdrs.at(i).open(outPrefix + "MCD_" + mcType::to_string(components) + "c_" + names[i] + ".nii.gz", Nifti::Modes::Write);
	}
	
	for (int i = 0; i < nB0; i++) {
		bounds(nP + i, 0) = -0.5 / data[i].TR;
		bounds(nP + i, 1) =  0.5 / data[i].TR;
		paramsData.at(nP + i).resize(voxelsPerSlice);
		paramsHdrs.at(nP + i).open(outPrefix + "MCD_" + mcType::to_string(components) + "c_B0_" + to_string(i) + ".nii.gz", Nifti::Modes::Write);
	}
	for (int i = 0; i < nPD; i++) {
		bounds(nP + nB0 + i, 0) = INFINITY; // These will be set properly later in the main loop
		bounds(nP + nB0 + i, 1) = 0.;
		paramsData.at(nP + nB0 + i).resize(voxelsPerSlice);
		paramsHdrs.at(nP + nB0 + i).open(outPrefix + "MCD_" + mcType::to_string(components) + "c_PD_" + to_string(i) + ".nii.gz", Nifti::Modes::Write);
	}
	
	if (verbose) {
		cout << "Low bounds: " << bounds.col(0).transpose() << endl;
		cout << "Hi bounds:  " << bounds.col(1).transpose() << endl;
	}
	//**************************************************************************
	#pragma mark Do the fitting
	//**************************************************************************
	if ((start_slice < 0) || (start_slice >= savedHeader.dim(3)))
		start_slice = 0;
	if ((end_slice < 0) || (end_slice > savedHeader.dim(3)))
		end_slice = savedHeader.dim(3);
	signal(SIGINT, int_handler);	// If we've got here there's actually allocated data to save
	
    auto procStart = chrono::system_clock::now();
	time_t c_time = chrono::system_clock::to_time_t(procStart); // Still have to convert to c to use IO functions
	char theTime[512];
	strftime(theTime, 512, "%H:%M:%S", localtime(&c_time));
	cout << "Started processing at " << theTime << endl;
	for (slice = start_slice; slice < end_slice; slice++)
	{
		if (verbose) cout << "Reading data for slice " << slice << "..." << flush;
		atomic<int> voxCount{0};
		const int sliceOffset = slice * voxelsPerSlice;
		
		// Read data for slices
		for (size_t i = 0; i < signalFiles.size(); i++) {
			signalFiles[i].readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, signalVolumes[i]);
			if (B1_files[i].isOpen()) B1_files[i].readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, B1Volumes[i]);
			if (B0_loFiles[i].isOpen()) B0_loFiles[i].readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, B0LoVolumes[i]);
			if (B0_hiFiles[i].isOpen()) B0_hiFiles[i].readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, B0HiVolumes[i]);
		}
		if (verbose) cout << "processing..." << endl;
		auto start = chrono::steady_clock::now();
		function<void (const int&)> processVox = [&] (const int &vox)
		{
			ArrayXd params(nP + nB0 + nPD), residuals(totalSignals);
			params.setZero();
			residuals.setZero();
			if ((maskData.size() == 0) || (maskData[sliceOffset + vox] > 0.)) {
				voxCount++;
				vector<VectorXd> signals(signalFiles.size());
				// Need local copies because of per-voxel changes
				vector<DESPOTData> localData = data;
				ArrayXXd localBounds = bounds;
				for (size_t i = 0; i < signalFiles.size(); i++) {
					VectorXd sig(localData[i].size());
					for (size_t j = 0; j < localData[i].size(); j++) {
						sig(j) = signalVolumes[i][voxelsPerSlice*j + vox];
					}
					switch (PD) {
						case (mcType::PDMode::Normalise): sig /= sig.mean(); break;
						case (mcType::PDMode::Global):
							localBounds(nP + nB0 + i, 0) = min(localBounds(nP + nB0 + i, 0), sig.maxCoeff());
							localBounds(nP + nB0 + i, 1) = max(localBounds(nP + nB0 + i, 0), 100 * sig.maxCoeff());
						case (mcType::PDMode::Individual):
							localBounds(nP + nB0 + i, 0) = sig.maxCoeff();
							localBounds(nP + nB0 + i, 1) = 100 * sig.maxCoeff();
							break;
					}
					localData[i].setSignal(sig);
					if (B0fit == mcType::OffResMode::Map) {
						localData[i].f0_off = B0_loFiles[i].isOpen() ? B0LoVolumes[i][vox] : 0.;
					}
					localData[i].B1 = B1_files[i].isOpen() ? B1Volumes[i][vox] : 1.;
				}
				// Add the voxel number to the time to get a decent random seed
				int rSeed = static_cast<int>(time(NULL)) + vox;
				if ((B0fit == mcType::OffResMode::Bounded) || (B0fit == mcType::OffResMode::MultiBounded)) {
					for (int b = 0; b < nB0; b++) {
						localBounds(nP + b, 0) = B0_loFiles[b].isOpen() ? B0LoVolumes[b][vox] : 0.;
						localBounds(nP + b, 1) = B0_hiFiles[b].isOpen() ? B0HiVolumes[b][vox] : 0.;
					}
				}
				mcType mcd(components, localData, B0fit, PD, (voxI > -1));
				residuals = regionContraction<mcType>(params, mcd, localBounds, weights,
													  samples, retain, contract, 0.05, expand, rSeed);
			}
			for (int p = 0; p < (nP + nB0 + nPD); p++) {
				paramsData.at(p).at(vox) = params[p];
			}
			for (int i = 0; i < totalSignals; i++) {
				residualData.at(i).at(slice * voxelsPerSlice + vox) = residuals[i];
			}
		};
		if (voxI == -1)
			threads.for_loop(processVox, voxelsPerSlice);
		else {
			int voxInd = savedHeader.dim(1) * voxJ + voxI;
			processVox(voxInd);
			exit(0);
		}
		
		if (verbose) {
		    auto end = chrono::steady_clock::now();
			if (voxCount > 0)
				cout << voxCount << " unmasked voxels, CPU time per voxel was "
				     << chrono::duration_cast<chrono::milliseconds>(end-start).count() / voxCount << " ms." << endl;
		}
		
		for (int p = 0; p < (nP + nB0 + nPD); p++)
			paramsHdrs.at(p).writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, paramsData[p]);
		if (interrupt_received)
			break;
	}
    auto procEnd = chrono::system_clock::now();
	c_time = chrono::system_clock::to_time_t(procEnd);
	strftime(theTime, 512, "%H:%M:%S", localtime(&c_time));
	cout << "Finished processing at " << theTime << ". Run-time was "
	          << chrono::duration_cast<chrono::minutes>(procEnd - procStart).count() << " minutes." << endl;
	
	// Clean up memory and close files (automatically done in destructor)
	// Residuals can only be written here if we want them to go in a 4D gzipped file
	for (int r = 0; r < totalSignals; r++) {
		residualHdr.writeSubvolume(0, 0, 0, r, -1, -1, -1, r+1, residualData[r]);
	}
	residualHdr.close();
	cout << "Finished writing data." << endl;
	return EXIT_SUCCESS;
}

