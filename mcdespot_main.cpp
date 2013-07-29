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
#include <iomanip>
#include <Eigen/Dense>

#include "Nifti.h"
#include "DESPOT.h"
#include "DESPOT_Functors.h"
#include "RegionContraction.h"

#ifdef HAVE_NRECON
	#include "procpar.h"
	using namespace Recon;
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
	--M0 file         : M0 Map file.\n\
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

static mcDESPOT::B0Mode B0fit = mcDESPOT::B0Mode::Single;
static int verbose = false, prompt = true,
           normalise = false, finiteRF = false,
           start_slice = -1, end_slice = -1, slice = 0,
		   samples = 5000, retain = 50, contract = 10,
		   components = 2, tesla = 3,
           nP = 0, nB0 = 0, voxI = -1, voxJ = -1;
static double expand = 0., weighting = 1.0;
static string outPrefix;
static struct option long_options[] =
{
	{"M0", required_argument, 0, 'M'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'p'},
	{"start_slice", required_argument, 0, 'S'},
	{"end_slice", required_argument, 0, 'E'},
	{"normalise", no_argument, &normalise, true},
	{"finite", no_argument, &finiteRF, true},
	{"tesla", required_argument, 0, 't'},
	{"samples", required_argument, 0, 's'},
	{"retain", required_argument, 0, 'r'},
	{"contract", required_argument, 0, 'c'},
	{"expand", required_argument, 0, 'e'},
	{"help", no_argument, 0, 'h'},
	{"1", no_argument, &components, 1},
	{"2", no_argument, &components, 2},
	{"3", no_argument, &components, 3},
	{"B0", no_argument, 0, 'b'},
	{0, 0, 0, 0}
};
//******************************************************************************
#pragma mark SIGTERM interrupt handler - for ensuring data is saved on a ctrl-c
//******************************************************************************
vector<Nifti::File> paramsHdrs;
Nifti::File residualHdr;
vector<vector<double>> paramsData;
vector<vector<double>> residualData;

void int_handler(int sig);
void int_handler(int sig)
{
	cout << "Processing terminated. Writing currently processed data.\n" << endl;
	for (int p = 0; p < nP; p++) {
		paramsHdrs[p].writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, paramsData[p]);
		paramsHdrs[p].close();
	}
	for (int b = 0; b < nB0; b++) {
		paramsHdrs[nP + b].writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, paramsData[nP + b]);
		paramsHdrs[nP + b].close();
	}
	for (int r = 0; r < residualData.size(); r++) {
		residualHdr.writeSubvolume(0, 0, 0, r, -1, -1, -1, r+1, residualData[r]);
	}
	residualHdr.close();
	exit(EXIT_FAILURE);
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
					   const mcDESPOT::B0Mode &B0fit, const bool &finiteRF);
Nifti::File parseInput(vector<DESPOTData> &data,
				       vector<Nifti::File > &signalFiles,
				       vector<Nifti::File > &B1_files,
				       vector<Nifti::File > &B0_loFiles,
					   vector<Nifti::File > &B0_hiFiles,
					   const mcDESPOT::B0Mode &B0fit, const bool &finiteRF) {
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
			if (finiteRF) {
				if (spoil)
					inTE = RealValue(pars, "te");
				inTrf = RealValue(pars, "p1") / 1e6; // p1 is in microseconds
			}
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
			}
			if (finiteRF) {
				if (spoil) {
					if (prompt) cout << "Enter TE (seconds): " << flush;
					cin >> inTE;
				}
				if (prompt) cout << "Enter RF Pulse Length (seconds): " << flush;
				cin >> inTrf;

			}
			getline(cin, path); // Just to eat the newline
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
		
		if ((!spoil) && ((B0fit == mcDESPOT::B0Mode::Map) || (B0fit == mcDESPOT::B0Mode::Bounded) || (B0fit == mcDESPOT::B0Mode::MultiBounded))) {
			if (prompt && (B0fit == mcDESPOT::B0Mode::Map))
				cout << "Enter path to B0 map: " << flush;
			else if (prompt)
				cout << "Enter path to low B0 bound map: " << flush;
			getline(cin, path);
			B0_loFiles.push_back(openAndCheck(path, savedHeader, "B0"));
		} else {
			B0_loFiles.push_back(Nifti::File());
		}
		
		if ((!spoil) && ((B0fit == mcDESPOT::B0Mode::Bounded) || (B0fit == mcDESPOT::B0Mode::MultiBounded))) {
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
	Nifti::File maskFile, PDFile, savedHeader;
	vector<double> maskData, PDData;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvpt:b:m:o:nfw:s:r:c:e:i:j:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'i': voxI = atoi(optarg); break;
			case 'j': voxJ = atoi(optarg); break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Modes::Read);
				maskData = maskFile.readVolume<double>(0);
				break;
			case 'M':
				cout << "Reading PD file " << optarg << endl;
				PDFile.open(optarg, Nifti::Modes::Read);
				PDData = PDFile.readVolume<double>(0);
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'v': verbose = true; break;
			case 'p': prompt = false; break;
			case 'S': start_slice = atoi(optarg); break;
			case 'E': end_slice = atoi(optarg); break;
			case 'n': normalise = true; break;
			case 'f': finiteRF = true; break;
			case 'w': weighting = atof(optarg); break;
			case 's': samples  = atoi(optarg); break;
			case 'r': retain   = atoi(optarg); break;
			case 'c': contract = atoi(optarg); break;
			case 'e': expand   = atof(optarg); break;
			case 'b':
				switch (*optarg) {
					case '0' : B0fit = mcDESPOT::B0Mode::Map; break;
					case '1' : B0fit = mcDESPOT::B0Mode::Single; break;
					case '2' : B0fit = mcDESPOT::B0Mode::Multi; break;
					case '3' : B0fit = mcDESPOT::B0Mode::Bounded; break;
					case '4' : B0fit = mcDESPOT::B0Mode::MultiBounded; break;
					default:
						cout << "Invalid B0 Mode." << endl;
						exit(EXIT_FAILURE);
						break;
				}
				break;
			case 't':
				switch (*optarg) {
					case '3': tesla = 3; break;
					case '7': tesla = 7; break;
					case 'u': tesla = -1; break;
					default:
						cout << "Unknown boundaries type " << optarg << endl;
						exit(EXIT_FAILURE);
						break;
				}
				cout << "Using " << tesla << "T boundaries." << endl;
				break;
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
	
	cout << "Using " << components << " component model." << endl;
	nP = mcDESPOT::nP(components);
	nB0 = mcDESPOT::nB0(B0fit, signalFiles.size());
	const vector<string> names = mcDESPOT::names(components);
	
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
	
	residualData.resize(totalSignals);
	for (int i = 0; i < residualData.size(); i ++)
		residualData[i].resize(voxelsPerVolume);
	residualHdr = savedHeader;
	residualHdr.setDim(4, totalSignals);
	residualHdr.setDatatype(NIFTI_TYPE_FLOAT32);
	residualHdr.open(outPrefix + "MCD_" + to_string(components) + "c_" + "Residual.nii.gz", Nifti::Modes::Write);
	
	paramsData.resize(nP + nB0);
	savedHeader.setDim(4, 1);
	savedHeader.setDatatype(DT_FLOAT32);
	paramsHdrs.resize(nP + nB0, savedHeader);
	
	ArrayXd loBounds(nP + nB0), hiBounds(nP + nB0);
	if (tesla > 0) {
		loBounds.head(nP) = mcDESPOT::defaultLo(components, tesla);
		hiBounds.head(nP) = mcDESPOT::defaultHi(components, tesla);
	} else if (prompt) {
		cout << "Enter parameter pairs (low then high)" << endl;
	}
	for (int i = 0; i < nP; i++) {
		if (tesla <= 0) {
			if (prompt) cout << names[i] << ": " << flush;
			cin >> loBounds[i] >> hiBounds[i];
		}
		paramsData[i].resize(voxelsPerSlice);
		paramsHdrs[i].open(outPrefix + "MCD_" + to_string(components) + "c_" + names[i] + ".nii.gz", Nifti::Modes::Write);
	}
	
	for (int i = 0; i < nB0; i++) {
		loBounds[nP + i] = -0.5 / data[i].TR;
		hiBounds[nP + i] =  0.5 / data[i].TR;
		paramsData[nP + i] .resize(voxelsPerSlice);
		paramsHdrs[nP + i].open(outPrefix + "MCD_" + to_string(components) + "c_B0_" + to_string(i) + ".nii.gz", Nifti::Modes::Write);
	}
	// If normalising, don't bother fitting for PD
	if (normalise) {
		loBounds[0] = 1.;
		hiBounds[0] = 1.;
	}
	
	if (verbose) {
		cout << "Low bounds: " << loBounds.transpose() << endl;
		cout << "Hi bounds:  " << hiBounds.transpose() << endl;
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
	cout << "Starting processing at " << put_time(localtime(&c_time), "%F %T") << endl;
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
		if (verbose) cout << "processing..." << flush;
		auto start = chrono::steady_clock::now();
		function<void (const int&)> processVox = [&] (const int &vox)
		{
			ArrayXd params(nP + nB0), residuals(totalSignals);
			params.setZero();
			residuals.setZero();
			if ((!maskFile.isOpen() || (maskData[sliceOffset + vox] > 0.)) &&
			    (!PDFile.isOpen() || (PDData[sliceOffset + vox] > 0.))) {
				voxCount++;
				
				vector<VectorXd> signals(signalFiles.size());
				// Need local copy because B1 changes per voxel.
				vector<DESPOTData> localData = data;
				for (size_t i = 0; i < signalFiles.size(); i++) {
					VectorXd sig(localData[i].size());
					for (size_t j = 0; j < localData[i].size(); j++) {
						sig(j) = signalVolumes[i][voxelsPerSlice*j + vox];
					}
					if (normalise)
						sig /= sig.mean();
					localData[i].setSignal(sig);
					if (B0fit == mcDESPOT::B0Mode::Map) {
						localData[i].B0 = B0_loFiles[i].isOpen() ? B0LoVolumes[i][vox] : 0.;
					}
					localData[i].B1 = B1_files[i].isOpen() ? B1Volumes[i][vox] : 1.;
				}
				// Add the voxel number to the time to get a decent random seed
				int rSeed = static_cast<int>(time(NULL)) + vox;
				ArrayXd localLo = loBounds, localHi = hiBounds;
				if (PDData.size()) {
					localLo(0) = (double)PDData[sliceOffset + vox];
					localHi(0) = (double)PDData[sliceOffset + vox];
				}
				if ((B0fit == mcDESPOT::B0Mode::Bounded) || (B0fit == mcDESPOT::B0Mode::MultiBounded)) {
					for (int b = 0; b < nB0; b++) {
						localLo(nP + b) = B0_loFiles[b].isOpen() ? B0LoVolumes[b][vox] : 0.;
						localHi(nP + b) = B0_hiFiles[b].isOpen() ? B0HiVolumes[b][vox] : 0.;
					}
				}
				
				if (!finiteRF) {
					mcDESPOT mcd(components, localData, B0fit, normalise, (voxI > -1));
					residuals = regionContraction<mcDESPOT>(params, mcd, localLo, localHi, weights,
															samples, retain, contract, 0.05, expand, rSeed);
				} else {
					mcFinite mcd(components, localData, B0fit, normalise, (voxI > -1));
					residuals = regionContraction<mcFinite>(params, mcd, localLo, localHi, weights,
															samples, retain, contract, 0.05, expand, rSeed);
				}
			}
			for (int p = 0; p < nP; p++) {
				paramsData[p][vox] = params[p];
			}
			for (int b = 0; b < nB0; b++) {
				paramsData[nP + b][vox] = params[nP + b];
			}
			for (int i = 0; i < totalSignals; i++) {
				residualData[i][slice * voxelsPerSlice + vox] = residuals[i];
			}
		};
		if (voxI == -1)
			apply_for(voxelsPerSlice, processVox);
		else {
			int voxInd = savedHeader.dim(1) * voxJ + voxI;
			processVox(voxInd);
			exit(0);
		}
		
		if (verbose) {
		    auto end = chrono::steady_clock::now();
			if (voxCount > 0)
				cout << voxCount << " unmasked voxels, CPU time per voxel was "
				     << chrono::duration_cast<chrono::milliseconds>(end-start).count() << " ms." << endl;
		}
		
		for (int p = 0; p < nP; p++)
			paramsHdrs[p].writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, paramsData[p]);
		for (int b = 0; b < nB0; b++)
			paramsHdrs[nP + b].writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, paramsData[nP + b]);
	}
    auto procEnd = chrono::system_clock::now();
	c_time = chrono::system_clock::to_time_t(procEnd);
	cout << "Finished processing at " << put_time(localtime(&c_time), "%F %T") << ". Run-time was "
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

