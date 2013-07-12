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

#include <time.h>
#include <getopt.h>
#include <signal.h>
#include <iostream>
#include <atomic>
#include <Eigen/Dense>

#include "NiftiImage.h"
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
	         1        : Fit one B0 value to all scans.\n\
	         2        : Fit a B0 value to each scan individually.\n\
	         3        : Fit one, bounded, B0 value to all scans.\n\
	         4        : Fit a bounded B0 value to each scan.\n\
	--tesla, -t 3     : Boundaries suitable for 3T (default)\n\
	            7     : Boundaries suitable for 7T \n\
	            u     : User specified boundaries from stdin.\n"
};

static int verbose = false, prompt = true,
           normalise = false,
		   B0Mode = mcDESPOT::B0_Single, finiteRF = false,
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
vector<NiftiImage> paramsHdrs;
NiftiImage residualHdr;
vector<vector<double>> paramsData;
vector<vector<double>> residualData;

void int_handler(int sig);
void int_handler(int sig)
{
	fprintf(stdout, "Processing terminated. Writing currently processed data.\n");
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
NiftiImage *openAndCheck(const string &path, const NiftiImage *saved, const string &type) {
	NiftiImage *in = new NiftiImage(path, NiftiImage::READ);
	if (!(in->matchesSpace(*saved))) {
		cerr << "Header for " << in->imagePath() << " does not match " << saved->imagePath() << endl;
		exit(EXIT_FAILURE);
	}
	if (verbose) cout << "Opened " << type << " image: " << in->imagePath() << endl;
	return in;
}

NiftiImage *parseInput(vector<mcDESPOT::SignalType> &signalTypes, vector<VectorXd> &angles,
                       vector<DESPOTConstants> &consts,
				       vector<NiftiImage *> &signalFiles,
				       vector<NiftiImage *> &B1_files,
				       vector<NiftiImage *> &B0_loFiles,
					   vector<NiftiImage *> &B0_hiFiles,
					   const int &B0Mode, const bool &finiteRF);
NiftiImage *parseInput(vector<mcDESPOT::SignalType> &signalTypes, vector<VectorXd> &angles,
                       vector<DESPOTConstants> &consts,
				       vector<NiftiImage *> &signalFiles,
				       vector<NiftiImage *> &B1_files,
				       vector<NiftiImage *> &B0_loFiles,
					   vector<NiftiImage *> &B0_hiFiles,
					   const int &B0Mode, const bool &finiteRF) {
	NiftiImage *inHdr, *savedHeader;
	string type, path;
	if (prompt) cout << "Specify next image type (SPGR/SSFP): " << flush;
	while (getline(cin, type) && (type != "END") && (type != "")) {
		
		if (type == "SPGR")
			signalTypes.push_back(mcDESPOT::SignalSPGR);
		else if (type == "SSFP")
			signalTypes.push_back(mcDESPOT::SignalSSFP);
		else {
			cerr << "Unknown signal type: " << type << endl;
			exit(EXIT_FAILURE);
		}
		if (prompt) cout << "Enter image path: " << flush;
		getline(cin, path);
		if (signalFiles.size() == 0) {
			savedHeader = new NiftiImage(path, NiftiImage::READ_HEADER);
		}
		inHdr = openAndCheck(path, savedHeader, type);
		signalFiles.push_back(inHdr);
		double inTR = 0., inTrf = 0., inPhase = 0., inTE = 0.;
		VectorXd inAngles(inHdr->dim(4));
		#ifdef HAVE_NRECON
		ParameterList pars;
		if (ReadProcpar(inHdr->basePath() + ".procpar", pars)) {
			inTR = RealValue(pars, "tr");
			for (int i = 0; i < inAngles.size(); i++)
				inAngles[i] = RealValue(pars, "flip1", i);
			if (signalTypes.back() == mcDESPOT::SignalSSFP) {
				inPhase = RealValue(pars, "rfphase");
			}
			if (finiteRF) {
				inTrf = RealValue(pars, "p1") / 1e6; // p1 is in microseconds
				inTE = RealValue(pars, "te");
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
			if (signalTypes.back() == mcDESPOT::SignalSSFP) {
				if (prompt) cout << "Enter SSFP Phase-Cycling (degrees): " << flush;
				cin >> inPhase;
			}
			if (finiteRF) {
				cout << "Enter RF Pulse Length (seconds): " << flush;
				cin >> inTrf;
				cout << "Enter TE (seconds): " << flush;
				cin >> inTE;
			}
			getline(cin, path); // Just to eat the newline
		}
		consts.push_back( { inTR, inTrf, inPhase * M_PI / 180., 0., 0. } );
		angles.push_back(inAngles * M_PI / 180.);
		
		inHdr = NULL;
		if (prompt) cout << "Enter B1 Map Path (Or NONE): " << flush;
		getline(cin, path);
		if ((path != "NONE") && (path != "")) {
			inHdr = openAndCheck(path, savedHeader, "B1");
		}
		B1_files.push_back(inHdr);
		
		inHdr = NULL;
		if ((signalTypes.back() == mcDESPOT::SignalSSFP) &&
		    ((B0Mode == mcDESPOT::B0_Map) || (B0Mode == mcDESPOT::B0_Bounded) || (B0Mode == mcDESPOT::B0_MultiBounded))) {
			if (prompt && (B0Mode == mcDESPOT::B0_Map)) cout << "Enter path to B0 map: " << flush;
			else if (prompt) cout << "Enter path to low B0 bound map: " << flush;
			getline(cin, path);
			inHdr = openAndCheck(path, savedHeader, "B0");		
		}
		B0_loFiles.push_back(inHdr);
		
		inHdr = NULL;
		if ((signalTypes.back() == mcDESPOT::SignalSSFP) &&
		    ((B0Mode == mcDESPOT::B0_Bounded) || (B0Mode == mcDESPOT::B0_MultiBounded))) {
			if (prompt) cout << "Enter path to high B0 bound map: " << flush;
			getline(cin, path);
			inHdr = openAndCheck(path, savedHeader, "B0");
		}
		B0_hiFiles.push_back(inHdr);
		// Print message ready for next loop
		if (prompt) cout << "Specify next image type (SPGR/SSFP, END to finish input): " << flush;
	}
	if (signalTypes.size() == 0) {
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
	NiftiImage inHeader, *savedHeader;
	vector<double> maskData, PDData;
	bool haveMask = false, havePD = false;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvpt:b:m:o:nfw:s:r:c:e:i:j:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'i': voxI = atoi(optarg); break;
			case 'j': voxJ = atoi(optarg); break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				if (!inHeader.open(optarg, NiftiImage::READ)) {
					exit(EXIT_FAILURE);
				}
				maskData = inHeader.readVolume<double>(0);
				inHeader.close();
				haveMask = true;
				break;
			case 'M':
				cout << "Reading PD file " << optarg << endl;
				if (!inHeader.open(optarg, NiftiImage::READ)) {
					exit(EXIT_FAILURE);
				}
				PDData = inHeader.readVolume<double>(0);
				inHeader.close();
				havePD = true;
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
					case '0' : B0Mode = mcDESPOT::B0_Map; break;
					case '1' : B0Mode = mcDESPOT::B0_Single; break;
					case '2' : B0Mode = mcDESPOT::B0_Multi; break;
					case '3' : B0Mode = mcDESPOT::B0_Bounded; break;
					case '4' : B0Mode = mcDESPOT::B0_MultiBounded; break;
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
	vector<mcDESPOT::SignalType> signalTypes;
	vector<DESPOTConstants> consts;
	vector<VectorXd> angles;
	vector<NiftiImage *> signalFiles, B1_files, B0_loFiles, B0_hiFiles;
	savedHeader = parseInput(signalTypes, angles, consts, signalFiles, B1_files, B0_loFiles, B0_hiFiles, B0Mode, finiteRF);
	//**************************************************************************
	#pragma mark Allocate memory and set up boundaries.
	// Use NULL to indicate that default values should be used -
	// 0 for B0, 1 for B1
	//**************************************************************************
	int voxelsPerSlice = savedHeader->voxelsPerSlice();
	int voxelsPerVolume = savedHeader->voxelsPerVolume();
	
	vector<vector<double>> signalVolumes(signalFiles.size()),
	                 B1Volumes(signalFiles.size()),
				     B0LoVolumes(signalFiles.size()),
					 B0HiVolumes(signalFiles.size());
	for (int i = 0; i < signalFiles.size(); i++) {
		signalVolumes[i].resize(voxelsPerSlice * angles[i].size());
		if (B1_files[i]) B1Volumes[i].resize(voxelsPerSlice);
		if (B0_loFiles[i]) B0LoVolumes[i].resize(voxelsPerSlice);
		if (B0_hiFiles[i]) B0HiVolumes[i].resize(voxelsPerSlice);
	}
	
	cout << "Using " << components << " component model." << endl;
	nP = mcDESPOT::nP(components);
	nB0 = mcDESPOT::nB0(B0Mode, signalFiles.size());
	const vector<string> names = mcDESPOT::names(components);
	
	int totalSignals = 0;
	for (int i = 0; i < angles.size(); i++)
		totalSignals += angles[i].size();
	ArrayXd weights(totalSignals);
	size_t index = 0;
	for (int i = 0; i < signalTypes.size(); i++) {
		if (signalTypes[i] == mcDESPOT::SignalSPGR)
			weights.segment(index, angles[i].size()).setConstant(weighting);
		else
			weights.segment(index, angles[i].size()).setConstant(1.0);
		index += angles[i].size();
	}
	
	residualData.resize(totalSignals);
	for (int i = 0; i < residualData.size(); i ++)
		residualData[i].resize(voxelsPerVolume);
	residualHdr = *savedHeader;
	residualHdr.setDim(4, totalSignals); residualHdr.setDatatype(NIFTI_TYPE_FLOAT32);
	residualHdr.open(outPrefix + "MCD_" + to_string(components) + "c_" + "Residual.nii.gz", NiftiImage::WRITE);
	
	paramsData.resize(nP + nB0);
	savedHeader->setDim(4, 1);
	savedHeader->setDatatype(DT_FLOAT32);
	paramsHdrs.resize(nP + nB0, *savedHeader);
	
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
		paramsHdrs[i].open(outPrefix + "MCD_" + to_string(components) + "c_" + names[i] + ".nii.gz", NiftiImage::WRITE);
	}
	
	for (int i = 0; i < nB0; i++) {
		loBounds[nP + i] = -0.5 / consts[i].TR;
		hiBounds[nP + i] =  0.5 / consts[i].TR;
		paramsData[nP + i] .resize(voxelsPerSlice);
		paramsHdrs[nP + i].open(outPrefix + "MCD_" + to_string(components) + "c_B0_" + to_string(i) + ".nii.gz", NiftiImage::WRITE);
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
    time_t procStart = time(NULL);
	if ((start_slice < 0) || (start_slice >= savedHeader->dim(3)))
		start_slice = 0;
	if ((end_slice < 0) || (end_slice > savedHeader->dim(3)))
		end_slice = savedHeader->dim(3);
	signal(SIGINT, int_handler);	// If we've got here there's actually allocated data to save
	cout << "Starting processing." << endl;
	
	for (slice = start_slice; slice < end_slice; slice++)
	{
		if (verbose) cout << "Reading data for slice " << slice << "..." << flush;
		atomic<int> voxCount{0};
		const int sliceOffset = slice * voxelsPerSlice;
		
		// Read data for slices
		for (int i = 0; i < signalFiles.size(); i++) {
			signalFiles[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, signalVolumes[i]);
			if (B1_files[i]) B1_files[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, B1Volumes[i]);
			if (B0_loFiles[i]) B0_loFiles[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, B0LoVolumes[i]);
			if (B0_hiFiles[i]) B0_hiFiles[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, B0HiVolumes[i]);
		}
		if (verbose) cout << "processing..." << flush;
		clock_t loopStart = clock();
		function<void (const int&)> processVox = [&] (const int &vox)
		{
			ArrayXd params(nP + nB0), residuals(totalSignals);
			params.setZero();
			residuals.setZero();
			if ((!haveMask || (maskData[sliceOffset + vox] > 0.)) &&
			    (!havePD || (PDData[sliceOffset + vox] > 0.))) {
				voxCount++;
				
				vector<VectorXd> signals(signalFiles.size());
				// Need local copy because B1 changes per voxel. If you assign
				// to the global 'consts' nothing happens but Clang doesn't give
				// an error
				vector<DESPOTConstants> localConsts = consts;
				for (int i = 0; i < signalFiles.size(); i++) {
					VectorXd temp(angles[i].size());
					for (int j = 0; j < angles[i].size(); j++) {
						temp[j] = signalVolumes[i][voxelsPerSlice*j + vox];
					}
					if (normalise)
						temp /= temp.mean();
					signals[i] = temp;
					if (B0Mode == mcDESPOT::B0_Map) {
						localConsts[i].B0 = B0_loFiles[i] ? B0LoVolumes[i][vox] : 0.;
					}
					localConsts[i].B1 = B1_files[i] ? B1Volumes[i][vox] : 1.;
				}
				// Add the voxel number to the time to get a decent random seed
				int rSeed = static_cast<int>(time(NULL)) + vox;
				ArrayXd localLo = loBounds, localHi = hiBounds;
				if (PDData.size()) {
					localLo(0) = (double)PDData[sliceOffset + vox];
					localHi(0) = (double)PDData[sliceOffset + vox];
				}
				if ((B0Mode == mcDESPOT::B0_Bounded) || (B0Mode == mcDESPOT::B0_MultiBounded)) {
					for (int b = 0; b < nB0; b++) {
						localLo(nP + b) = B0_loFiles[b] ? B0LoVolumes[b][vox] : 0.;
						localHi(nP + b) = B0_hiFiles[b] ? B0HiVolumes[b][vox] : 0.;
					}
				}
				
				if (!finiteRF) {
					mcDESPOT mcd(components, signalTypes, angles, signals, localConsts, B0Mode, normalise, (voxI > -1));
					residuals = regionContraction<mcDESPOT>(params, mcd, localLo, localHi, weights,
															samples, retain, contract, 0.05, expand, rSeed);
				} else {
					mcFinite mcd(components, signalTypes, angles, signals, localConsts, B0Mode, normalise, (voxI > -1));
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
			int voxInd = savedHeader->dim(1) * voxJ + voxI;
			processVox(voxInd);
			exit(0);
		}
		
		if (verbose) {
			clock_t loopEnd = clock();
			if (voxCount > 0)
				cout << voxCount << " unmasked voxels, CPU time per voxel was "
				          << ((loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC)) << " s, ";
			cout << "finished." << endl;
		}
		
		for (int p = 0; p < nP; p++)
			paramsHdrs[p].writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, paramsData[p]);
		for (int b = 0; b < nB0; b++)
			paramsHdrs[nP + b].writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, paramsData[nP + b]);
	}
    time_t procEnd = time(NULL);
    struct tm *localEnd = localtime(&procEnd);
	char theTime[512];
    strftime(theTime, 512, "%H:%M:%S", localEnd);
	cout << "Finished processing at " << theTime << ". Run-time was " 
	          << difftime(procEnd, procStart) << " s." << endl;
	
	// Clean up memory and close files (automatically done in destructor)
	// Residuals can only be written here if we want them to go in a 4D gzipped file
	for (int r = 0; r < totalSignals; r++) {
		residualHdr.writeSubvolume(0, 0, 0, r, -1, -1, -1, r+1, residualData[r]);
	}
	residualHdr.close();
	cout << "Finished writing data." << endl;
	return EXIT_SUCCESS;
}

