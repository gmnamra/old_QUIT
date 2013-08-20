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
#include <signal.h>
#include <time.h>
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

static auto B0fit = mcType::OffResMode::Single;
static auto components = mcType::Components::Two;
static auto tesla = mcType::FieldStrength::Three;
static auto PD = mcType::PDMode::Normalise;
static int verbose = false, prompt = true, debug = false,
           start_slice = -1, end_slice = -1,
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

Nifti::File parseInput(vector<Info> &info,
				       vector<Nifti::File > &signalFiles,
				       vector<Nifti::File > &B1_files,
				       vector<Nifti::File > &B0_loFiles,
					   vector<Nifti::File > &B0_hiFiles,
					   const mcDESPOT::OffResMode &B0fit);
Nifti::File parseInput(vector<Info> &info,
				       vector<Nifti::File > &signalFiles,
				       vector<Nifti::File > &B1_files,
				       vector<Nifti::File > &B0_loFiles,
					   vector<Nifti::File > &B0_hiFiles,
					   const mcDESPOT::OffResMode &B0fit) {
	Nifti::File templateFile;
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
			templateFile.open(path, Nifti::Modes::ReadHeader);
		}
		signalFiles.push_back(openAndCheck(path, templateFile, type));
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
				inTrf = RealValue(pars, "p1") / 1.e6; // p1 is in microseconds
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
		info.emplace_back(inAngles * M_PI / 180., spoil, inTR, inTrf, inTE, inPhase * M_PI / 180.);
		
		if (prompt) cout << "Enter B1 Map Path (Or NONE): " << flush;
		getline(cin, path);
		if ((path != "NONE") && (path != "")) {
			B1_files.push_back(openAndCheck(path, templateFile, "B1"));
		} else {
			B1_files.push_back(Nifti::File());
		}
		
		if ((!spoil) && ((B0fit == mcDESPOT::OffResMode::Map) || (B0fit == mcDESPOT::OffResMode::Bounded) || (B0fit == mcDESPOT::OffResMode::MultiBounded))) {
			if (prompt && (B0fit == mcDESPOT::OffResMode::Map))
				cout << "Enter path to B0 map: " << flush;
			else if (prompt)
				cout << "Enter path to low B0 bound map: " << flush;
			getline(cin, path);
			B0_loFiles.push_back(openAndCheck(path, templateFile, "B0"));
		} else {
			B0_loFiles.push_back(Nifti::File());
		}
		
		if ((!spoil) && ((B0fit == mcDESPOT::OffResMode::Bounded) || (B0fit == mcDESPOT::OffResMode::MultiBounded))) {
			if (prompt) cout << "Enter path to high B0 bound map: " << flush;
			getline(cin, path);
			B0_hiFiles.push_back(openAndCheck(path, templateFile, "B0"));
		} else {
			B0_hiFiles.push_back(Nifti::File());
		}
		// Print message ready for next loop
		if (prompt) cout << "Specify next image type (SPGR/SSFP, END to finish input): " << flush;
	}
	if (info.size() == 0) {
		cerr << "No input images specified." << endl;
		exit(EXIT_FAILURE);
	}
	return templateFile;
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
	Nifti::File maskFile, templateFile;
	vector<double> maskData(0);
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "123hvpt:b:m:o:P:w:s:r:c:e:i:j:d", long_options, &indexptr)) != -1) {
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
			case 'd': debug = true; break;
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
	vector<Info> info;
	vector<Nifti::File > signalFiles, B1_files, B0_loFiles, B0_hiFiles;
	templateFile = parseInput(info, signalFiles, B1_files, B0_loFiles, B0_hiFiles, B0fit);
	if ((maskData.size() > 0) && !(maskFile.matchesSpace(templateFile))) {
		cerr << "Mask file has different dimensions/transform to input data." << endl;
		exit(EXIT_FAILURE);
	}
	//**************************************************************************
	#pragma mark Allocate memory and set up boundaries.
	// Use if files are open to indicate default values should be used -
	// 0 for B0, 1 for B1
	//**************************************************************************
	int voxelsPerSlice = templateFile.voxelsPerSlice();
	int voxelsPerVolume = templateFile.voxelsPerVolume();
	
	vector<vector<double>> signalVolumes(signalFiles.size()),
	                 B1Volumes(signalFiles.size()),
				     B0LoVolumes(signalFiles.size()),
					 B0HiVolumes(signalFiles.size());
	for (int i = 0; i < signalFiles.size(); i++) {
		signalVolumes[i].resize(voxelsPerSlice * info.at(i).nAngles());
		if (B1_files[i].isOpen()) B1Volumes[i].resize(voxelsPerSlice);
		if (B0_loFiles[i].isOpen()) B0LoVolumes[i].resize(voxelsPerSlice);
		if (B0_hiFiles[i].isOpen()) B0HiVolumes[i].resize(voxelsPerSlice);
	}
	
	// Build a Functor here so we can query number of parameters etc.
	cout << "Using " << mcType::to_string(components) << " component model." << endl;
	mcType mcd(components, info, tesla, B0fit, PD);
	outPrefix = outPrefix + mcType::to_string(components) + "CD_";
	
	ArrayXd weights(mcd.values());
	size_t index = 0;
	for (int i = 0; i < info.size(); i++) {
		if (info.at(i).spoil)
			weights.segment(index, info.at(i).nAngles()).setConstant(weighting);
		else
			weights.segment(index, info.at(i).nAngles()).setConstant(1.0);
		index += info.at(i).nAngles();
	}
	templateFile.setDim(4, 1);
	templateFile.setDatatype(DT_FLOAT32);
	
	vector<Nifti::File> paramsFiles(mcd.inputs(), templateFile);
	vector<Nifti::File> midpFiles(mcd.inputs(), templateFile);
	vector<Nifti::File> widthFiles(mcd.inputs(), templateFile);
	Nifti::File contractFile(templateFile);
	vector<vector<float>> paramsData(mcd.inputs());
	vector<vector<float>> midpData(mcd.inputs());
	vector<vector<float>> widthData(mcd.inputs());
	vector<size_t> contractData(voxelsPerSlice);
	for (int i = 0; i < mcd.inputs(); i++) {
		paramsData.at(i).resize(voxelsPerSlice);
		paramsFiles.at(i).open(outPrefix + mcd.names()[i] + ".nii.gz", Nifti::Modes::Write);
		if (debug) {
			midpData.at(i).resize(voxelsPerSlice);
			midpFiles.at(i).open(outPrefix + mcd.names()[i] + "_mid.nii.gz", Nifti::Modes::Write);
			widthData.at(i).resize(voxelsPerSlice);
			widthFiles.at(i).open(outPrefix + mcd.names()[i] + "_width.nii.gz", Nifti::Modes::Write);
		}
	}
	if (debug)
		contractFile.open(outPrefix + "_n_contract.nii.gz", Nifti::Modes::Write);
	
	vector<vector<double>> residualData(mcd.values());
	for (int i = 0; i < residualData.size(); i ++)
		residualData.at(i).resize(voxelsPerVolume);
	Nifti::File residualFile(templateFile);
	residualFile.setDim(4, static_cast<int>(mcd.values()));
	residualFile.open(outPrefix + mcType::to_string(components) + "CD_residuals.nii.gz", Nifti::Modes::Write);
	
	ArrayXXd bounds = mcd.defaultBounds();
	if (prompt && tesla == mcType::FieldStrength::Unknown) {
		cout << "Enter parameter pairs (low then high)" << endl;
		for (int i = 0; i < mcd.nP(); i++) {
			if (prompt) cout << mcd.names()[i] << ": " << flush;
			cin >> bounds(i, 0) >> bounds(i, 1);
		}
	}
	
	if (verbose) {
		cout << "Low bounds: " << bounds.col(0).transpose() << endl;
		cout << "Hi bounds:  " << bounds.col(1).transpose() << endl;
	}
	//**************************************************************************
	#pragma mark Do the fitting
	//**************************************************************************
	if ((start_slice < 0) || (start_slice >= templateFile.dim(3)))
		start_slice = 0;
	if ((end_slice < 0) || (end_slice > templateFile.dim(3)))
		end_slice = templateFile.dim(3);
	signal(SIGINT, int_handler);	// If we've got here there's actually allocated data to save
	
    time_t procStart = time(NULL);
	char theTime[512];
	strftime(theTime, 512, "%H:%M:%S", localtime(&procStart));
	cout << "Started processing at " << theTime << endl;
	for (int slice = start_slice; slice < end_slice; slice++)
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
		clock_t loopStart = clock();
		function<void (const int&)> processVox = [&] (const int &vox)
		{
			mcType localf(mcd); // Take a thread local copy so we can change info/signals
			ArrayXd params(localf.inputs()), residuals(localf.values()),
					width(localf.inputs()), midp(localf.inputs());
			size_t c = 0;
			width.setZero(); midp.setZero(); params.setZero(); residuals.setZero();
			if ((maskData.size() == 0) || (maskData[sliceOffset + vox] > 0.)) {
				voxCount++;
				vector<VectorXd> signals(signalFiles.size());
				// Need local copies because of per-voxel changes
				ArrayXXd localBounds = bounds;
				for (size_t i = 0; i < signalFiles.size(); i++) {
					for (size_t j = 0; j < localf.signal(i).size(); j++) {
						localf.signal(i)(j) = signalVolumes[i][voxelsPerSlice*j + vox];
					}
					if (PD == mcType::PDMode::Normalise) {
						localf.signal(i) /= localf.signal(i).mean();
					}
					if (B0fit == mcType::OffResMode::Map) {
						localf.info(i).f0_off = B0_loFiles[i].isOpen() ? B0LoVolumes[i][vox] : 0.;
					}
					localf.info(i).B1 = B1_files[i].isOpen() ? B1Volumes[i][vox] : 1.;
				}
				if ((B0fit == mcType::OffResMode::Bounded) || (B0fit == mcType::OffResMode::MultiBounded)) {
					for (int b = 0; b < localf.nOffRes(); b++) {
						localBounds(localf.nP() + b, 0) = B0_loFiles[b].isOpen() ? B0LoVolumes[b][vox] : 0.;
						localBounds(localf.nP() + b, 1) = B0_hiFiles[b].isOpen() ? B0HiVolumes[b][vox] : 0.;
					}
				}
				// Add the voxel number to the time to get a decent random seed
				int rSeed = static_cast<int>(time(NULL)) + vox;
				RegionContraction<mcType> rc(localf, localBounds, weights,
											 samples, retain, contract, 0.05, expand);
				rc.optimise(params, rSeed);
				residuals = rc.residuals();
				if (debug) {
					c = rc.contractions();
					width = rc.width();
					midp = rc.midPoint();
				}
			}
			if (debug)
				contractData.at(sliceOffset + vox) = c;
			for (int p = 0; p < paramsData.size(); p++) {
				paramsData.at(p).at(vox) = params[p];
				if (debug) {
					widthData.at(p).at(sliceOffset + vox) = width(p);
					midpData.at(p).at(sliceOffset + vox) = midp(p);
				}
			}
			for (int i = 0; i < residuals.size(); i++) {
				residualData.at(i).at(slice * voxelsPerSlice + vox) = residuals[i];
			}
		};
		if (voxI == -1)
			threads.for_loop(processVox, voxelsPerSlice);
		else {
			int voxInd = templateFile.dim(1) * voxJ + voxI;
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
		
		for (int p = 0; p < paramsFiles.size(); p++) {
			paramsFiles.at(p).writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, paramsData.at(p));
			if (debug) {
				midpFiles.at(p).writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, midpData.at(p));
				widthFiles.at(p).writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, widthData.at(p));
			}
		}
		if (debug)
			contractFile.writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, contractData);
		if (interrupt_received)
			break;
	}
    time_t procEnd = time(NULL);
    strftime(theTime, 512, "%H:%M:%S", localtime(&procEnd));
	cout << "Finished processing at " << theTime << ". Run-time was " 
	          << difftime(procEnd, procStart) << " s." << endl;
	// Residuals can only be written here if we want them to go in a 4D gzipped file
	for (int r = 0; r < residualData.size(); r++)
		residualFile.writeSubvolume(0, 0, 0, r, -1, -1, -1, r+1, residualData.at(r));
	cout << "Finished writing data." << endl;
	return EXIT_SUCCESS;
}

