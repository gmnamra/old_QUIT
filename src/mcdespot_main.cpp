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

#include "Nifti/Nifti.h"
#include "DESPOT_Functors.h"
#include "ThreadPool.h"
#include "RegionContraction.h"

#ifdef AGILENT
	#include "procpar.h"
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
	--PD, -P n/i/g    : Proton Density mode. Default is n (normalise).\n\
	--finite, -f      : Use Finite Pulse Length correction.\n\
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

static auto B0fit = mcDESPOT::OffResMode::SingleSymmetric;
static auto components = Components::Two;
static auto tesla = mcDESPOT::FieldStrength::Three;
static auto scale = mcDESPOT::Scaling::NormToMean;
static size_t start_slice = 0, end_slice = numeric_limits<size_t>::max();
static int verbose = false, prompt = true, debug = false, converge_f = false,
           use_finite = false,
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
	{"scaling", required_argument, 0, 's'},
	{"tesla", required_argument, 0, 't'},
	{"B0", required_argument, 0, 'b'},
	{"cnv", no_argument, 0, 'C'},
	{"finite", no_argument, 0, 'f'},
	{"samples", required_argument, 0, 'm'},
	{"retain", required_argument, 0, 'r'},
	{"contract", required_argument, 0, 'c'},
	{"expand", required_argument, 0, 'e'},
	{"help", no_argument, 0, 'h'},
	{"1", no_argument, 0, '1'},
	{"2", no_argument, 0, '2'},
	{"3", no_argument, 0, '3'},
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
Nifti openAndCheck(const string &path, const Nifti &saved, const string &type);
Nifti openAndCheck(const string &path, const Nifti &saved, const string &type) {
	Nifti in(path, Nifti::Mode::Read);
	if (!(in.matchesSpace(saved))) {
		cerr << "Header for " << in.imagePath() << " does not match " << saved.imagePath() << endl;
		exit(EXIT_FAILURE);
	}
	if (verbose) cout << "Opened " << type << " image: " << in.imagePath() << endl;
	return in;
}

Nifti parseInput(vector<shared_ptr<SignalFunctor>> &sigs,
				       vector<Nifti > &signalFiles,
				       vector<Nifti > &B1Files,
				       vector<Nifti > &B0Files,
					   const mcDESPOT::OffResMode &B0fit,
					   const bool &use_finite);
Nifti parseInput(vector<shared_ptr<SignalFunctor>> &sigs,
				       vector<Nifti > &signalFiles,
				       vector<Nifti > &B1Files,
				       vector<Nifti > &B0Files,
					   const mcDESPOT::OffResMode &B0fit,
					   const bool &use_finite) {
	Nifti templateFile;
	string type, path;
	if (prompt) cout << "Specify next image type (SPGR/SSFP): " << flush;
	while (getline(cin, type) && (type != "END") && (type != "")) {
		if (type != "SPGR" && type != "SSFP") {
			cerr << "Unknown signal type: " << type << endl;
			exit(EXIT_FAILURE);
		}
		if (prompt) cout << "Enter image path: " << flush;
		getline(cin, path);
		if (signalFiles.size() == 0) {
			signalFiles.emplace_back(path, Nifti::Mode::Read);
			templateFile = Nifti(signalFiles.back(), 1); // Save header info for later
		} else {
			signalFiles.push_back(openAndCheck(path, templateFile, type));
		}
		if (type == "SPGR")
			sigs.emplace_back(parseSPGR(signalFiles.back(), prompt, components, use_finite));
		else
			sigs.emplace_back(parseSSFP(signalFiles.back(), prompt, components, use_finite));
		
		if (prompt) cout << "Enter B1 Map Path (Or NONE): " << flush;
		getline(cin, path);
		if ((path != "NONE") && (path != "")) {
			B1Files.push_back(openAndCheck(path, templateFile, "B1"));
		} else {
			B1Files.push_back(Nifti());
		}
		
		if ((type == "SSFP") && (B0fit == mcDESPOT::OffResMode::Map)) {
			if (prompt)
				cout << "Enter path to B0 map: " << flush;
			getline(cin, path);
			B0Files.push_back(openAndCheck(path, templateFile, "B0"));
		} else {
			B0Files.push_back(Nifti());
		}
		// Print message ready for next loop
		if (prompt) cout << "Specify next image type (SPGR/SSFP, END to finish input): " << flush;
	}
	if (sigs.size() == 0) {
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
	Nifti maskFile, templateFile;
	vector<double> maskData(0);
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "123hvpCft:b:m:o:s:w:n:r:c:e:i:j:d", long_options, &indexptr)) != -1) {
		switch (c) {
			case '1': components = Components::One; break;
			case '2': components = Components::Two; break;
			case '3': components = Components::Three; break;
			case 'C': converge_f = true; break;
			case 'f': use_finite = true; break;
			case 'i': voxI = atoi(optarg); break;
			case 'j': voxJ = atoi(optarg); break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
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
			case 's':
				switch (atoi(optarg)) {
					case 1 : scale = DESPOT2FM::Scaling::Global; break;
					case 2 : scale = DESPOT2FM::Scaling::NormToMean; break;
					default:
						cout << "Invalid scaling mode: " + to_string(atoi(optarg)) << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 'w': weighting = atof(optarg); break;
			case 'n': samples  = atoi(optarg); break;
			case 'r': retain   = atoi(optarg); break;
			case 'c': contract = atoi(optarg); break;
			case 'e': expand   = atof(optarg); break;
			case 'b':
				switch (*optarg) {
					case '0' : B0fit = mcDESPOT::OffResMode::Map; break;
					case '1' : B0fit = mcDESPOT::OffResMode::Single; break;
					default:
						cout << "Invalid B0 Mode." << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 't':
				switch (*optarg) {
					case '3': tesla = mcDESPOT::FieldStrength::Three; break;
					case '7': tesla = mcDESPOT::FieldStrength::Seven; break;
					case 'u': tesla = mcDESPOT::FieldStrength::Unknown; break;
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
	vector<shared_ptr<SignalFunctor>> sigs;
	vector<Nifti > signalFiles, B1Files, B0Files;
	templateFile = parseInput(sigs, signalFiles, B1Files, B0Files, B0fit, use_finite);
	if ((maskData.size() > 0) && !(maskFile.matchesSpace(templateFile))) {
		cerr << "Mask file has different dimensions/transform to input data." << endl;
		exit(EXIT_FAILURE);
	}
	//**************************************************************************
	#pragma mark Allocate memory and set up boundaries.
	// Use if files are open to indicate default values should be used -
	// 0 for B0, 1 for B1
	//**************************************************************************
	size_t voxelsPerSlice = templateFile.voxelsPerSlice();
	size_t voxelsPerVolume = templateFile.voxelsPerVolume();
	
	vector<vector<double>> signalVolumes(signalFiles.size()),
	                 B1Volumes(signalFiles.size()),
				     B0Volumes(signalFiles.size());
	for (size_t i = 0; i < signalFiles.size(); i++) {
		signalVolumes[i].resize(voxelsPerSlice * sigs.at(i)->size());
		if (B1Files[i].isOpen()) B1Volumes[i].resize(voxelsPerSlice);
		if (B0Files[i].isOpen()) B0Volumes[i].resize(voxelsPerSlice);
	}
	
	// Build a Functor here so we can query number of parameters etc.
	cout << "Using " << mcDESPOT::to_string(components) << " component model." << endl;
	mcDESPOT mcd(components, sigs, tesla, B0fit, scale, use_finite);
	outPrefix = outPrefix + mcDESPOT::to_string(components) + "C_";
	ArrayXd weights(mcd.values()); weights.setConstant(1.0);
	ArrayXd threshes(mcd.inputs()); threshes.setConstant(0.05);
	if (converge_f)
		threshes = mcd.defaultThresholds();
	templateFile.setDim(4, 1);
	templateFile.setDatatype(Nifti::DataType::FLOAT32);
	
	vector<Nifti> paramsFiles(mcd.inputs(), templateFile);
	vector<Nifti> midpFiles(mcd.inputs(), templateFile);
	vector<Nifti> widthFiles(mcd.inputs(), templateFile);
	Nifti contractFile(templateFile);
	vector<vector<double>> paramsData(mcd.inputs());
	vector<vector<double>> midpData(mcd.inputs());
	vector<vector<double>> widthData(mcd.inputs());
	vector<size_t> contractData(voxelsPerSlice);
	for (int i = 0; i < mcd.inputs(); i++) {
		paramsData.at(i).resize(voxelsPerSlice);
		paramsFiles.at(i).open(outPrefix + mcd.names()[i] + ".nii.gz", Nifti::Mode::Write);
		if (debug) {
			midpData.at(i).resize(voxelsPerSlice);
			midpFiles.at(i).open(outPrefix + mcd.names()[i] + "_mid.nii.gz", Nifti::Mode::Write);
			widthData.at(i).resize(voxelsPerSlice);
			widthFiles.at(i).open(outPrefix + mcd.names()[i] + "_width.nii.gz", Nifti::Mode::Write);
		}
	}
	if (debug)
		contractFile.open(outPrefix + "n_contract.nii.gz", Nifti::Mode::Write);
	
	vector<vector<double>> residualData(mcd.values());
	for (size_t i = 0; i < residualData.size(); i ++)
		residualData.at(i).resize(voxelsPerVolume);
	Nifti residualFile(templateFile);
	residualFile.setDim(4, static_cast<int>(mcd.values()));
	residualFile.open(outPrefix + mcDESPOT::to_string(components) + "_residuals.nii.gz", Nifti::Mode::Write);
	
	ArrayXXd bounds = mcd.defaultBounds();
	if (tesla == mcDESPOT::FieldStrength::Unknown) {
		if (prompt) cout << "Enter parameter pairs (low then high)" << endl;
		for (size_t i = 0; i < mcd.nP(); i++) {
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
	if (end_slice > templateFile.dim(3))
		end_slice = templateFile.dim(3);
	
	signal(SIGINT, int_handler);	// If we've got here there's actually allocated data to save
	
    time_t procStart = time(NULL);
	char theTime[512];
	strftime(theTime, 512, "%H:%M:%S", localtime(&procStart));
	cout << "Started processing at " << theTime << endl;
	for (size_t slice = start_slice; slice < end_slice; slice++)
	{
		if (verbose) cout << "Reading data for slice " << slice << "..." << flush;
		atomic<int> voxCount{0};
		const size_t sliceOffset = slice * voxelsPerSlice;
		
		// Read data for slices
		for (size_t i = 0; i < signalFiles.size(); i++) {
			signalFiles[i].readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, signalVolumes[i]);
			if (B1Files[i].isOpen()) B1Files[i].readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, B1Volumes[i]);
			if (B0Files[i].isOpen()) B0Files[i].readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, B0Volumes[i]);
		}
		if (verbose) cout << "processing..." << endl;
		clock_t loopStart = clock();
		function<void (const size_t&)> processVox = [&] (const size_t &vox)
		{
			mcDESPOT localf(mcd); // Take a thread local copy so we can change info/signals
			ArrayXd params(localf.inputs()), residuals(localf.values()),
					width(localf.inputs()), midp(localf.inputs());
			size_t c = 0;
			width.setZero(); midp.setZero(); params.setZero(); residuals.setZero();
			if ((maskData.size() == 0) || (maskData[sliceOffset + vox] > 0.)) {
				voxCount++;
				vector<VectorXd> signals(signalFiles.size());
				for (size_t i = 0; i < signalFiles.size(); i++) {
					for (size_t j = 0; j < localf.signal(i)->size(); j++) {
						localf.actual(i)(j) = signalVolumes[i][voxelsPerSlice*j + vox];
					}
					if (scale == mcDESPOT::Scaling::NormToMean)
						localf.actual(i) /= localf.actual(i).mean();
					
					if (B0fit == mcDESPOT::OffResMode::Map) {
						localf.m_f0 = B0Files[i].isOpen() ? B0Volumes[i][vox] : 0.;
					}
					localf.m_B1 = B1Files[i].isOpen() ? B1Volumes[i][vox] : 1.;
				}
				// Add the voxel number to the time to get a decent random seed
				size_t rSeed = time(NULL) + vox;
				RegionContraction<mcDESPOT> rc(localf, bounds, weights, threshes,
											 samples, retain, contract, expand, (voxI != -1));
				rc.optimise(params, rSeed);
				residuals = rc.residuals();
				if (debug) {
					c = rc.contractions();
					width = rc.width();
					midp = rc.midPoint();
				}
			}
			if (debug)
				contractData.at(vox) = c;
			for (size_t p = 0; p < paramsData.size(); p++) {
				paramsData.at(p).at(vox) = params[p];
				if (debug) {
					widthData.at(p).at(vox) = width(p);
					midpData.at(p).at(vox) = midp(p);
				}
			}
			for (int i = 0; i < residuals.size(); i++) {
				residualData.at(i).at(slice * voxelsPerSlice + vox) = residuals[i];
			}
			if (voxI != -1) {
				cout << "Final: " << params.transpose() << endl;
			}
		};
		if (voxI == -1)
			threads.for_loop(processVox, voxelsPerSlice);
		else {
			size_t voxInd = templateFile.dim(1) * voxJ + voxI;
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
		
		for (size_t p = 0; p < paramsFiles.size(); p++) {
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
	for (size_t r = 0; r < residualData.size(); r++)
		residualFile.writeSubvolume(0, 0, 0, r, -1, -1, -1, r+1, residualData.at(r));
	cout << "Finished writing data." << endl;
	return EXIT_SUCCESS;
}

