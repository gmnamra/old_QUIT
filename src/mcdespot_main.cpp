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
	--no-prompt, -n   : Don't print prompts for input.\n\
	--1, --2, --3     : Use 1, 2 or 3 component model (default 3).\n\
	--mask, -m file   : Mask input with specified file.\n\
	--out, -o path    : Add a prefix to the output filenames.\n\
	--f0, -f 0        : Read f0 values from map files.\n\
			 1        : (Default) Fit one symmetric f0 value to all scans.\n\
			 2        : Fit an unsymmetric f0 value to all scans.\n\
	--start, -s n     : Only start processing at slice n.\n\
	--stop, -p n      : Finish at slice n-1.\n\
	--scale, -S 0     : Normalise signals to mean (default).\n\
	            1     : Fit a scaling factor/proton density.\n\
	--tesla, -t 3     : Boundaries suitable for 3T (default)\n\
	            7     : Boundaries suitable for 7T \n\
	            u     : User specified boundaries from stdin.\n\
	--finite, -F      : Use Finite Pulse Length correction.\n\
	--contract, -c n  : Read contraction settings from stdin (Will prompt).\n"
};

static auto f0fit = mcDESPOT::OffResMode::SingleSymmetric;
static auto components = Components::Two;
static auto tesla = mcDESPOT::FieldStrength::Three;
static auto scale = mcDESPOT::Scaling::NormToMean;
static size_t start_slice = 0, stop_slice = numeric_limits<size_t>::max();
static int verbose = false, prompt = true, extra = false, early_finish = false,
           use_finite = false,
		   samples = 5000, retain = 50, contract = 10,
           voxI = -1, voxJ = -1;
static double expand = 0., weighting = 1.0;
static string outPrefix;
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'p'},
	{"1", no_argument, 0, '1'},
	{"2", no_argument, 0, '2'},
	{"3", no_argument, 0, '3'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"f0", required_argument, 0, 'O'},
	{"start", required_argument, 0, 'S'},
	{"stop", required_argument, 0, 'P'},
	{"scale", required_argument, 0, 's'},
	{"tesla", required_argument, 0, 't'},
	{"finite", no_argument, 0, 'f'},
	{"contract", no_argument, 0, 'c'},
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
				       vector<Nifti > &f0Files,
					   const mcDESPOT::OffResMode &f0fit,
					   const bool &use_finite);
Nifti parseInput(vector<shared_ptr<SignalFunctor>> &sigs,
				       vector<Nifti > &signalFiles,
				       vector<Nifti > &B1Files,
				       vector<Nifti > &f0Files,
					   const mcDESPOT::OffResMode &f0fit,
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
		
		if ((type == "SSFP") && (f0fit == mcDESPOT::OffResMode::Map)) {
			if (prompt)
				cout << "Enter path to f0 map: " << flush;
			getline(cin, path);
			f0Files.push_back(openAndCheck(path, templateFile, "f0"));
		} else {
			f0Files.push_back(Nifti());
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
	while ((c = getopt_long(argc, argv, "hvn123m:o:f:s:p:S:t:Fceijw", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case '1': components = Components::One; break;
			case '2': components = Components::Two; break;
			case '3': components = Components::Three; break;
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
			case 'f':
				switch (*optarg) {
					case '0' : f0fit = mcDESPOT::OffResMode::Map; break;
					case '1' : f0fit = mcDESPOT::OffResMode::SingleSymmetric; break;
					case '2' : f0fit = mcDESPOT::OffResMode::Single; break;
					default:
						cout << "Invalid Off Resonance Mode." << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 's': start_slice = atoi(optarg); break;
			case 'p': stop_slice = atoi(optarg); break;
			case 'S':
				switch (atoi(optarg)) {
					case 1 : scale = DESPOT2FM::Scaling::Global; break;
					case 2 : scale = DESPOT2FM::Scaling::NormToMean; break;
					default:
						cout << "Invalid scaling mode: " + to_string(atoi(optarg)) << endl;
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
			case 'F': use_finite = true; break;
			case 'c':
				cout << "Enter max number of contractions: " << flush; cin >> contract;
				cout << "Enter number of samples per contraction: " << flush; cin >> samples;
				cout << "Enter number of samples to retain: " << flush; cin >> retain;
				cout << "Enter fraction to expand region by: " << flush; cin >> expand;
				break;
			case 'e': extra = true; break;
			case 'E': early_finish = true; break;
			case 'i': voxI = atoi(optarg); break;
			case 'j': voxJ = atoi(optarg); break;
			case 'w': weighting = atof(optarg); break;
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
	vector<shared_ptr<SignalFunctor>> sigs;
	vector<Nifti > signalFiles, B1Files, f0Files;
	templateFile = parseInput(sigs, signalFiles, B1Files, f0Files, f0fit, use_finite);
	if ((maskData.size() > 0) && !(maskFile.matchesSpace(templateFile))) {
		cerr << "Mask file has different dimensions/transform to input data." << endl;
		exit(EXIT_FAILURE);
	}
	//**************************************************************************
	#pragma mark Allocate memory and set up boundaries.
	// Use if files are open to indicate default values should be used -
	// 0 for f0, 1 for B1
	//**************************************************************************
	size_t voxelsPerSlice = templateFile.voxelsPerSlice();
	size_t voxelsPerVolume = templateFile.voxelsPerVolume();
	
	vector<vector<double>> signalVolumes(signalFiles.size()),
	                 B1Volumes(signalFiles.size()),
				     f0Volumes(signalFiles.size());
	for (size_t i = 0; i < signalFiles.size(); i++) {
		signalVolumes[i].resize(voxelsPerSlice * sigs.at(i)->size());
		if (B1Files[i].isOpen()) B1Volumes[i].resize(voxelsPerSlice);
		if (f0Files[i].isOpen()) f0Volumes[i].resize(voxelsPerSlice);
	}
	
	// Build a Functor here so we can query number of parameters etc.
	cout << "Using " << mcDESPOT::to_string(components) << " component model." << endl;
	mcDESPOT mcd(components, sigs, tesla, f0fit, scale, use_finite);
	outPrefix = outPrefix + mcDESPOT::to_string(components) + "C_";
	ArrayXd weights(mcd.values()); weights.setConstant(1.0);
	ArrayXd threshes(mcd.inputs()); threshes.setConstant(0.05);
	if (early_finish)
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
		if (extra) {
			midpData.at(i).resize(voxelsPerSlice);
			midpFiles.at(i).open(outPrefix + mcd.names()[i] + "_mid.nii.gz", Nifti::Mode::Write);
			widthData.at(i).resize(voxelsPerSlice);
			widthFiles.at(i).open(outPrefix + mcd.names()[i] + "_width.nii.gz", Nifti::Mode::Write);
		}
	}
	if (extra)
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
	if (stop_slice > templateFile.dim(3))
		stop_slice = templateFile.dim(3);
	
	signal(SIGINT, int_handler);	// If we've got here there's actually allocated data to save
	
    time_t procStart = time(NULL);
	char theTime[512];
	strftime(theTime, 512, "%H:%M:%S", localtime(&procStart));
	cout << "Started processing at " << theTime << endl;
	for (size_t slice = start_slice; slice < stop_slice; slice++)
	{
		if (verbose) cout << "Reading data for slice " << slice << "..." << flush;
		atomic<int> voxCount{0};
		const size_t sliceOffset = slice * voxelsPerSlice;
		
		// Read data for slices
		for (size_t i = 0; i < signalFiles.size(); i++) {
			signalFiles[i].readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, signalVolumes[i]);
			if (B1Files[i].isOpen()) B1Files[i].readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, B1Volumes[i]);
			if (f0Files[i].isOpen()) f0Files[i].readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, f0Volumes[i]);
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
					
					if (f0fit == mcDESPOT::OffResMode::Map) {
						localf.m_f0 = f0Files[i].isOpen() ? f0Volumes[i][vox] : 0.;
					}
					localf.m_B1 = B1Files[i].isOpen() ? B1Volumes[i][vox] : 1.;
				}
				// Add the voxel number to the time to get a decent random seed
				size_t rSeed = time(NULL) + vox;
				RegionContraction<mcDESPOT> rc(localf, bounds, weights, threshes,
											 samples, retain, contract, expand, (voxI != -1));
				rc.optimise(params, rSeed);
				residuals = rc.residuals();
				if (extra) {
					c = rc.contractions();
					width = rc.width();
					midp = rc.midPoint();
				}
			}
			if (extra)
				contractData.at(vox) = c;
			for (size_t p = 0; p < paramsData.size(); p++) {
				paramsData.at(p).at(vox) = params[p];
				if (extra) {
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
			if (extra) {
				midpFiles.at(p).writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, midpData.at(p));
				widthFiles.at(p).writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, widthData.at(p));
			}
		}
		if (extra)
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

