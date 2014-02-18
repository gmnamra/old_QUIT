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
#include <fstream>
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
const string usage {
"Usage is: mcdespot [options]\n\
\n\
The program will prompt for input (unless --no-prompt specified)\n\
\n\
All times (TR) are in SECONDS. All angles are in degrees.\n\
\n\
Options:\n\
	--help, -h        : Print this message.\n\
	--verbose, -v     : Print writeResiduals information.\n\
	--no-prompt, -p   : Don't print prompts for input.\n\
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
	--model, -M s     : Use simple model (default).\n\
	            f     : Use Finite Pulse Length correction.\n\
	--contract, -c n  : Read contraction settings from stdin (Will prompt).\n"
};

static auto f0fit = OffRes::FitSym;
static auto components = Signal::Components::Three;
static auto modelType = ModelTypes::Simple;
static auto tesla = Model::FieldStrength::Three;
static auto scale = Model::Scaling::NormToMean;
static size_t start_slice = 0, stop_slice = numeric_limits<size_t>::max();
static int verbose = false, prompt = true, writeResiduals = false,
           early_finish = false, use_weights = false,
		   samples = 5000, retain = 50, contract = 10,
           voxI = -1, voxJ = -1;
static double expand = 0.;
static string outPrefix;
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"1", no_argument, 0, '1'},
	{"2", no_argument, 0, '2'},
	{"3", no_argument, 0, '3'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"f0", required_argument, 0, 'f'},
	{"start", required_argument, 0, 's'},
	{"stop", required_argument, 0, 'p'},
	{"scale", required_argument, 0, 'S'},
	{"tesla", required_argument, 0, 't'},
	{"model", no_argument, 0, 'M'},
	{"contract", no_argument, 0, 'c'},
	{0, 0, 0, 0}
};
//******************************************************************************
#pragma mark SIGTERM interrupt handler and Threads
//******************************************************************************
ThreadPool threads;
bool interrupt_received = false;
void int_handler(int); // Need the int to conform to handler definition but we don't use it
void int_handler(int) {
	cout << endl << "Stopping processing early." << endl;
	threads.stop();
	interrupt_received = true;
}

//******************************************************************************
#pragma mark Read in all required files and data from cin
//******************************************************************************
//Utility function
Nifti openAndCheck(const string &path, const Nifti &saved);
Nifti openAndCheck(const string &path, const Nifti &saved) {
	Nifti in(path, Nifti::Mode::Read);
	if (!(in.matchesSpace(saved))) {
		cerr << "Header for " << in.imagePath() << " does not match " << saved.imagePath() << endl;
		exit(EXIT_FAILURE);
	}
	if (verbose) cout << "Opened: " << in.imagePath() << endl;
	return in;
}

Nifti parseInput(shared_ptr<Model> &mdl, vector<Nifti> &signalFiles, Nifti &B1File, Nifti &f0File);
Nifti parseInput(shared_ptr<Model> &mdl, vector<Nifti> &signalFiles, Nifti &B1File, Nifti &f0File)
{
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
			if (verbose) cout << "Opened: " << signalFiles.back().imagePath() << endl;
			templateFile = Nifti(signalFiles.back(), 1); // Save header info for later
		} else {
			signalFiles.push_back(openAndCheck(path, templateFile));
		}
		#ifdef AGILENT
		Agilent::ProcPar pp;
		if (ReadPP(signalFiles.back(), pp)) {
			if (type == "SPGR") {
				mdl->procparseSPGR(pp);
			} else {
				mdl->procparseSSFP(pp);
			}
		} else
		#endif
		{
			if (type == "SPGR") {
				mdl->parseSPGR(signalFiles.back().dim(4), prompt);
			} else {
				size_t nPhases;
				if (prompt) cout << "Enter number of phase-cycling patterns: " << flush;
				cin >> nPhases;
				mdl->parseSSFP(signalFiles.back().dim(4) / nPhases, nPhases, prompt);
			}
		}
		// Print message ready for next loop
		if (prompt) cout << "Specify next image type (SPGR/SSFP, END to finish input): " << flush;
	}
	if (prompt) cout << "Enter B1 Map Path (Or NONE): " << flush;
	getline(cin, path);
	if ((path != "NONE") && (path != "")) {
		B1File = openAndCheck(path, templateFile);
	}
	if (f0fit == OffRes::Map) {
		if (prompt)
			cout << "Enter path to f0 map: " << flush;
		getline(cin, path);
		f0File = openAndCheck(path, templateFile);
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
	cout << version << endl << credit_me << endl;
	Eigen::initParallel();
	
	try { // To fix uncaught exceptions on Mac
	
	Nifti maskFile, templateFile;
	vector<double> maskData(0);
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvn123m:o:f:s:p:S:t:M:cei:j:w", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case '1': components = Signal::Components::One; break;
			case '2': components = Signal::Components::Two; break;
			case '3': components = Signal::Components::Three; break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskData.resize(maskFile.dims().head(3).prod());
				maskFile.readVolumes(0, 1, maskData);
				maskFile.close();
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'f':
				switch (*optarg) {
					case '0' : f0fit = OffRes::Map; break;
					case '1' : f0fit = OffRes::FitSym; break;
					case '2' : f0fit = OffRes::Fit; break;
					default:
						cout << "Invalid Off Resonance Mode." << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 's': start_slice = atoi(optarg); break;
			case 'p': stop_slice = atoi(optarg); break;
			case 'S':
				switch (atoi(optarg)) {
					case 1 : scale = Model::Scaling::None; break;
					case 2 : scale = Model::Scaling::NormToMean; break;
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
						cout << "Unknown boundaries type " << *optarg << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 'M':
				switch (*optarg) {
					case 's': modelType = ModelTypes::Simple; if (prompt) cout << "Simple model selected." << endl; break;
					case 'f': modelType = ModelTypes::Finite; if (prompt) cout << "Finite pulse correction selected." << endl; break;
					default:
						cout << "Unknown model type " << *optarg << endl;
						exit(EXIT_FAILURE);
						break;
				}
				break;
			case 'c':
				cout << "Enter max number of contractions: " << flush; cin >> contract;
				cout << "Enter number of samples per contraction: " << flush; cin >> samples;
				cout << "Enter number of samples to retain: " << flush; cin >> retain;
				cout << "Enter fraction to expand region by: " << flush; cin >> expand;
				{ string dummy; getline(cin, dummy); } // Eat newlines
				break;
			case 'e': writeResiduals = true; break;
			case 'i': voxI = atoi(optarg); break;
			case 'j': voxJ = atoi(optarg); break;
			case 'w': use_weights = true; break;
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
	shared_ptr<Model> model;
	// Build a Functor here so we can query number of parameters etc.
	cout << "Using " << Signal::to_string(components) << " component model." << endl;
	switch (modelType) {
		case ModelTypes::Simple : model = make_shared<SimpleModel>(components, scale);
		case ModelTypes::Finite : model = make_shared<FiniteModel>(components, scale);
	}
	vector<Nifti> signalFiles;
	Nifti B1File, f0File;
	templateFile = parseInput(model, signalFiles, B1File, f0File);
	if ((maskData.size() > 0) && !(maskFile.matchesSpace(templateFile))) {
		cerr << "Mask file has different dimensions/transform to input data." << endl;
		exit(EXIT_FAILURE);
	}
	//**************************************************************************
	#pragma mark Allocate memory and set up boundaries.
	// Use if files are open to indicate default values should be used -
	// 0 for f0, 1 for B1
	//**************************************************************************
	outPrefix = outPrefix + Signal::to_string(components) + "C_";
	ArrayXd threshes(model->nParameters()); threshes.setConstant(0.05);
	
	templateFile.setDim(4, 1);
	templateFile.setDatatype(Nifti::DataType::FLOAT32);
	templateFile.description = version;
	
	vector<Nifti> paramsFiles(model->nParameters(), templateFile);
	Nifti SoSFile(templateFile);
	Nifti residualsFile(templateFile, static_cast<int>(model->size()));

	size_t voxelsPerSlice = templateFile.dims().head(2).prod();
	size_t voxelsPerVolume = templateFile.dims().head(3).prod();
	vector<double> B1Slice(voxelsPerSlice);
	vector<double> f0Slice(voxelsPerSlice);
	vector<double> SoSSlice(voxelsPerSlice, 0.);
	vector<vector<double>> paramsSlice(model->nParameters());
	vector<vector<double>> residualsVolume(model->size());
	vector<vector<double>> sigSlices(signalFiles.size());
	for (int i = 0; i < model->nParameters(); i++) {
		paramsSlice.at(i).resize(voxelsPerSlice, 0.);
		paramsFiles.at(i).open(outPrefix + model->names()[i] + ".nii.gz", Nifti::Mode::Write);
	}
	SoSFile.open(outPrefix + "SoS.nii.gz", Nifti::Mode::Write);
	for (size_t i = 0; i < signalFiles.size(); i++) {
		sigSlices.at(i).resize(voxelsPerSlice * signalFiles.at(i).dim(4));
	}
	if (writeResiduals) {
		for (size_t i = 0; i < residualsVolume.size(); i ++)
			residualsVolume.at(i).resize(voxelsPerVolume, 0.);
	}
	
	ArrayXXd bounds = model->bounds(tesla);
	if (tesla == Model::FieldStrength::User) {
		if (prompt) cout << "Enter parameter pairs (low then high)" << endl;
		for (size_t i = 0; i < model->nParameters() - 1; i++) {
			if (prompt) cout << model->names()[i] << ": " << flush;
			cin >> bounds(i, 0) >> bounds(i, 1);
		}
	}
	if (f0fit == OffRes::FitSym) {
		bounds(model->nParameters() - 1, 0) = 0.;
	}
	ArrayXd weights(model->size()); weights.setOnes();
	if (verbose) {
		cout << "Low bounds: " << bounds.col(0).transpose() << endl;
		cout << "Hi bounds:  " << bounds.col(1).transpose() << endl;
	}
	ofstream boundsFile(outPrefix + "bounds.txt");
	for (int p = 0; p < model->nParameters(); p++) {
		boundsFile << model->names()[p] << "\t" << bounds.row(p) << endl;
	}
	boundsFile.close();
	
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
	for (size_t slice = start_slice; slice < stop_slice; slice++) {
		if (verbose) cout << "Reading data for slice " << slice << "..." << flush;
		atomic<int> voxCount{0};
		const size_t sliceOffset = slice * voxelsPerSlice;
		
		// Read data for slices
		Nifti::ArrayXs sliceStart(4), sliceSize(4);
		sliceStart << 0, 0, slice, 0;
		sliceSize << 0, 0, 1, 0; // Zeros will be replaced with dimension
		for (size_t i = 0; i < signalFiles.size(); i++) {
			signalFiles[i].readVoxels<double>(sliceStart, sliceSize, sigSlices[i]);
		}
		if (B1File.isOpen()) B1File.readVoxels<double>(sliceStart, sliceSize, B1Slice);
		if (f0File.isOpen()) f0File.readVoxels<double>(sliceStart, sliceSize, f0Slice);
		if (verbose) cout << "processing..." << endl;
		clock_t loopStart = clock();
		
		function<void (const size_t&)> processVox = [&] (const size_t &vox) {
			if ((maskData.size() == 0) || (maskData[sliceOffset + vox] > 0.)) {
				voxCount++;
				ArrayXd signal = model->loadSignals(sigSlices, voxelsPerSlice, vox);
				ArrayXXd localBounds = bounds;
				if (f0fit == OffRes::Map) {
					localBounds.row(model->nParameters() - 1).setConstant(f0Slice[vox]);
				}
				double B1 = B1File.isOpen() ? B1Slice[vox] : 1.;
				size_t rSeed = time(NULL) + vox; // Add the voxel number to the time to get a decent random seed
				DESPOTFunctor func(model, signal, B1, false);
				RegionContraction<DESPOTFunctor> rc(func, localBounds, weights, threshes,
											        samples, retain, contract, expand, (voxI != -1));
				ArrayXd params(model->nParameters());
				rc.optimise(params, rSeed);
				if (voxI != -1)
				if (verbose && (rc.status() == RegionContraction<DESPOTFunctor>::Status::ErrorResidual)) {
					cerr << "Thread ID: " << this_thread::get_id() << endl;
					cerr << "RC address: " << &rc << endl;
					cerr << "Slice: " << slice << "\tVoxel: " << vox << endl;
					cerr << "B1: " << B1 << endl;
					cerr << "nContract: " << rc.contractions() << endl;
					cerr << "Params: " << params.transpose() << endl;
					cerr << "Theory: " << model->signal(params, B1).transpose() << endl;
				}
				for (size_t p = 0; p < paramsSlice.size(); p++) {
					paramsSlice.at(p).at(vox) = params[p];
				}
				SoSSlice.at(vox) = rc.SoS();
				if (writeResiduals) {
					for (int i = 0; i < model->size(); i++) {
						residualsVolume.at(i).at(slice * voxelsPerSlice + vox) = rc.residuals()[i];
					}
				}
			}
		};
		if (voxI == -1)
			threads.for_loop(processVox, voxelsPerSlice);
		else {
			size_t voxInd = templateFile.dim(1) * voxJ + voxI;
			processVox(voxInd);
			exit(0);
		}
		for (size_t p = 0; p < paramsFiles.size(); p++) {
			paramsFiles.at(p).writeVoxels(sliceStart, sliceSize, paramsSlice.at(p));
		}
		SoSFile.writeVoxels(sliceStart, sliceSize, SoSSlice);
		if (verbose) {
			clock_t loopEnd = clock();
			if (voxCount > 0)
				cout << voxCount << " unmasked voxels, CPU time per voxel was "
				          << ((loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC)) << " s, ";
			cout << "finished." << endl;
		}
		if (interrupt_received)
			break;
	}
	time_t procEnd = time(NULL);
	strftime(theTime, 512, "%H:%M:%S", localtime(&procEnd));
	cout << "Finished processing at " << theTime << ". Run-time was " 
		 << difftime(procEnd, procStart) << " s." << endl;
	// Residuals can only be written here if we want them to go in a 4D gzipped file
	if (writeResiduals) {
		residualsFile.open(outPrefix + "residuals.nii.gz", Nifti::Mode::Write);
		for (size_t r = 0; r < residualsVolume.size(); r++)
			residualsFile.writeVolumes(r, 1, residualsVolume.at(r));
		residualsFile.close();
	}
	cout << "Finished writing data." << endl;
	
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

