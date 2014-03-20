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
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--mask, -m file   : Mask input with specified file\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--f0, -f SYM     : Fit symmetric f0 map (default)\n\
	         ASYM    : Fit asymmetric f0 map\n\
	         file    : Use f0 Map file (in Hertz)\n\
	--B1, -b file     : B1 Map file (ratio)\n\
	--start, -s n     : Only start processing at slice n.\n\
	--stop, -p n      : Finish at slice n-1\n\
	--scale, -S 0     : Normalise signals to mean (default)\n\
	            1     : Fit a scaling factor/proton density\n\
	--tesla, -t 3     : Boundaries suitable for 3T (default)\n\
	            7     : Boundaries suitable for 7T \n\
	            u     : User specified boundaries from stdin\n\
	--model, -M s     : Use simple model (default)\n\
	            f     : Use Finite Pulse Length correction\n\
	--complex, -x    : Fit to complex data\n\
	--contract, -c n  : Read contraction settings from stdin (Will prompt)\n\
	--resid, -r       : Write out per-flip angle residuals\n\
	--no-prompt, -n   : Don't print prompts for input\n\
	--1, --2, --3     : Use 1, 2 or 3 component model (default 3)\n"
};

static auto components = Signal::Components::Three;
static auto modelType = ModelTypes::Simple;
static auto scale = Model::Scaling::NormToMean;
static auto tesla = Model::FieldStrength::Three;
static auto f0fit = OffRes::FitSym;
static size_t start_slice = 0, stop_slice = numeric_limits<size_t>::max();
static int verbose = false, prompt = true, writeResiduals = false, fitComplex = false,
           samples = 5000, retain = 50, contract = 10,
           voxI = -1, voxJ = -1;
static double expand = 0.;
static string outPrefix;
static struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
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
	{"no-prompt", no_argument, 0, 'n'},
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
	return in;
}

Nifti parseInput(shared_ptr<Model> &mdl, vector<Series<complex<float>>> &signalVols);
Nifti parseInput(shared_ptr<Model> &mdl, vector<Series<complex<float>>> &signalVols)
{
	Nifti templateFile, inFile;
	string type, path;
	if (prompt) cout << "Specify next image type (SPGR/SSFP): " << flush;
	while (getline(cin, type) && (type != "END") && (type != "")) {
		if (type != "SPGR" && type != "SSFP") {
			cerr << "Unknown signal type: " << type << endl;
			exit(EXIT_FAILURE);
		}
		if (prompt) cout << "Enter image path: " << flush;
		getline(cin, path);
		if (signalVols.size() == 0) {
			inFile.open(path, Nifti::Mode::Read);
			templateFile = Nifti(inFile, 1); // Save header info for later
		} else {
			inFile = openAndCheck(path, templateFile);
		}
		if (verbose) cout << "Opened: " << inFile.imagePath() << endl;
		#ifdef AGILENT
		Agilent::ProcPar pp;
		if (ReadPP(inFile, pp)) {
			if (type == "SPGR") {
				mdl->procparseSPGR(pp);
			} else {
				mdl->procparseSSFP(pp);
			}
		} else
		#endif
		{
			if (type == "SPGR") {
				mdl->parseSPGR(inFile.dim(4), prompt);
			} else {
				size_t nPhases;
				if (prompt) cout << "Enter number of phase-cycling patterns: " << flush;
				cin >> nPhases;
				mdl->parseSSFP(inFile.dim(4) / nPhases, nPhases, prompt);
			}
		}
		signalVols.emplace_back(Series<complex<float>>(inFile));
		inFile.close();
		// Print message ready for next loop
		if (prompt) cout << "Specify next image type (SPGR/SSFP, END to finish input): " << flush;
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
	
	Nifti maskFile, f0File, B1File, templateFile;
	Volume<int8_t> maskVol;
	Volume<float> f0Vol, B1Vol;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvm:o:f:b:s:p:S:t:M:xcrn123i:j:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case '1': components = Signal::Components::One; break;
			case '2': components = Signal::Components::Two; break;
			case '3': components = Signal::Components::Three; break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskVol.readFrom(maskFile);
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
					f0Vol.readFrom(f0File);
					f0fit = OffRes::Map;
				}
				break;
			case 'b':
				cout << "Reading B1 file: " << optarg << endl;
				B1File.open(optarg, Nifti::Mode::Read);
				B1Vol.readFrom(B1File);
				break;
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
			case 'r': writeResiduals = true; break;
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
	shared_ptr<Model> model;
	// Build a Functor here so we can query number of parameters etc.
	cout << "Using " << Signal::to_string(components) << " component model." << endl;
	switch (modelType) {
		case ModelTypes::Simple : model = make_shared<SimpleModel>(components, scale); break;
		case ModelTypes::Finite : model = make_shared<FiniteModel>(components, scale); break;
	}
	vector<Series<complex<float>>> signalVols;
	templateFile = parseInput(model, signalVols);
	if ((maskFile.isOpen() && !templateFile.matchesSpace(maskFile)) ||
		(f0File.isOpen() && !templateFile.matchesSpace(f0File)) ||
		(B1File.isOpen() && !templateFile.matchesSpace(B1File))){
		cerr << "Dimensions/transforms do not match in input files." << endl;
		exit(EXIT_FAILURE);
	}
	//**************************************************************************
	#pragma mark Allocate memory and set up boundaries.
	//**************************************************************************
	Series<float> paramsVols(templateFile.dims().head(3), model->nParameters());
	Series<float> residualVols(templateFile.dims().head(3), model->size());;
	Volume<float> SoSVol(templateFile.dims().head(3));
	
	ArrayXd threshes(model->nParameters()); threshes.setConstant(0.05);
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
		cout << *model;
		cout << "Bounds:" << endl <<  bounds.transpose() << endl;
		ofstream boundsFile(outPrefix + "bounds.txt");
		for (size_t p = 0; p < model->nParameters(); p++) {
			boundsFile << model->names()[p] << "\t" << bounds.row(p) << endl;
		}
		boundsFile.close();
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
	for (size_t k = start_slice; k < stop_slice; k++) {
		if (verbose) cout << "Processing slice " << k << "..." << flush;
		atomic<int> voxCount{0};
		clock_t loopStart = clock();
		
		for (size_t j = voxJ; j < templateFile.dim(2); j++) {
			function<void (const size_t&)> processVox = [&] (const size_t &i) {
				const Volume<float>::IndexArray vox{i, j, k};
				if (maskFile.isOpen() || maskVol[vox]) {
					voxCount++;
					ArrayXcd signal = model->loadSignals(signalVols, vox);
					ArrayXXd localBounds = bounds;
					if (f0fit == OffRes::Map) {
						localBounds.row(model->nParameters() - 1).setConstant(f0Vol[vox]);
					}
					double B1 = B1File.isOpen() ? B1Vol[vox] : 1.;
					DESPOTFunctor func(model, signal, B1, fitComplex, false);
					RegionContraction<DESPOTFunctor> rc(func, localBounds, weights, threshes,
														samples, retain, contract, expand, (voxI != -1));
					ArrayXd params(model->nParameters());
					rc.optimise(params, time(NULL) + i); // Add the voxel number to the time to get a decent random seed
					paramsVols.series(vox) = params.cast<float>();
					residualVols.series(vox) = rc.residuals().cast<float>();
					SoSVol[vox] = static_cast<float>(rc.SoS());
				}
			};
			if (voxI == 0)
				threads.for_loop(processVox, templateFile.dim(1));
			else {
				processVox(voxI);
				exit(0);
			}
		}
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
	outPrefix = outPrefix + Signal::to_string(components) + "C_";
	templateFile.setDim(4, 1);
	templateFile.setDatatype(Nifti::DataType::FLOAT32);
	templateFile.description = version;
	for (size_t p = 1; p < model->nParameters(); p++) { // Skip PD for now
		templateFile.open(outPrefix + model->names().at(p) + ".nii.gz", Nifti::Mode::Write);
		paramsVols.view(p).writeTo(templateFile);
		templateFile.close();
	}
	templateFile.open(outPrefix + "SoS.nii.gz", Nifti::Mode::Write);
	SoSVol.writeTo(templateFile);
	templateFile.close();
	if (writeResiduals) {
		templateFile.setDim(4, static_cast<int>(model->size()));
		templateFile.open(outPrefix + "residuals.nii.gz", Nifti::Mode::Write);
		residualVols.writeTo(templateFile);
		templateFile.close();
	}
	cout << "Finished writing data." << endl;
	
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

