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
#include "DESPOT.h"
#include "DESPOT_Functors.h"
#include "RegionContraction.h"
#include "ThreadPool.h"

#ifdef AGILENT
	#include "procpar.h"
#endif

using namespace std;
using namespace Eigen;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: despot2-fm [options] T1_map ssfp_files\n\
\
Options:\n\
	--help, -h          : Print this message.\n\
	--verbose, -v       : Print slice processing times.\n\
	--mask, -m file     : Mask input with specified file.\n\
	--out, -o path      : Add a prefix to the output filenames.\n\
	--f0, -f file/ASYM  : f0 Map file in Hertz, or asymmetric fitting.\n\
	--B1, -b file       : B1 Map file.\n\
	--start, -s N       : Start processing from slice N.\n\
	--stop, -p  N       : Stop processing at slice N.\n\
	--scale, -S 0       : Normalise signals to mean (default).\n\
	            1       : Fit a scaling factor/proton density.\n\
	--tesla, -t 3       : Use boundaries suitable for 3T (default)\n\
	            7       : Boundaries suitable for 7T\n\
	            u       : User specified boundaries from stdin.\n\
	--model, -M s       : Use simple model (default).\n\
	            e       : Use echo-time correction.\n\
				f       : Use finite pulse length correction.\n\
	--contract, -c n    : Read contraction settings from stdin (Will prompt).\n"
};

static auto scale = Model::Scaling::NormToMean;
static auto tesla = Model::FieldStrength::Three;
static auto f0Fit = OffRes::FitSym;
static int verbose = false, extra = false,
		   samples = 2000, retain = 20, contract = 10,
           voxI = -1, voxJ = -1;
static auto modelType = ModelTypes::Simple;
static size_t start_slice = 0, stop_slice = numeric_limits<size_t>::max();
static double expand = 0., weighting = 1.0;
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
	{"contract", no_argument, 0, 'c'},
	{0, 0, 0, 0}
};

//******************************************************************************
#pragma mark SIGTERM interrupt handler and Threads
//******************************************************************************
ThreadPool threads;
bool interrupt_received = false;
void int_handler(int sig);
void int_handler(int) {
	cout << endl << "Stopping processing early." << endl;
	threads.stop();
	interrupt_received = true;
}

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	cout << version << endl << credit_me << endl;
	Eigen::initParallel();
	
	try { // To fix uncaught exceptions on Mac
	
	Nifti maskFile, f0File, B1File;
	vector<double> maskData, f0Data, B1Data;
	string procPath;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvm:o:f:b:s:p:S:t:M:cei:j:w", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskData.resize(maskFile.dims().head(3).prod());
				maskFile.readVolumes(0, 1, maskData);
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'f':
				if (string(optarg) == "ASYM") {
					f0Fit = OffRes::Fit;
				} else {
					cout << "Reading f0 file: " << optarg << endl;
					f0File.open(optarg, Nifti::Mode::Read);
					f0Data.resize(f0File.dims().head(3).prod());
					f0File.readVolumes(0, 1, f0Data);
					f0Fit = OffRes::Map;
				}
				break;
			case 'b':
				cout << "Reading B1 file: " << optarg << endl;
				B1File.open(optarg, Nifti::Mode::Read);
				B1Data.resize(B1File.dims().head(3).prod());
				B1File.readVolumes(0, 1, B1Data);
				break;
			case 's': start_slice = atoi(optarg); break;
			case 'p': stop_slice = atoi(optarg); break;
			case 'S':
				switch (atoi(optarg)) {
					case 0 : scale = Model::Scaling::None; break;
					case 1 : scale = Model::Scaling::NormToMean; break;
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
						cout << "Unknown boundaries type " << optarg << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 'M':
				switch (*optarg) {
					case 's': modelType = ModelTypes::Simple; cout << "Simple model selected." << endl; break;
					case 'f': modelType = ModelTypes::Finite; cout << "Finite pulse correction selected." << endl; break;
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
				break;
			case 'e': extra = true; break;
			case 'i': voxI = atoi(optarg); break;
			case 'j': voxJ = atoi(optarg); break;
			case 'w': weighting = atof(optarg); break;
			case '?': // getopt will print an error message
			case 'h':
			default:
				cout << usage << endl;
				exit(EXIT_SUCCESS);
				break;
		}
	}
	if ((argc - optind) < 2) {
		cout << "Wrong number of arguments. Need at least a T1 map and 1 SSFP file." << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Reading T1 Map from: " << argv[optind] << endl;
	Nifti inFile(argv[optind++], Nifti::Mode::Read);
	size_t voxelsPerSlice = inFile.dims().head(2).prod();
	size_t voxelsPerVolume = inFile.dims().head(3).prod();
	vector<double> T1Data(voxelsPerVolume);
	inFile.readVolumes(0, 1, T1Data);
	inFile.close();
	if ((maskFile.isOpen() && !inFile.matchesSpace(maskFile)) ||
	    (f0File.isOpen() && !inFile.matchesSpace(f0File)) ||
		(B1File.isOpen() && !inFile.matchesSpace(B1File))){
		cerr << "Dimensions/transforms do not match in input files." << endl;
		exit(EXIT_FAILURE);
	}
	Nifti templateFile(inFile, 1); // Save header data to write out results
	//**************************************************************************
	// Gather SSFP Data
	//**************************************************************************
	size_t nFiles = argc - optind;
	vector<vector<double>> ssfpData(nFiles);
	shared_ptr<Model> model;
	switch (modelType) {
		case ModelTypes::Simple: model = make_shared<SimpleModel>(Signal::Components::One, scale); break;
		case ModelTypes::Finite: model = make_shared<FiniteModel>(Signal::Components::One, scale); break;
	}
	VectorXd inFlip;
	for (size_t p = 0; p < nFiles; p++) {
		cout << "Reading SSFP header from " << argv[optind] << endl;
		inFile.open(argv[optind], Nifti::Mode::Read);
		if (p == 0)
			templateFile = Nifti(inFile, 1);
		if (!inFile.matchesSpace(templateFile)) {
			cerr << "Input file dimensions and/or transforms do not match." << endl;
			exit(EXIT_FAILURE);
		}
		#ifdef AGILENT
		Agilent::ProcPar pp;
		if (ReadPP(inFile, pp)) {
			model->procparseSPGR(pp);
		} else
		#endif
		{
			size_t nPhases;
			cout << "Enter number of phase-cycling patterns: " << flush;
			cin >> nPhases;
			model->parseSSFP(inFile.dim(4) / nPhases, nPhases, true);
		}
		cout << "Reading data..." << endl;
		ssfpData.at(p).resize(inFile.dims().head(4).prod());
		inFile.readVolumes(0, inFile.dim(4), ssfpData.at(p));
		inFile.close();
		optind++;
	}
	size_t totalSize = model->size();
	if (optind != argc) {
		cerr << "Unprocessed arguments supplied.\n" << usage;
		exit(EXIT_FAILURE);
	}
	
	ArrayXd thresh(model->nParameters()); thresh.setConstant(0.05);
	ArrayXXd bounds = model->bounds(tesla);
	if (tesla == Model::FieldStrength::User) {
		cout << "Enter parameter pairs (low then high)" << endl;
		for (size_t i = 0; i < model->nParameters() - 1; i++) {
			cout << model->names()[i] << ": " << flush;
			cin >> bounds(i, 0) >> bounds(i, 1);
		}
	}
	if (f0Fit == OffRes::FitSym) {
		bounds(model->nParameters() - 1, 0) = 0.;
	}
	ArrayXd weights(totalSize); weights.setOnes();
	
	if (verbose) {
		cout << "Low bounds: " << bounds.col(0).transpose() << endl
		     << "Hi bounds:  " << bounds.col(1).transpose() << endl;
	}
	//**************************************************************************
	// Set up results data
	//**************************************************************************
	vector<vector<double>> paramsData(model->nParameters());
	for (auto &p : paramsData)
		p.resize(voxelsPerVolume, 0.);
	vector<vector<double>> residuals(totalSize);
	for (auto &r : residuals)
		r.resize(voxelsPerVolume, 0.);
	vector<size_t> contractData;
	vector<vector<double>> midpData(totalSize);
	vector<vector<double>> widthData(totalSize);
	if (extra) {
		contractData.resize(voxelsPerVolume);
		for (auto &m : midpData)
			m.resize(voxelsPerVolume, 0.);
		for (auto &r : widthData)
			r.resize(voxelsPerVolume, 0.);
	}
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	if (stop_slice > templateFile.dim(3))
		stop_slice = templateFile.dim(3);
    time_t procStart = time(NULL);
	char theTime[512];
	strftime(theTime, 512, "%H:%M:%S", localtime(&procStart));
	cout << "Started processing at " << theTime << endl;
	signal(SIGINT, int_handler);	// If we've got here there's actually allocated data to save
	for (size_t slice = start_slice; slice < stop_slice; slice++) {
		// Read in data
		if (verbose)
			cout << "Starting slice " << slice << "..." << flush;
		
		atomic<int> voxCount{0};
		const size_t sliceOffset = slice * voxelsPerSlice;
		clock_t loopStart = clock();
		function<void (const size_t&)> processVox = [&] (const size_t &vox) {
			// Set up parameters and constants
			ArrayXd params(model->nParameters()); params.setZero();
			ArrayXd resid(model->size()); resid.setZero();
			// extra stuff
			size_t c = 0;
			ArrayXd width(model->size()); width.setZero();
			ArrayXd midp(model->size()); midp.setZero();
			if (!maskFile.isOpen() || ((maskData[sliceOffset + vox] > 0.) && (T1Data[sliceOffset + vox] > 0.)))
			{	// -ve T1 is non-sensical, no point fitting
				voxCount++;
				ArrayXd signal = model->loadSignals(ssfpData, voxelsPerVolume, vox);
				ArrayXXd localBounds = bounds;
				localBounds.row(0).setConstant(T1Data[sliceOffset + vox]);
				if (f0Fit == OffRes::Map) {
					localBounds.row(2).setConstant(f0Data[vox]);
				}
				double B1 = B1File.isOpen() ? B1Data[vox] : 1.;
				size_t rSeed = time(NULL) + vox; // Add the voxel number to the time to get a decent random seed
				DESPOTFunctor func(model, signal, B1, false);
				RegionContraction<DESPOTFunctor> rc(func, localBounds, weights, thresh,
											        samples, retain, contract, expand, (voxI != -1));
				rc.optimise(params);
				resid = rc.residuals();
				if (extra) {
					c = rc.contractions();
					width = rc.width();
					midp = rc.midPoint();
				}
			}
			for (ArrayXd::Index p = 0; p < params.size(); p++) {
				paramsData.at(p).at(sliceOffset + vox) = params(p);
			}
			for (ArrayXd::Index i = 0; i < resid.size(); i++) {
				residuals.at(i).at(sliceOffset + vox) = resid(i);
			}
			if (extra) {
				contractData.at(sliceOffset + vox) = c;
				for (size_t i = 0; i < widthData.size(); i++) {
					widthData.at(i).at(sliceOffset + vox) = width(i);
					midpData.at(i).at(sliceOffset + vox) = midp(i);
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
	
	outPrefix = outPrefix + "FM_";
	templateFile.description = version;
	for (int p = 1; p < model->nParameters(); p++) { // Skip T1
		templateFile.open(outPrefix + model->names().at(p) + ".nii.gz", Nifti::Mode::Write);
		templateFile.writeVolumes(0, 1, paramsData.at(p));
		templateFile.close();
	}
	templateFile.setDim(4, static_cast<int>(residuals.size()));
	templateFile.open(outPrefix + "residuals.nii.gz", Nifti::Mode::Write);
	for (size_t i = 0; i < residuals.size(); i++)
		templateFile.writeVolumes(i, 1, residuals[i]);
	templateFile.close();
	if (extra) {
		templateFile.setDim(4, 1);
		templateFile.setDatatype(Nifti::DataType::INT16);
		templateFile.open(outPrefix + "n_contract.nii.gz", Nifti::Mode::Write);
		templateFile.writeVolumes(0, 1, contractData);
		templateFile.close();
		templateFile.setDatatype(Nifti::DataType::FLOAT32);
		for (int p = 0; p < model->nParameters(); p++) {
			templateFile.open(outPrefix + model->names().at(p) + "_width.nii.gz", Nifti::Mode::Write);
			templateFile.writeVolumes(0, 1, widthData.at(p));
			templateFile.close();
			templateFile.open(outPrefix + model->names().at(p) + "_mid.nii.gz", Nifti::Mode::Write);
			templateFile.writeVolumes(0, 1, midpData.at(p));
			templateFile.close();
		}
	}
	
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
