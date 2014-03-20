/*
 *  mcsignal_main.cpp
 *
 *  Created by Tobias Wood on 12/11/2012.
 *  Copyright (c) 2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <string>
#include <iostream>
#include <getopt.h>
#include <exception>
#include <Eigen/Dense>
#include <Nifti/Nifti.h>
#include "ThreadPool.h"
#include "Model.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: mcsignal [options]\n\
\n\
Calculates multi-component DESPOT signals (mainly for debugging mcdespot).\n\
The program will prompt for input (unless --no-prompt specified)\n\
\n\
All times (TR) are in SECONDS. All angles are in degrees.\n\
\n\
Options:\n\
	--help, -h        : Print this message.\n\
	--verbose, -v     : Print extra information.\n\
	--mask, -m file   : Only calculate inside the mask.\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--B1, -b file     : B1 Map file (ratio)\n\
	--no-prompt, -n   : Don't print prompts for input.\n\
	--1, --2, --3     : Use 1, 2 or 3 component model (default 3).\n\
	--model, -M s     : Use simple model (default).\n\
	            f     : Use Finite Pulse Length correction.\n"
};

static auto components = Signal::Components::Three;
static auto modelType = ModelTypes::Simple;
static bool verbose = false, prompt = true;
static string outPrefix = "";
static struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"no-prompt", no_argument, 0, 'n'},
	{"1", no_argument, 0, '1'},
	{"2", no_argument, 0, '2'},
	{"3", no_argument, 0, '3'},
	{"model", no_argument, 0, 'M'},
	{"B1", required_argument, 0, 'b'},
	{0, 0, 0, 0}
};

//******************************************************************************
#pragma mark Read in all required files and data from cin
//******************************************************************************
void parseInput(shared_ptr<Model> &mdl);
void parseInput(shared_ptr<Model> &mdl) {
	string type;
	size_t nFlip, nPhase;
	if (prompt) cout << "Specify next signal type (SPGR/SSFP): " << flush;
	while (getline(cin, type) && (type != "END") && (type != "")) {
		if (type != "SPGR" && type != "SSFP") {
			cerr << "Unknown signal type: " << type << endl;
			exit(EXIT_FAILURE);
		}
		if (prompt) cout << "Number of Flip-Angles: " << flush;
		cin >> nFlip;
		if (type == "SPGR") {
			mdl->parseSPGR(nFlip, prompt);
		} else {
			if (prompt) cout << "Number of phase-cycles: " << flush;
			cin >> nPhase;
			mdl->parseSSFP(nFlip, nPhase, prompt);
		}
		// Print message ready for next loop
		if (prompt) cout << "Specify next image type (SPGR/SSFP, END to finish input): " << flush;
	}
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
	
	Nifti maskFile, B1File;
	Volume<int8_t> maskVol;
	Volume<float> B1Vol;
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvnm:o:b:123M:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskVol.readFrom(maskFile);
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'b':
				cout << "Reading B1 file: " << optarg << endl;
				B1File.open(optarg, Nifti::Mode::Read);
				B1Vol.readFrom(B1File);
				break;
			case '1': components = Signal::Components::One; break;
			case '2': components = Signal::Components::Two; break;
			case '3': components = Signal::Components::Three; break;
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
	#pragma mark  Set up model
	//**************************************************************************
	shared_ptr<Model> model;
	switch (modelType) {
		case ModelTypes::Simple: model = make_shared<SimpleModel>(components, Model::Scaling::NormToMean); break;
		case ModelTypes::Finite: model = make_shared<FiniteModel>(components, Model::Scaling::NormToMean); break;
	}
	parseInput(model);
	//**************************************************************************
	#pragma mark Read in parameter files
	//**************************************************************************
	// Build a Functor here so we can query number of parameters etc.
	cout << "Using " << Signal::to_string(components) << " component model." << endl;
	Series<float> paramsVols;
	Nifti saveFile;
	size_t numVoxels;
	if (prompt) cout << "Loading parameters." << endl;
	for (size_t i = 0; i < model->nParameters(); i++) {
		if (prompt) cout << "Enter path to " << model->names()[i] << " file: " << flush;
		string filename; cin >> filename;
		cout << "Reading " << filename << endl;
		Nifti input(filename, Nifti::Mode::Read);

		if (i == 0) {
			saveFile = Nifti(input, model->size());
			paramsVols = Series<float>(input.dims().head(3), model->nParameters());
			numVoxels = input.dims().head(3).prod();
		} else {
			if (!input.matchesSpace(saveFile)) {
				cout << "Mismatched input volumes" << endl;
				exit(EXIT_FAILURE);
			}
		}
		paramsVols.view(i).readFrom(input);
	}
	Series<float> signalVols(saveFile.dims().head(3), model->size());
	
	cout << "Started calculating." << endl;
	function<void (const size_t&)> calcVox = [&] (const size_t &v) {
		if ((maskFile.isOpen() == 0) || (maskVol[v])) {
			ArrayXd params = paramsVols.series(v).cast<double>();
			double B1 = B1File.isOpen() ? B1Vol[v] : 1.;
			signalVols.series(v) = model->signal(params, B1).abs().cast<float>();
		}
	};
	ThreadPool threads(1);
	threads.for_loop(calcVox, numVoxels);
	
	cout << "Finished calculating." << endl;
	saveFile.open(outPrefix + "mcsigout.nii.gz", Nifti::Mode::Write);
	signalVols.writeTo(saveFile);
	saveFile.close();
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

