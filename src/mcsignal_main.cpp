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
	--no-prompt, -n   : Don't print prompts for input.\n\
	--1, --2, --3     : Use 1, 2 or 3 component model (default 3).\n\
	--model, -M s     : Use simple model (default).\n\
	            e     : Use echo-time correction.\n\
	            f     : Use Finite Pulse Length correction.\n"
};

static auto components = Signal::Components::Three;
static auto modelType = ModelTypes::Simple;
static bool verbose = false, prompt = true;
static struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"1", no_argument, 0, '1'},
	{"2", no_argument, 0, '2'},
	{"3", no_argument, 0, '3'},
	{"model", no_argument, 0, 'M'},
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
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvn123M:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
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
	#pragma mark Allocate memory and set up boundaries.
	//**************************************************************************
	// Build a Functor here so we can query number of parameters etc.
	cout << "Using " << Signal::to_string(components) << " component model." << endl;
	VectorXd params(model->nParameters());
	if (prompt) cout << "Enter parameters." << endl;
	for (VectorXd::Index i = 0; i < params.rows(); i++) {
		if (prompt) cout << model->names()[i] << ": " << flush;
		cin >> params(i);
	}
	double B1;
	if (prompt) cout << "Enter B1: " << flush; cin >> B1;
	cout << model->signal(params, B1).transpose() << endl;
	
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

