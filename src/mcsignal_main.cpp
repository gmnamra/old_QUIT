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
#include <atomic>
#include <getopt.h>
#include <signal.h>
#include <time.h>
#include <Eigen/Dense>

#include "DESPOT_Functors.h"

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

static auto components = Components::Three;
static auto model = Model::Simple;
static int verbose = false, prompt = true;
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
void parseInput(vector<shared_ptr<SignalFunctor>> &sigs);
void parseInput(vector<shared_ptr<SignalFunctor>> &sigs) {
	string type;
	size_t nFlip;
	if (prompt) cout << "Specify next signal type (SPGR/SSFP): " << flush;
	while (getline(cin, type) && (type != "END") && (type != "")) {
		if (type != "SPGR" && type != "SSFP") {
			cerr << "Unknown signal type: " << type << endl;
			exit(EXIT_FAILURE);
		}
		if (prompt) cout << "Number of Flip-Angles: " << flush;
		cin >> nFlip;
		if (type == "SPGR")
			sigs.emplace_back(parseSPGR(components, model, nFlip, prompt, false));
		else
			sigs.emplace_back(parseSSFP(components, model, nFlip, prompt, false));
		// Print message ready for next loop
		if (prompt) cout << "Specify next image type (SPGR/SSFP, END to finish input): " << flush;
	}
	if (sigs.size() == 0) {
		cerr << "No signals specified." << endl;
		exit(EXIT_FAILURE);
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
	
	Nifti maskFile, templateFile;
	vector<double> maskData(0);
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvn123M:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case '1': components = Components::One; break;
			case '2': components = Components::Two; break;
			case '3': components = Components::Three; break;
			case 'M':
				switch (*optarg) {
					case 's': model = Model::Simple; if (prompt) cout << "Simple model selected." << endl; break;
					case 'e': model = Model::Echo; if (prompt) cout << "TE correction selected." << endl; break;
					case 'f': model = Model::Finite; if (prompt) cout << "Finite pulse correction selected." << endl; break;
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
	#pragma mark  Set up SPGR & SSFP lists
	//**************************************************************************
	vector<shared_ptr<SignalFunctor>> sigs;
	parseInput(sigs);
	//**************************************************************************
	#pragma mark Allocate memory and set up boundaries.
	//**************************************************************************
	// Build a Functor here so we can query number of parameters etc.
	cout << "Using " << mcDESPOT::to_string(components) << " component model." << endl;
	mcDESPOT mcd(components, sigs, mcDESPOT::FieldStrength::Unknown, mcDESPOT::OffRes::Fit, mcDESPOT::Scaling::NormToMean, model == Model::Finite, true);
	VectorXd params(mcd.inputs());
	if (prompt) cout << "Enter parameters." << endl;
	for (size_t i = 0; i < mcd.nP(); i++) {
		if (prompt) cout << mcd.names()[i] << ": " << flush;
		cin >> params(i);
	}
	mcd.theory(params);
	
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

