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
#include "QUIT/QUIT.h"
#include "Model.h"

using namespace std;
using namespace Eigen;
using namespace QUIT;

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
static bool verbose = false, prompt = true, finiteModel = false;
static string outPrefix = "signal_";
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
void parseInput(Model &mdl);
void parseInput(Model &mdl) {
	string type;
	if (prompt) cout << "Specify next signal type (SPGR/SSFP): " << flush;
	while (getline(cin, type) && (type != "END") && (type != "")) {
		if (type != "SPGR" && type != "SSFP") {
			cerr << "Unknown signal type: " << type << endl;
			exit(EXIT_FAILURE);
		}
		if ((type == "SPGR") && !finiteModel) {
			mdl.addSignal(SignalType::SPGR, prompt);
		} else if ((type == "SPGR" && finiteModel)) {
			mdl.addSignal(SignalType::SPGR_Finite, prompt);
		} else if ((type == "SSFP" && !finiteModel)) {
			mdl.addSignal(SignalType::SSFP, prompt);
		} else if ((type == "SSFP" && finiteModel)) {
			mdl.addSignal(SignalType::SSFP_Finite, prompt);
		}
		// Print message ready for next loop
		string temp; getline(cin, temp); // Just to eat the newline
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
	
	Nifti::File maskFile, B1File;
	MultiArray<int8_t, 3> maskVol;
	MultiArray<float, 3> B1Vol;
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvnm:o:b:123M:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskVol.resize(maskFile.dims());
				maskFile.readVolumes(maskVol.begin(), maskVol.end(), 0, 1);
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'b':
				cout << "Reading B1 file: " << optarg << endl;
				B1File.open(optarg, Nifti::Mode::Read);
				B1Vol.resize(B1File.dims());
				B1File.readVolumes(B1Vol.begin(), B1Vol.end(), 0, 1);
				break;
			case '1': components = Signal::Components::One; break;
			case '2': components = Signal::Components::Two; break;
			case '3': components = Signal::Components::Three; break;
			case 'M':
				switch (*optarg) {
					case 's': finiteModel = false; if (prompt) cout << "Simple model selected." << endl; break;
					case 'f': finiteModel = true; if (prompt) cout << "Finite pulse correction selected." << endl; break;
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
	Model model(components, Model::Scaling::None);
	parseInput(model);
	cout << model << endl;
	//**************************************************************************
	#pragma mark Read in parameter files
	//**************************************************************************
	// Build a Functor here so we can query number of parameters etc.
	cout << "Using " << Signal::to_string(components) << " component model." << endl;
	MultiArray<float, 4> paramsVols;
	Nifti::Header templateHdr;
	if (prompt) cout << "Loading parameters." << endl;
	for (size_t i = 0; i < model.nParameters(); i++) {
		if (prompt) cout << "Enter path to " << model.names()[i] << " file: " << flush;
		string filename; cin >> filename;
		cout << "Opening " << filename << endl;
		Nifti::File input(filename, Nifti::Mode::Read);

		if (i == 0) {
			templateHdr = input.header();
			paramsVols = MultiArray<float, 4>(input.header().fulldims().head(3), model.nParameters());
		} else {
			if (!input.header().matchesSpace(templateHdr)) {
				cout << "Mismatched input volumes" << endl;
				exit(EXIT_FAILURE);
			}
		}
		auto inVol = paramsVols.slice<3>({0,0,0,i},{-1,-1,-1,0});
		cout << "Reading data." << endl;
		input.readVolumes(inVol.begin(), inVol.end(), 0, 1);
	}
	auto d = paramsVols.dims();
	MultiArray<complex<float>, 4> signalVols(d.head(3), model.size());
	cout << "Calculating..." << endl;
	function<void (const size_t&)> calcVox = [&] (const size_t &k) {
		for (size_t j = 0; j < d[1]; j++) {
			for (size_t i = 0; i < d[0]; i++) {
				if ((maskFile.isOpen() == 0) || (maskVol[{i,j,k}])) {
					ArrayXd params = paramsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray().cast<double>();
					double B1 = B1File.isOpen() ? B1Vol[{i,j,k}] : 1.;
					signalVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = model.signal(params, B1).cast<complex<float>>();
				}
			}
		}
	};
	ThreadPool threads;
	threads.for_loop(calcVox, d[2]);
	
	cout << "Finished calculating." << endl;
	cout << "Saving data." << endl;
	templateHdr.setDatatype(Nifti::DataType::COMPLEX64);
	size_t startVol = 0;
	for (size_t i = 0; i < model.m_signals.size(); i++) {
		size_t thisSize = model.m_signals[i]->size();
		templateHdr.setDim(4, thisSize);
		Nifti::File saveFile(templateHdr, outPrefix + to_string(i) + ".nii.gz");
		auto thisSignal = signalVols.slice<4>({0,0,0,startVol},{-1,-1,-1,thisSize});
		saveFile.writeVolumes(thisSignal.begin(), thisSignal.end());
		saveFile.close();
	}
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

