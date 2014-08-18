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
	--noise, -N val   : Add complex noise with std=val.\n\
	--1, --2, --3     : Use 1, 2 or 3 component sequences (default 3).\n\
	--sequences, -M s     : Use simple sequences (default).\n\
	            f     : Use Finite Pulse Length correction.\n"
};

static auto components = Components::Three;
static bool verbose = false, prompt = true, finitesequences = false;
static string outPrefix = "";
static double sigma = 0.;
static struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"no-prompt", no_argument, 0, 'n'},
	{"noise", required_argument, 0, 'N'},
	{"1", no_argument, 0, '1'},
	{"2", no_argument, 0, '2'},
	{"3", no_argument, 0, '3'},
	{"sequences", no_argument, 0, 'M'},
	{"B1", required_argument, 0, 'b'},
	{0, 0, 0, 0}
};

//******************************************************************************
#pragma mark Read in all required files and data from cin
//******************************************************************************
void parseInput(Sequences &cs, vector<string> &names);
void parseInput(Sequences &cs, vector<string> &names) {
	string type;
	if (prompt) cout << "Specify next signal type (SPGR/SSFP): " << flush;
	while (getline(cin, type) && (type != "END") && (type != "")) {
		if (type != "SPGR" && type != "SSFP") {
			throw(std::runtime_error("Unknown signal type: " + type));
		}
		if (prompt) cout << "Enter output filename: " << flush;
		string filename;
		getline(cin, filename);
		names.push_back(filename);
		if ((type == "SPGR") && !finitesequences) {
			cs.addSequence(SequenceType::SPGR, prompt);
		} else if ((type == "SPGR" && finitesequences)) {
			cs.addSequence(SequenceType::SPGR_Finite, prompt);
		} else if ((type == "SSFP" && !finitesequences)) {
			cs.addSequence(SequenceType::SSFP, prompt);
		} else if ((type == "SSFP" && finitesequences)) {
			cs.addSequence(SequenceType::SSFP_Finite, prompt);
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
	while ((c = getopt_long(argc, argv, "hvnN:m:o:b:123M:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'N': sigma = atof(optarg); break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskVol.resize(maskFile.matrix());
				maskFile.readVolumes(maskVol.begin(), maskVol.end(), 0, 1);
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'b':
				cout << "Reading B1 file: " << optarg << endl;
				B1File.open(optarg, Nifti::Mode::Read);
				B1Vol.resize(B1File.matrix());
				B1File.readVolumes(B1Vol.begin(), B1Vol.end(), 0, 1);
				break;
			case '1': components = Components::One; break;
			case '2': components = Components::Two; break;
			case '3': components = Components::Three; break;
			case 'M':
				switch (*optarg) {
					case 's': finitesequences = false; if (prompt) cout << "Simple sequences selected." << endl; break;
					case 'f': finitesequences = true; if (prompt) cout << "Finite pulse correction selected." << endl; break;
					default:
						cout << "Unknown sequences type " << *optarg << endl;
						return EXIT_FAILURE;
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
	#pragma mark  Set up sequences
	//**************************************************************************
	Sequences sequences(components, Scale::None);
	vector<string> filenames;
	parseInput(sequences, filenames);
	cout << sequences << endl;
	//**************************************************************************
	#pragma mark Read in parameter files
	//**************************************************************************
	// Build a Functor here so we can query number of parameters etc.
	cout << "Using " << to_string(components) << " component sequences." << endl;
	MultiArray<float, 4> paramsVols;
	Nifti::Header templateHdr;
	if (prompt) cout << "Loading parameters." << endl;
	for (size_t i = 0; i < sequences.nParameters(); i++) {
		if (prompt) cout << "Enter path to " << sequences.names()[i] << " file: " << flush;
		string filename; cin >> filename;
		cout << "Opening " << filename << endl;
		Nifti::File input(filename);

		if (i == 0) {
			templateHdr = input.header();
			paramsVols = MultiArray<float, 4>(input.matrix(), sequences.nParameters());
		} else {
			if (!input.header().matchesSpace(templateHdr)) {
				cout << "Mismatched input volumes" << endl;
				return EXIT_FAILURE;
			}
		}
		auto inVol = paramsVols.slice<3>({0,0,0,i},{-1,-1,-1,0});
		cout << "Reading data." << endl;
		input.readVolumes(inVol.begin(), inVol.end(), 0, 1);
	}
	const auto d = paramsVols.dims();
	vector<MultiArray<complex<float>, 4>> signalVols(sequences.count()); //d.head(3), sequences.combinedSize());
	for (size_t s = 0; s < sequences.count(); s++) {
		signalVols[s] = MultiArray<complex<float>, 4>(d.head(3), sequences.sequence(s)->size());
	}
	cout << "Calculating..." << endl;
	function<void (const size_t&)> calcVox = [&] (const size_t &k) {
		for (size_t j = 0; j < d[1]; j++) {
			for (size_t i = 0; i < d[0]; i++) {
				if (!maskFile || (maskVol[{i,j,k}])) {
					ArrayXd params = paramsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray().cast<double>();
					double B1 = B1File ? B1Vol[{i,j,k}] : 1.;
					for (size_t s = 0; s < sequences.count(); s++) {
						ArrayXcd signal = sequences.signal(s, params, B1);
						ArrayXcd noise(sequences.combinedSize());
						noise.real() = (ArrayXd::Ones(sequences.combinedSize()) * sigma).unaryExpr(function<double(double)>(randNorm<double>));
						noise.imag() = (ArrayXd::Ones(sequences.combinedSize()) * sigma).unaryExpr(function<double(double)>(randNorm<double>));
						signalVols[s].slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = (signal + noise).cast<complex<float>>();
					}
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
	for (size_t i = 0; i < sequences.count(); i++) {
		size_t thisSize = sequences.sequence(i)->size();
		templateHdr.setDim(4, thisSize);
		Nifti::File saveFile(templateHdr, outPrefix + filenames[i] + OutExt());
		auto thisSignal = signalVols[i].slice<4>({0,0,0,startVol},{-1,-1,-1,thisSize});
		saveFile.writeVolumes(thisSignal.begin(), thisSignal.end());
		saveFile.close();
	}
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

