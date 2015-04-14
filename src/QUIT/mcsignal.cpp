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
#include "Models.h"
#include "Sequence.h"

using namespace std;
using namespace Eigen;
using namespace QUIT;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: mcsignal [options]\n\
\n\
Calculates multi-component DESPOT signals (mainly for testing purposes).\n\
The program will prompt for input (unless --no-prompt specified)\n\
\n\
All times (TR) are in SECONDS. All angles are in degrees.\n\
\n\
Options:\n\
	--help, -h        : Print this message.\n\
	--verbose, -v     : Print extra information.\n\
	--mask, -m file   : Only calculate inside the mask.\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--no-prompt, -n   : Don't print prompts for input.\n\
	--noise, -N val   : Add complex noise with std=val.\n\
	--1, --2, --3     : Use 1, 2 or 3 component sequences (default 3).\n\
	--sequences, -M s : Use simple sequences (default).\n\
	            f     : Use Finite Pulse Length correction.\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};

static shared_ptr<Model> model = make_shared<SCD>();
static bool verbose = false, prompt = true, finitesequences = false;
static string outPrefix = "";
static double sigma = 0.;
static struct option long_opts[] = {
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
	{"threads", required_argument, 0, 'T'},
	{0, 0, 0, 0}
};
static const char *short_opts = "hvnN:m:o:123M:T:";
//******************************************************************************
#pragma mark Read in all required files and data from cin
//******************************************************************************
void parseInput(vector<shared_ptr<SequenceBase>> &cs, vector<string> &names);
void parseInput(vector<shared_ptr<SequenceBase>> &cs, vector<string> &names) {
	string type;
	if (prompt) cout << "Specify next signal type (SPGR/SSFP): " << flush;
	while (Read(cin, type) && (type != "END") && (type != "")) {
		if (type == "SPGR") {
			cs.push_back(make_shared<SPGRSimple>(prompt));
		} else if (type == "SPGRFinite") {
			cs.push_back(make_shared<SPGRFinite>(prompt));
		} else if (type == "SSFP") {
			cs.push_back(make_shared<SSFPSimple>(prompt));
		} else if (type == "SSFPFinite") {
			cs.push_back(make_shared<SSFPFinite>(prompt));
		} else if (type == "SSFPEllipse") {
			cs.push_back(make_shared<SSFPEllipse>(prompt));
		} else if (type == "IRSPGR") {
			cs.push_back(make_shared<IRSPGR>(prompt));
		} else if (type == "MPRAGE") {
			cs.push_back(make_shared<MPRAGE>(prompt));
		} else if (type == "SPINECHO") {
			cs.push_back(make_shared<MultiEcho>(prompt));
		} else {
			throw(std::runtime_error("Unknown signal type: " + type));
		}
		string filename;
		if (prompt) cout << "Enter output filename: " << flush;
		Read(cin, filename);
		names.push_back(filename);
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
	ThreadPool threads;

	try { // To fix uncaught exceptions on Mac
	
	Nifti::File maskFile, B1File;
	MultiArray<int8_t, 3> maskVol;
	MultiArray<float, 3> B1Vol;
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
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
			case '1': model = make_shared<SCD>(); break;
			case '2': model = make_shared<MCD2>(); break;
			case '3': model = make_shared<MCD3>(); break;
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
			case 'T': threads.resize(atoi(optarg)); break;
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

	/***************************************************************************
	 * mark Read in parameter files
	 **************************************************************************/
	cout << "Using " << model->Name() << " model." << endl;
	MultiArray<float, 4> paramsVols;
	Nifti::Header templateHdr;
	if (prompt) cout << "Loading parameters." << endl;
	for (size_t i = 0; i < model->nParameters(); i++) {
		if (prompt) cout << "Enter path to " << model->Names()[i] << " file: " << flush;
		string filename;
		getline(cin, filename);
		cout << "Opening " << filename << endl;
		Nifti::File input(filename);

		if (i == 0) {
			templateHdr = input.header();
			paramsVols = MultiArray<float, 4>(input.matrix(), model->nParameters());
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

	/***************************************************************************
	 * Set up sequences
	 **************************************************************************/
	vector<shared_ptr<SequenceBase>> sequences;
	vector<string> filenames;
	parseInput(sequences, filenames);
	for (auto& s : sequences) {
		cout << *s << endl;
	}

	vector<MultiArray<complex<float>, 4>> signalVols(sequences.size()); //d.head(3), sequences.combinedSize());
	for (size_t s = 0; s < sequences.size(); s++) {
		signalVols[s] = MultiArray<complex<float>, 4>(d.head(3), sequences.at(s)->size());
	}
	cout << "Calculating..." << endl;
	function<void (const size_t&)> calcVox = [&] (const size_t &k) {
		for (size_t j = 0; j < d[1]; j++) {
			for (size_t i = 0; i < d[0]; i++) {
				if (!maskFile || (maskVol[{i,j,k}])) {
					ArrayXd params = paramsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray().cast<double>();
					for (size_t s = 0; s < sequences.size(); s++) {
						const size_t sigsize = sequences.at(s)->size();
						ArrayXcd signal = sequences.at(s)->signal(model, params);
						ArrayXcd noise(sigsize);
						noise.real() = (ArrayXd::Ones(sigsize) * sigma).unaryExpr(function<double(double)>(randNorm<double>));
						noise.imag() = (ArrayXd::Ones(sigsize) * sigma).unaryExpr(function<double(double)>(randNorm<double>));
						signalVols[s].slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = (signal + noise).cast<complex<float>>();
					}
				}
			}
		}
	};
	threads.for_loop(calcVox, d[2]);
	
	cout << "Finished calculating." << endl;
	cout << "Saving data." << endl;
	templateHdr.setDatatype(Nifti::DataType::COMPLEX64);
	size_t startVol = 0;
	for (size_t i = 0; i < sequences.size(); i++) {
		size_t thisSize = sequences.at(i)->size();
		templateHdr.setDim(4, thisSize);
		Nifti::File saveFile(templateHdr, outPrefix + filenames[i]);
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

