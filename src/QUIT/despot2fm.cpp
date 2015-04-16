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
#include <iostream>
#include <atomic>
#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include <unsupported/Eigen/NumericalDiff>

#include "Nifti/Nifti.h"
#include "QUIT/QUIT.h"
#include "Models.h"
#include "Sequence.h"
#include "RegionContraction.h"

using namespace std;
using namespace Eigen;
using namespace QUIT;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: despot2-fm [options] T1_map ssfp_files\n\
\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print slice processing times\n\
	--no-prompt, -n   : Suppress input prompts\n\
	--mask, -m file   : Mask input with specified file\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--f0, -f SYM      : Fit symmetric f0 map (default)\n\
	         ASYM     : Fit asymmetric f0 map\n\
	         file     : Use f0 Map file (in Hertz)\n\
	--B1, -b file     : B1 Map file (ratio)\n\
	--start, -s N     : Start processing from slice N\n\
	--stop, -p  N     : Stop processing at slice N\n\
	--scale, -S 0     : Normalise signals to mean (default)\n\
	            1     : Fit a scaling factor/proton density\n\
	--flip, -F        : Data order is phase, then flip-angle (default opposite)\n\
	--sequences, -M s : Use simple sequences (default)\n\
	            f     : Use finite pulse length correction\n\
	--complex, -x     : Fit to complex data\n\
	--contract, -c n  : Read contraction settings from stdin (Will prompt)\n\
	--resids, -r      : Write out per flip-angle residuals\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};

static auto tesla = FieldStrength::Three;
static auto f0fit = OffRes::FitSym;
static size_t start_slice = 0, stop_slice = numeric_limits<size_t>::max();
static int verbose = false, prompt = true, all_residuals = false,
           fitFinite = false, fitComplex = false, flipData = false,
           samples = 2000, retain = 20, contract = 10,
           voxI = 0, voxJ = 0, seed = -1;
static double expand = 0.;
static string outPrefix;
static struct option long_opts[] = {
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"f0", required_argument, 0, 'f'},
	{"B1", required_argument, 0, 'b'},
	{"start", required_argument, 0, 's'},
	{"stop", required_argument, 0, 'p'},
	{"scale", required_argument, 0, 'S'},
	{"flip", required_argument, 0, 'F'},
	{"threads", required_argument, 0, 'T'},
	{"sequences", no_argument, 0, 'M'},
	{"complex", no_argument, 0, 'x'},
	{"contract", no_argument, 0, 'c'},
	{"resids", no_argument, 0, 'r'},
	{0, 0, 0, 0}
};
static const char* short_opts = "hvnm:o:f:b:s:p:S:FT:M:xcri:j:d:";

//******************************************************************************
// T2 Only Functor
//******************************************************************************
class FMFunctor : public DenseFunctor<double> {
	public:
		SequenceBase &m_sequence;
		shared_ptr<SCD> m_model;
		ArrayXcd m_data;
		const double m_T1, m_B1;
		const bool m_complex, m_debug;

		FMFunctor(const shared_ptr<SCD> m, const double T1, SequenceBase &s, const ArrayXcd &d, const double B1, const bool fitComplex, const bool debug = false) :
			DenseFunctor<double>(3, s.size()),
			m_model(m), m_sequence(s), m_data(d), m_complex(fitComplex), m_debug(debug),
			m_T1(T1), m_B1(B1)
		{
			assert(static_cast<size_t>(m_data.rows()) == values());
		}

		const bool constraint(const VectorXd &params) const {
			Array4d fullparams;
			fullparams << params(0), m_T1, params(1), params(2);
			return m_model->ValidParameters(fullparams);
		}

		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
			eigen_assert(diffs.size() == values());

			ArrayXd fullparams(5);
			fullparams << params(0), m_T1, params(1), params(2), m_B1;
			ArrayXcd s = m_sequence.signal(m_model, fullparams);
			if (m_complex) {
				diffs = (s - m_data).abs();
			} else {
				diffs = s.abs() - m_data.abs();
			}
			if (m_debug) {
				cout << endl << __PRETTY_FUNCTION__ << endl;
				cout << "p:     " << params.transpose() << endl;
				cout << "s:     " << s.transpose() << endl;
				cout << "data:  " << m_data.transpose() << endl;
				cout << "diffs: " << diffs.transpose() << endl;
			}
			return 0;
		}
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	try { // To fix uncaught exceptions on Mac
	cout << version << endl << credit_me << endl;
	Eigen::initParallel();
	Nifti::File maskFile, f0File, B1File;
	MultiArray<int8_t, 3> maskVol;
	MultiArray<float, 3> f0Vol, B1Vol;
	string procPath;
	ThreadPool threads;
	shared_ptr<SCD> model = make_shared<SCD>();
	//ThreadPool::EnableDebug = true;

	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
			case 'm':
				if (verbose) cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskVol.resize(maskFile.matrix());
				maskFile.readVolumes(maskVol.begin(), maskVol.end(), 0, 1);
				break;
			case 'o':
				outPrefix = optarg;
				if (verbose) cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'f':
				if (string(optarg) == "SYM") {
					f0fit = OffRes::FitSym;
				} else if (string(optarg) == "ASYM") {
					f0fit = OffRes::Fit;
				} else {
					if (verbose) cout << "Reading f0 file: " << optarg << endl;
					f0File.open(optarg, Nifti::Mode::Read);
					f0Vol.resize(f0File.dims());
					f0File.readVolumes(f0Vol.begin(), f0Vol.end(), 0, 1);
					f0fit = OffRes::Map;
				}
				break;
			case 'b':
				if (verbose) cout << "Reading B1 file: " << optarg << endl;
				B1File.open(optarg, Nifti::Mode::Read);
				B1Vol.resize(B1File.matrix());
				B1File.readVolumes(B1Vol.begin(), B1Vol.end(), 0, 1);
				break;
			case 's': start_slice = atoi(optarg); break;
			case 'p': stop_slice = atoi(optarg); break;
			case 'S':
				switch (atoi(optarg)) {
					case 0 : model->setScaling(Model::Scale::ToMean); break;
					case 1 : model->setScaling(Model::Scale::None); break;
					default:
						cout << "Invalid scaling mode: " + to_string(atoi(optarg)) << endl;
						return EXIT_FAILURE;
						break;
				} break;
			case 'F':
				flipData = true;
				break;
			case 'T':
				threads.resize(atoi(optarg));
				break;
			case 'd':
				seed = atoi(optarg);
				break;
			case 'M':
				switch (*optarg) {
					case 's': fitFinite = false; cout << "Simple sequences selected." << endl; break;
					case 'f': fitFinite = true; cout << "Finite pulse correction selected." << endl; break;
					default:
						cout << "Unknown sequences type " << *optarg << endl;
						return EXIT_FAILURE;
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
				break;
			case 'r': all_residuals = true; break;
			case 'i': voxI = atoi(optarg); break;
			case 'j': voxJ = atoi(optarg); break;
			case '?': // getopt will print an error message
			case 'h':
			default:
				cout << usage << endl;
				return EXIT_SUCCESS;
				break;
		}
	}
	if ((argc - optind) < 2) {
		cout << "Wrong number of arguments. Need at least a T1 map and 1 SSFP file." << endl;
		return EXIT_FAILURE;
	}
	if (verbose) cout << "Reading T1 Map from: " << argv[optind] << endl;
	Nifti::File T1File(argv[optind++]);
	const auto dims = T1File.matrix();
	MultiArray<float, 3> T1Vol(dims);
	T1File.readVolumes(T1Vol.begin(), T1Vol.end(), 0, 1);
	checkHeaders(T1File.header(), {maskFile, f0File, B1File});
	//**************************************************************************
	// Gather SSFP Data
	//**************************************************************************
	size_t nFiles = argc - optind;
	vector<MultiArray<complex<float>, 4>> ssfpData(nFiles);
	SequenceGroup sequences;
	for (size_t p = 0; p < nFiles; p++) {
		if (verbose) cout << "Reading SSFP header from " << argv[optind] << endl;
		Nifti::File inFile(argv[optind]);
		checkHeaders(inFile.header(), {T1File});
		Agilent::ProcPar pp; ReadPP(inFile, pp);
		if (fitFinite) {
			sequences.addSequence(make_shared<SSFPFinite>(prompt, pp));
		} else {
			sequences.addSequence(make_shared<SSFPSimple>(prompt, pp));
		}
		if (sequences.sequence(sequences.count() - 1)->size() != inFile.dim(4)) {
			throw(std::runtime_error("Number of volumes in file " + inFile.imagePath() + " does not match input."));
		}
		if (verbose) cout << "Reading data." << endl;
		ssfpData.at(p).resize(inFile.dims().head(4));
		inFile.readVolumes(ssfpData.at(p).begin(), ssfpData.at(p).end());
		inFile.close();
		optind++;
	}
	if (optind != argc) {
		cerr << "Unprocessed arguments supplied.\n" << usage;
		return EXIT_FAILURE;
	}

	ArrayXd thresh(3); thresh.setConstant(0.05);
	ArrayXd weights(sequences.size()); weights.setOnes();
	Array2d f0Bounds(-0.5/sequences.minTR(),0.5/sequences.minTR());
	if (f0fit == OffRes::FitSym) {
		f0Bounds(0) = 0.;
	}
	
	if (verbose) {
		cout << sequences;
		if (flipData) {
			cout << "Data order is phase, then flip-angle." << endl;
		} else {
			cout << "Data order is flip-angle, then phase." << endl;
		}
	}
	//**************************************************************************
	// Set up results data
	//**************************************************************************
	size_t nParams = 2;
	if (model->scaling() == Model::Scale::None)
		nParams = 3;
	MultiArray<float, 4> paramsVols(dims, nParams);
	MultiArray<float, 4> ResidsVols(dims, sequences.size());
	MultiArray<float, 3> ResVol(dims);
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	if (stop_slice > dims[2])
		stop_slice = dims[2];
	time_t startTime;
	if (verbose) startTime = printStartTime();
	clock_t startClock = clock();
	int voxCount = 0;
	for (size_t k = start_slice; k < stop_slice; k++) {
		if (verbose) cout << "Starting slice " << k << "..." << flush;
		atomic<int> sliceCount{0};
		clock_t loopStart = clock();

		function<void (const size_t, const size_t)>
		processVox = [&] (const size_t i, const size_t j) {
			const MultiArray<float, 3>::Index idx{i,j,k};
			if (!maskFile || (maskVol[idx] && T1Vol[idx] > 0.)) {
				// -ve T1 is nonsensical, no point fitting
				ArrayXcd signal = sequences.loadSignals(ssfpData, i, j, k, flipData);
				ArrayXXd bounds(3, 2);
				bounds.setZero();
				if (model->scaling() == Model::Scale::None) {
					bounds(0, 0) = 0.;
					bounds(0, 1) = signal.abs().maxCoeff() * 25;
				} else {
					bounds.row(0).setConstant(1.);
				}
				bounds(1,0) = 0.001;
				bounds(1,1) = T1Vol[idx];
				if (f0fit == OffRes::Map) {
					bounds.row(2).setConstant(f0Vol[idx]);
				} else {
					bounds.row(2) = f0Bounds;
				}
				double B1 = B1File ? B1Vol[{i,j,k}] : 1.;
				FMFunctor func(model, T1Vol[idx], sequences, signal, B1, fitComplex, false);
				RegionContraction<FMFunctor> rc(func, bounds, weights, thresh,
				                                samples, retain, contract, expand, (voxI > 0), seed);
				ArrayXd params(3); params.setZero();
				rc.optimise(params);
				double res = sqrt(rc.SoS() / sequences.size());
				if (model->scaling() == Model::Scale::None) {
					paramsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = params.cast<float>();
					res /= params(0); // Divide residual by PD to make it a fraction
				} else {
					paramsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = params.tail(2).cast<float>(); // Skip PD
				}
				ResVol[{i,j,k}] = static_cast<float>(res);
				if (all_residuals) {
					ResidsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = rc.residuals().cast<float>();
				}
				if ((rc.status() == RCStatus::Converged) || (rc.status() == RCStatus::IterationLimit)) {
					sliceCount++;
				}
			}
		};
		if (voxI == 0) {
			threads.for_loop2(processVox, dims[0], dims[1]);
			if (threads.interrupted())
				break;
		} else {
			processVox(voxI, voxJ);
			voxCount = 1;
			break;
		}
		
		if (verbose) printLoopTime(loopStart, sliceCount);
		voxCount += sliceCount;

		if (threads.interrupted())
			break;
	}
	if (verbose) printElapsedTime(startTime);
	printElapsedClock(startClock, voxCount);
	if (voxI != 0)
		return EXIT_SUCCESS;

	outPrefix = outPrefix + "FM_";
	Nifti::Header hdr = T1File.header();
	hdr.setDim(4, 1);
	hdr.setDatatype(Nifti::DataType::FLOAT32);
	hdr.description = version;
	hdr.intent = Nifti::Intent::Estimate;
	if (model->scaling() == Model::Scale::None) {
		hdr.intent_name = model->Names().at(0);
		Nifti::File out(hdr, outPrefix + model->Names().at(0) + OutExt());
		auto p = paramsVols.slice<3>({0,0,0,0},{-1,-1,-1,0});
		out.writeVolumes(p.begin(), p.end());
		out.close();
		hdr.intent_name = model->Names().at(2);
		out.open(outPrefix + model->Names().at(2) + OutExt(), Nifti::Mode::Write);
		p = paramsVols.slice<3>({0,0,0,1},{-1,-1,-1,0});
		out.writeVolumes(p.begin(), p.end());
		out.close();
		hdr.intent_name = model->Names().at(3);
		out.open(outPrefix + model->Names().at(3) + OutExt(), Nifti::Mode::Write);
		p = paramsVols.slice<3>({0,0,0,2},{-1,-1,-1,0});
		out.writeVolumes(p.begin(), p.end());
		out.close();
	} else {
		hdr.intent_name = model->Names().at(2);
		Nifti::File out(hdr, outPrefix + model->Names().at(2) + OutExt());
		auto p = paramsVols.slice<3>({0,0,0,0},{-1,-1,-1,0});
		out.writeVolumes(p.begin(), p.end());
		out.close();
		hdr.intent_name = model->Names().at(3);
		out.open(outPrefix + model->Names().at(3) + OutExt(), Nifti::Mode::Write);
		p = paramsVols.slice<3>({0,0,0,1},{-1,-1,-1,0});
		out.writeVolumes(p.begin(), p.end());
		out.close();
	}
	hdr.intent_name = "Sum of Squared Residuals";
	Nifti::File SoS(hdr, outPrefix + "residual" + OutExt());
	SoS.writeVolumes(ResVol.begin(), ResVol.end());
	SoS.close();
	if (all_residuals) {
		hdr.setDim(4, static_cast<int>(sequences.size()));
		Nifti::File res(hdr, outPrefix + "residuals" + OutExt());
		res.writeVolumes(ResidsVols.begin(), ResidsVols.end());
		res.close();
	}
	
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
