/*
 *  despot1_main.cpp
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright (c) 2011-2013 Tobias Wood.
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
#include "Sequence.h"

using namespace std;
using namespace Eigen;
using namespace QUIT;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: despot1 [options] spgr_input \n\
\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--no-prompt, -n   : Suppress input prompts\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--mask, -m file   : Mask input with specified file\n\
	--B1, -b file     : B1 Map file (ratio)\n\
	--thresh, -t n    : Threshold maps at PD < n\n\
	--clamp, -c n     : Clamp T1 between 0 and n\n\
	--algo, -a l      : LLS algorithm (default)\n\
	           w      : WLLS algorithm\n\
	           n      : NLLS (Levenberg-Marquardt)\n\
	--its, -i N       : Max iterations for WLLS (default 4)\n\
	--resids, -r      : Write out per flip-angle residuals\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};

enum class Algos { LLS, WLLS, NLLS };
static bool verbose = false, prompt = true, all_residuals = false;
static size_t nIterations = 4;
static string outPrefix;
static double thresh = -numeric_limits<double>::infinity();
static double clamp_lo = -numeric_limits<double>::infinity(), clamp_hi = numeric_limits<double>::infinity();
static Algos algo;
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"B1", required_argument, 0, 'b'},
	{"thresh", required_argument, 0, 't'},
	{"clamp", required_argument, 0, 'c'},
	{"algo", required_argument, 0, 'a'},
	{"its", required_argument, 0, 'i'},
	{"threads", required_argument, 0, 'T'},
	{"resids", no_argument, 0, 'r'},
	{0, 0, 0, 0}
};
static const char *short_opts = "hvnm:o:b:t:c:a:i:T:r";

// T1 only Functor
class T1Functor : public DenseFunctor<double> {
	protected:
		const SequenceBase &m_sequence;
		const ArrayXd m_data;
		const bool m_debug;
		const double m_B1;
		const shared_ptr<SCD> m_model;

	public:
		T1Functor(SequenceBase &cs, const ArrayXd &data,
		          const double B1, const bool debug) :
			DenseFunctor<double>(2, cs.size()),
			m_sequence(cs), m_data(data),
			m_B1(B1), m_debug(debug), m_model()
		{
			assert(static_cast<size_t>(m_data.rows()) == values());
		}

		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
			eigen_assert(diffs.size() == values());
			ArrayXcd s = m_sequence.signal(m_model, params, m_B1);
			diffs = s.abs() - m_data;
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
int main(int argc, char **argv) {
	try { // To fix uncaught exceptions on Mac
	cout << version << endl << credit_shared << endl;
	Eigen::initParallel();
	Nifti::File spgrFile, B1File, maskFile;
	MultiArray<float, 3> B1Vol;
	MultiArray<int8_t, 3> maskVol;
	ThreadPool threads;
	shared_ptr<SCD> model = make_shared<SCD>();

	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'n': prompt = false; break;
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
				B1File.readVolumes(B1Vol.begin(), B1Vol.end());
				break;
			case 't': thresh = atof(optarg); break;
			case 'c':
				clamp_lo = 0;
				clamp_hi = atof(optarg);
				break;
			case 'a':
				switch (*optarg) {
					case 'l': algo = Algos::LLS;  cout << "LLS algorithm selected." << endl; break;
					case 'w': algo = Algos::WLLS; cout << "WLLS algorithm selected." << endl; break;
					case 'n': algo = Algos::NLLS; cout << "NLLS algorithm selected." << endl; break;
					default:
						cout << "Unknown algorithm type " << optarg << endl;
						return EXIT_FAILURE;
						break;
				} break;
			case 'i':
				nIterations = atoi(optarg);
				break;
			case 'T':
				threads.resize(atoi(optarg));
				break;
			case 'r': all_residuals = true; break;
			case 'h':
			case '?': // getopt will print an error message
				cout << usage << endl;
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != 1) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		return EXIT_FAILURE;
	}
	//**************************************************************************
	#pragma mark Gather SPGR data
	//**************************************************************************
	cout << "Opening SPGR file: " << argv[optind] << endl;
	spgrFile.open(argv[optind], Nifti::Mode::Read);
	checkHeaders(spgrFile.header(), {maskFile, B1File});
	Agilent::ProcPar pp; ReadPP(spgrFile, pp);
	SPGRSimple spgrSequence(prompt, pp);
	if (verbose) {
		cout << spgrSequence;
		cout << "Ouput prefix will be: " << outPrefix << endl;
		cout << "Clamp: " << clamp_lo << " " << clamp_hi << endl;
		cout << "Thresh: " << thresh << endl;
	}
	if (spgrSequence.size() != spgrFile.header().dim(4)) {
		throw(std::runtime_error("Specified number of flip-angles does not match number of volumes in file: " + spgrFile.imagePath()));
	}
	double TR = spgrSequence.m_TR;
	cout << "Reading SPGR data..." << flush;
	MultiArray<complex<float>, 4> spgrVols(spgrFile.dims().head(4));
	spgrFile.readVolumes(spgrVols.begin(), spgrVols.end());
	cout << "done." << endl;
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	const auto dims = spgrFile.matrix();
	MultiArray<float, 3> T1Vol(dims), PDVol(dims), ResVol(dims);
	MultiArray<float, 4> ResidsVols;
	if (all_residuals) {
		ResidsVols = MultiArray<float, 4>(dims, spgrSequence.size());
	}
	time_t startTime;
	if (verbose)
		startTime = printStartTime();
	clock_t startClock = clock();
	int voxCount = 0;
	for (size_t k = 0; k < spgrFile.dim(3); k++) {
		clock_t loopStart = clock();
		if (verbose) cout << "Starting slice " << k << "..." << flush;
		atomic<int> sliceCount{0};
		function<void (const size_t, const size_t)> process = [&] (const size_t i, const size_t j) {
			const MultiArray<float, 3>::Index idx{i,j,k};
			if (!maskFile || (maskVol[idx])) {
				sliceCount++;
				double B1 = B1File ? B1Vol[idx] : 1.;
				ArrayXd localAngles(spgrSequence.B1flip(B1));
				double T1, PD;
				ArrayXd signal = spgrVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray().abs().cast<double>();
				VectorXd Y = signal / localAngles.sin();
				MatrixXd X(Y.rows(), 2);
				X.col(0) = signal / localAngles.tan();
				X.col(1).setOnes();
				VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
				T1 = -TR / log(b[0]);
				PD = b[1] / (1. - b[0]);
				if (algo == Algos::WLLS) {
					VectorXd W(spgrSequence.size());
					for (size_t n = 0; n < nIterations; n++) {
						W = (localAngles.sin() / (1. - (exp(-TR/T1)*localAngles.cos()))).square();
						b = (X.transpose() * W.asDiagonal() * X).partialPivLu().solve(X.transpose() * W.asDiagonal() * Y);
						T1 = -TR / log(b[0]);
						PD = b[1] / (1. - b[0]);
					}
				} else if (algo == Algos::NLLS) {
					T1Functor f(spgrSequence, signal, B1, false);
					NumericalDiff<T1Functor> nDiff(f);
					LevenbergMarquardt<NumericalDiff<T1Functor>> lm(nDiff);
					lm.setMaxfev(nIterations);
					VectorXd p(4);
					p << PD, T1, 0., 0.; // Don't need T2 of f0 for this (yet)
					lm.lmder1(p);
					PD = p(0); T1 = p(1);
				}
				if (PD < thresh) {
					PD = 0.;
					T1 = 0.;
				}
				T1 = clamp(T1, clamp_lo, clamp_hi);
				ArrayXd theory = spgrSequence.signal(model, Vector4d(PD, T1, 0., 0.), B1).abs();
				ArrayXd resids = (signal - theory);
				if (all_residuals) {
					ResidsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = resids.cast<float>();
				}
				T1Vol[idx]  = static_cast<float>(T1);
				PDVol[idx]  = static_cast<float>(PD);
				ResVol[idx] = static_cast<float>(sqrt(resids.square().sum() / resids.rows()) / PD);
			}
		};
		threads.for_loop2(process, spgrFile.dim(1), spgrFile.dim(2));
		if (verbose) printLoopTime(loopStart, sliceCount);
		voxCount += sliceCount;
		if (threads.interrupted())
			break;
	}
	if (verbose) {
		printElapsedTime(startTime);
		printElapsedClock(startClock, voxCount);
		cout << "Writing results." << endl;
	}
	outPrefix = outPrefix + "D1_";
	Nifti::Header outHdr = spgrFile.header();
	outHdr.description = version;
	outHdr.setDim(4, 1);
	outHdr.setDatatype(Nifti::DataType::FLOAT32);
	outHdr.intent = Nifti::Intent::Estimate;
	outHdr.intent_name = "T1 (seconds)";
	Nifti::File outFile(outHdr, outPrefix + "T1" + OutExt());
	outFile.writeVolumes(T1Vol.begin(), T1Vol.end());
	outFile.close();
	outHdr.intent_name = "PD (au)";
	outFile.setHeader(outHdr);
	outFile.open(outPrefix + "PD" + OutExt(), Nifti::Mode::Write);
	outFile.writeVolumes(PDVol.begin(), PDVol.end());
	outFile.close();
	outHdr.intent_name = "Fractional Residual";
	outFile.setHeader(outHdr);
	outFile.open(outPrefix + "residual" + OutExt(), Nifti::Mode::Write);
	outFile.writeVolumes(ResVol.begin(), ResVol.end());
	outFile.close();
	if (all_residuals) {
		outHdr.intent_name = "Residuals";
		outHdr.setDim(4, spgrSequence.size());
		outFile.setHeader(outHdr);
		outFile.open(outPrefix + "residuals" + OutExt(), Nifti::Mode::Write);
		outFile.writeVolumes(ResidsVols.begin(), ResidsVols.end(), 0, spgrSequence.size());
		outFile.close();
	}
	cout << "All done." << endl;
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
