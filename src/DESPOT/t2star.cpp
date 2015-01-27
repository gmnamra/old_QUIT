/*
 *  despot1_main.cpp
 *
 *  Created by Tobias Wood on 27/01/2015.
 *  Copyright (c) 2015 Tobias Wood.
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

#include "Nifti/Nifti.h"
#include "QUIT/QUIT.h"
#include "DESPOT.h"
#include "DESPOT_Functors.h"
#include "unsupported/Eigen/NonLinearOptimization"
#include "unsupported/Eigen/NumericalDiff"

using namespace std;
using namespace Eigen;
using namespace QUIT;

class T2starFunctor : public Functor<double> {
	protected:
		const ArrayXd &m_echotimes;
		const ArrayXd &m_data;
	public:
		const long inputs() const override { return 2; }
		const long values() const override { return m_data.rows(); }

		T2starFunctor(const ArrayXd &echos, const ArrayXd &data) : m_echotimes(echos), m_data(data) {};

		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const override {
			double T2star = params[0];
			double PD = params[1];
			diffs = m_data - PD * (-m_echotimes / T2star).exp();
			return 0;
		};
};

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: t2star [options] input_file \n\
\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--no-prompt, -n   : Suppress input prompts\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--mask, -m file   : Mask input with specified file\n\
	--thresh, -t n    : Threshold maps at PD < n\n\
	--clamp, -c n     : Clamp T2* between 0 and n\n\
	--resids, -r      : Write out per flip-angle residuals\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};

static bool verbose = false, prompt = true, all_residuals = false;
static string outPrefix;
static double thresh = -numeric_limits<double>::infinity();
static double clamp_lo = -numeric_limits<double>::infinity(), clamp_hi = numeric_limits<double>::infinity();
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"thresh", required_argument, 0, 't'},
	{"clamp", required_argument, 0, 'c'},
	{"threads", required_argument, 0, 'T'},
	{"resids", no_argument, 0, 'r'},
	{0, 0, 0, 0}
};
static const char *short_opts = "hvnm:o:b:t:c:T:r";
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv) {
	try { // To fix uncaught exceptions on Mac
	cout << version << endl << credit_shared << endl;
	Eigen::initParallel();
	Nifti::File inputFile, maskFile;
	MultiArray<int8_t, 3> maskVol;
	ThreadPool threads;
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
			case 't': thresh = atof(optarg); break;
			case 'c':
				clamp_lo = 0;
				clamp_hi = atof(optarg);
				break;
			case 'T':
				threads.resize(atoi(optarg));
				break;
			case 'r': all_residuals = true; break;
			case 'h':
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != 1) {
		cout << "Incorrect number of arguments." << endl << usage << endl;
		return EXIT_FAILURE;
	}
	// Gather input data
	cout << "Opening input file: " << argv[optind] << endl;
	inputFile.open(argv[optind], Nifti::Mode::Read);
	checkHeaders(inputFile.header(), {maskFile});
	Agilent::ProcPar pp; ReadPP(inputFile, pp);
	double TE1, ESP;
	if (pp) {
		TE1 = pp.realValue("te");
		ESP = pp.realValue("te2");
	} else {
		if (prompt) cout << "Enter first echo-time: " << flush;
		QUIT::Read<double>::FromLine(cin, TE1);
		if (prompt) cout << "Enter echo spacing: " << flush;
		QUIT::Read<double>::FromLine(cin, ESP);
	}
	// Set up echo times array
	MatrixXd X(inputFile.dim(4), 2);
	X(0, 0) = TE1;
	for (int i = 1; i < X.rows(); i++) {
		X(i, 0) = X(i-1, 0) + ESP;
	}
	X.col(1).setOnes();
	if (verbose) {
		cout << "Ouput prefix will be: " << outPrefix << endl;
		cout << "Echo times are: " << X.col(0).transpose() << endl;
		cout << "Clamp: " << clamp_lo << " " << clamp_hi << endl;
		cout << "Thresh: " << thresh << endl;
	}
	cout << "Reading input data..." << flush;
	MultiArray<complex<float>, 4> inputVols(inputFile.dims().head(4));
	inputFile.readVolumes(inputVols.begin(), inputVols.end());
	cout << "done." << endl;
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	const auto dims = inputFile.matrix();
	MultiArray<float, 3> T2starVol(dims), PDVol(dims), ResVol(dims);
	MultiArray<float, 4> ResidsVols;
	if (all_residuals) {
		ResidsVols = MultiArray<float, 4>(dims, inputFile.dim(4));
	}
	for (size_t k = 0; k < dims[2]; k++) {
		clock_t loopStart;
		if (verbose) cout << "Starting slice " << k << "..." << flush;
		loopStart = clock();
		atomic<int> voxCount{0};
		function<void (const size_t, const size_t)> process = [&] (const size_t i, const size_t j) {
			const MultiArray<float, 3>::Index idx{i,j,k};
			if (!maskFile || (maskVol[idx])) {
				voxCount++;
				double T2star, PD;
				ArrayXd signal = inputVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray().abs().cast<double>();
				VectorXd Y = signal.log();
				VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
				T2star = -1 / b[0];
				PD = exp(b[1]);

				/*T2starFunctor f(X.col(0).array(), signal);
				NumericalDiff<T2starFunctor> nDiff(f);
				LevenbergMarquardt<NumericalDiff<T2starFunctor>> lm(nDiff);
				lm.parameters.maxfev = 20;
				VectorXd p(2);
				p << T2star, PD;
				lm.lmder1(p);
				T2star = p(0); PD = p(1);
				if (PD < thresh) {
					PD = 0.;
					T2star = 0.;
				}*/
				if (PD < thresh) {
					PD = 0.;
					T2star = 0.;
				}
				T2star = clamp(T2star, clamp_lo, clamp_hi);
				ArrayXd theory = PD * (-X.col(0).array() / T2star).exp();
				ArrayXd resids = (signal - theory);
				if (all_residuals) {
					ResidsVols.slice<1>({i,j,k,0},{0,0,0,-1}).asArray() = resids.cast<float>();
				}
				T2starVol[idx]  = static_cast<float>(T2star);
				PDVol[idx]  = static_cast<float>(PD);
				ResVol[idx] = static_cast<float>(sqrt(resids.square().sum() / resids.rows()) / PD);
			}
		};

		threads.for_loop2(process, dims[0], dims[1]);

		if (verbose) {
			clock_t loopEnd = clock();
			if (voxCount > 0)
				cout << voxCount << " unmasked voxels, CPU time per voxel was "
				          << ((loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC)) << " s, ";
			cout << "finished." << endl;
		}
	}

	if (verbose)
		cout << "Writing results." << endl;
	outPrefix = outPrefix + "ME_";
	Nifti::Header outHdr = inputFile.header();
	outHdr.description = version;
	outHdr.setDim(4, 1);
	outHdr.setDatatype(Nifti::DataType::FLOAT32);
	outHdr.intent = Nifti::Intent::Estimate;
	outHdr.intent_name = "T2* (seconds)";
	Nifti::File outFile(outHdr, outPrefix + "T2star" + OutExt());
	outFile.writeVolumes(T2starVol.begin(), T2starVol.end());
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
		outHdr.setDim(4, inputFile.dim(4));
		outFile.setHeader(outHdr);
		outFile.open(outPrefix + "residuals" + OutExt(), Nifti::Mode::Write);
		outFile.writeVolumes(ResidsVols.begin(), ResidsVols.end(), 0, inputFile.dim(4));
		outFile.close();
	}
	cout << "All done." << endl;
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

