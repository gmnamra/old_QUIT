/*
 *  multiecho.cpp
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
"Usage is: multiecho [options] input_file \n\
\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--no-prompt, -n   : Suppress input prompts\n\
	--out, -o path    : Add a prefix to the output filenames\n\
	--mask, -m file   : Mask input with specified file\n\
	--star, -S        : Data is T2*, not T2\n\
	--sum, -s         : Output sum images\n\
	--weighted, -w ?t : Output weighted sum (can fix T2, default average)\n\
	--thresh, -t n    : Threshold maps at PD < n\n\
	--clamp, -c n     : Clamp T2 between 0 and n\n\
	--algo, -a l      : LLS algorithm (default)\n\
	           a      : ARLO algorithm\n\
	           n      : Non-linear (Levenberg-Marquardt)\n\
	--its, -i N       : Max iterations for non-linear (default 10)\n\
	--resids, -r      : Write out per flip-angle residuals\n\
	--threads, -T N   : Use N threads (default=hardware limit)\n"
};

enum class Algo { LogLin, ARLO, Nonlin };
static Algo algo = Algo::LogLin;
static int NE = 0, nIterations = 10;
static bool verbose = false, prompt = true, all_residuals = false, sum = false, weightedSum = false;
static string outPrefix, suffix;
static double weightT2 = 0;
static double thresh = -numeric_limits<double>::infinity();
static double clamp_lo = -numeric_limits<double>::infinity(), clamp_hi = numeric_limits<double>::infinity();
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'n'},
	{"out", required_argument, 0, 'o'},
	{"mask", required_argument, 0, 'm'},
	{"star", no_argument, 0, 'S'},
	{"sum", no_argument, 0, 's'},
	{"weighted", optional_argument, 0, 'w'},
	{"thresh", required_argument, 0, 't'},
	{"clamp", required_argument, 0, 'c'},
	{"algo", required_argument, 0, 'a'},
	{"its", required_argument, 0, 'i'},
	{"threads", required_argument, 0, 'T'},
	{"resids", no_argument, 0, 'r'},
	{0, 0, 0, 0}
};
static const char *short_opts = "hvnm:Ssw::e:o:b:t:c:a:i:T:r";

class RelaxFunctor : public DenseFunctor<double> {
	protected:
		const SequenceBase &m_sequence;
		const ArrayXd m_data;
		const bool m_debug;
		const shared_ptr<SCD> m_model = make_shared<SCD>();

	public:
		RelaxFunctor(SequenceBase &cs, const ArrayXd &data,
		          const double B1, const bool debug) :
			DenseFunctor<double>(2, cs.size()),
			m_sequence(cs), m_data(data),
			m_debug(debug)
		{
			assert(static_cast<size_t>(m_data.rows()) == values());
		}

		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
			eigen_assert(diffs.size() == values());
			VectorXd fullp(5);
			fullp << params(0), 0, params(1), 0, 1.0; // Fix B1 to 1.0 for now
			ArrayXcd s = m_sequence.signal(m_model, fullp);
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
	cout << version << endl << credit_me << endl;
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
			case 'S': suffix = "star"; break;
			case 's': sum = true; break;
			case 'w':
				weightedSum = true;
				if (optarg)
					weightT2 = atof(optarg);
				break;
			case 'i': nIterations = atoi(optarg); break;
			case 't': thresh = atof(optarg); break;
			case 'c':
				clamp_lo = 0;
				clamp_hi = atof(optarg);
				break;
			case 'a':
				switch (*optarg) {
					case 'l': algo = Algo::LogLin; if (verbose) cout << "LogLin algorithm selected." << endl; break;
					case 'a': algo = Algo::ARLO; if (verbose) cout << "ARLO algorithm selected." << endl; break;
					case 'n': algo = Algo::Nonlin; if (verbose) cout << "Non-linear algorithm (Levenberg Marquardt) selected." << endl; break;
					default:
						cout << "Unknown algorithm type " << optarg << endl;
						return EXIT_FAILURE;
						break;
				} break;
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
	MultiEcho multiecho(prompt, pp);
	// Check that NE makes sense
	int NE = multiecho.size();
	if ((inputFile.dim(4) % NE) != 0) {
		throw(runtime_error("Number of volumes is not a multiple of NE."));
	}
	int NVols = inputFile.dim(4) / NE;

	// Set up echo times array
	MatrixXd X(multiecho.size(), 2);
	X.col(0) = multiecho.m_TE;
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
	MultiArray<float, 4> T2Vol(dims, NVols), PDVol(dims, NVols), ResVol(dims, NVols);
	MultiArray<float, 4> sumVol, weightedVol, ResidsVols;
	if (all_residuals) {
		ResidsVols = MultiArray<float, 4>(dims, inputFile.dim(4));
	}
	if (weightedSum) {
		sumVol = MultiArray<float, 4>(dims, NVols);
		weightedVol = MultiArray<float, 4>(dims, NVols);
	}
	for (size_t k = 0; k < dims[2]; k++) {
		clock_t loopStart;
		if (verbose) cout << "Starting slice " << k << "..." << flush;
		loopStart = clock();
		atomic<int> voxCount{0};
		function<void (const size_t, const size_t)> process = [&] (const size_t i, const size_t j) {
			if (!maskFile || (maskVol[{i,j,k}])) {
				voxCount++;
				for (size_t outVol = 0; outVol < NVols; outVol++) {
					double T2 = 0, PD = 0;
					ArrayXd signal = inputVols.slice<1>({i,j,k,outVol*NE},{0,0,0,NE}).asArray().abs().cast<double>();
					switch (algo) {
					case Algo::LogLin: {
						VectorXd Y = signal.log();
						VectorXd b = (X.transpose() * X).partialPivLu().solve(X.transpose() * Y);
						T2 = -1 / b[0];
						PD = exp(b[1]);
					} break;
					case Algo::ARLO: {
						double si2sum = 0, sidisum = 0;
						for (int i = 0; i < NE - 2; i++) {
							double si = (multiecho.m_ESP / 3) * (signal(i) + 4*signal(i+1) + signal(i+2));
							double di = signal(i) - signal(i+2);
							si2sum += si*si;
							sidisum = si*di;
						}
						T2 = (si2sum + (multiecho.m_ESP/3)*sidisum) / ((multiecho.m_ESP/3)*si2sum + sidisum);
						PD = (signal / (-X.col(0).array() / T2).exp()).mean();
					} break;
					case Algo::Nonlin: {
						RelaxFunctor f(multiecho, signal, 1, false);
						NumericalDiff<RelaxFunctor> nDiff(f);
						LevenbergMarquardt<NumericalDiff<RelaxFunctor>> lm(nDiff);
						lm.setMaxfev(nIterations * (NE + 1));
						VectorXd p(2);
						// Just PD & T2 for now
						// Basic guess of T2=50ms
						p << signal(0), 0.05;
						lm.minimize(p);
						PD = p(0); T2 = p(1);
					}
					}

					if (PD < thresh) {
						PD = 0.;
						T2 = 0.;
					}
					T2 = clamp(T2, clamp_lo, clamp_hi);
					ArrayXd theory = PD * (-X.col(0).array() / T2).exp();
					ArrayXd resids = (signal - theory);
					if (all_residuals) {
						ResidsVols.slice<1>({i,j,k,outVol*NE},{0,0,0,NE}).asArray() = resids.cast<float>();
					}
					const MultiArray<float, 4>::Index idx{i,j,k,outVol};
					T2Vol[idx]  = static_cast<float>(T2);
					PDVol[idx]  = static_cast<float>(PD);
					ResVol[idx] = static_cast<float>(sqrt(resids.square().sum() / resids.rows()) / PD);
					if (sum)
						sumVol[idx] = static_cast<float>(signal.sum());
				}
				if (weightedSum) {
					double avT2 = weightT2;
					if (weightT2 == 0)
						avT2 = T2Vol.slice<1>({i,j,k,0},{0,0,0,NVols}).asArray().mean();
					auto weights = (multiecho.m_TE / avT2) * (-multiecho.m_TE / avT2).exp();
					for (size_t outVol = 0; outVol < NVols; outVol++) {
						ArrayXd signal = inputVols.slice<1>({i,j,k,outVol*NE},{0,0,0,NE}).asArray().abs().cast<double>();
						auto sum = (signal * weights).sum();
						weightedVol[{i,j,k,outVol}] = sum;
					}
				}
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
	outHdr.setDim(4, NVols);
	outHdr.setDatatype(Nifti::DataType::FLOAT32);
	outHdr.intent = Nifti::Intent::Estimate;
	outHdr.intent_name = "T2* (seconds)";
	Nifti::File outFile(outHdr, outPrefix + "T2" + suffix + OutExt());
	outFile.writeVolumes(T2Vol.begin(), T2Vol.end());
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
		outFile.writeVolumes(ResidsVols.begin(), ResidsVols.end());
		outFile.close();
	}
	if (sum) {
		outHdr.intent_name = "Sum";
		outHdr.setDim(4, NVols);
		outFile.setHeader(outHdr);
		outFile.open(outPrefix + "sum" + OutExt(), Nifti::Mode::Write);
		outFile.writeVolumes(sumVol.begin(), sumVol.end());
		outFile.close();
	}
	if (weightedSum) {
		outHdr.intent_name = "Weighted Sum";
		outHdr.setDim(4, NVols);
		outFile.setHeader(outHdr);
		outFile.open(outPrefix + "wsum" + OutExt(), Nifti::Mode::Write);
		outFile.writeVolumes(weightedVol.begin(), weightedVol.end());
		outFile.close();
	}
	cout << "All done." << endl;
	} catch (exception &e) {
		cerr << e.what() << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

