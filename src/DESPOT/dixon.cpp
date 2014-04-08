//
//  main.cpp
//  dixon
//
//  Created by Tobias Wood on 25/03/2014.
//
//

#include <time.h>
#include <getopt.h>
#include <iostream>
#include <atomic>
#include "Eigen/Dense"

#include "Nifti/Nifti.h"
#include "Nifti/Volume.h"
#include "ThreadPool.h"
#include "DESPOT.h"

#ifdef AGILENT
#include "procpar.h"
using namespace Agilent;
#endif

using namespace std;
using namespace Eigen;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: dixon [options] magnitude phase \n\
\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--out, -o path    : Add a prefix to the output filenames\n"
};

static bool verbose = false;
static string outPrefix;
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{0, 0, 0, 0}
};
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	Nifti maskFile;
	Volume<uint8_t> maskVol;
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvm:o:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				maskFile.open(optarg, Nifti::Mode::Read);
				maskVol.readFrom(maskFile);
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'h':
			case '?': // getopt will print an error message
				exit(EXIT_FAILURE);
		}
	}
	//**************************************************************************
	#pragma mark Gather data
	//**************************************************************************
	if ((argc - optind) != 2) {
		cout << "Requires 1 magnitude file and 1 phase file with 3 echos each as input." << endl << usage << endl;
		exit(EXIT_FAILURE);
	}

	Series<float> mag, phase;
	cout << "Opening magnitude file: " << argv[optind] << endl;
	Nifti inputFile;
	inputFile.open(argv[optind++], Nifti::Mode::Read);
	Nifti templateFile(inputFile, 1);
	mag.readFrom(inputFile);
	inputFile.close();

	cout << "Opening magnitude file: " << argv[optind] << endl;
	inputFile.open(argv[optind++], Nifti::Mode::Read);
	phase.readFrom(inputFile);

	if (!templateFile.matchesSpace(inputFile) || (maskFile.isOpen() && !templateFile.matchesSpace(maskFile))) {
		cerr << "Input file dimensions or orientations do not match." << endl;
		exit(EXIT_FAILURE);
	}
	inputFile.close();

	Nifti::ArrayXs dims = templateFile.dims().head(3);
	Volume<float> Wv(dims), Fv(dims), Av(dims);
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	ThreadPool pool;
	cout << "Starting processing..." << endl;
	for (size_t k = 0; k < templateFile.dim(3); k++) {
		clock_t loopStart;
		// Read in data
		if (verbose)
			cout << "Starting slice " << k << "..." << flush;
		loopStart = clock();
		atomic<int> voxCount{0};
		
		auto S0s = mag.viewSlice(0).viewSlice(k),
		     S1s = mag.viewSlice(1).viewSlice(k),
			 S2s = mag.viewSlice(2).viewSlice(k),
			 phi0s = phase.viewSlice(0).viewSlice(k),
			 phi1s = phase.viewSlice(1).viewSlice(k),
			 phi2s = phase.viewSlice(2).viewSlice(k);
		auto Ms = maskVol.viewSlice(k);
		auto Ws = Wv.viewSlice(k),
		     Fs = Fv.viewSlice(k),
			 As = Av.viewSlice(k);
		//cout << endl << I0s << endl << I1s << endl << I2s << endl;
		//cout << Ws << endl << Fs << endl << As << endl;
		function<void (const size_t)> processVox = [&] (const size_t i) {
			if (!maskFile.isOpen() || Ms[i]) {
				// From Ma et al JMR 1997
				float S0 = S0s[i], S1 = S1s[i], S2 = S2s[i];
				float phi0 = phi0s[i], phi1 = phi1s[i], phi2 = phi2s[i];
				As[i] = sqrt(S2 / S0);
				float phi = (phi2 - phi0) / 2.;
				float psi = cos((phi1 - phi0) - phi);
				float frac = S1 / sqrt(S0*S2);
				Ws[i] = (1 + psi * frac) * S0 / 2.;
				Fs[i] = (1 - psi * frac) * S0 / 2.;
			}
		};
		pool.for_loop(processVox, S0s.size());
		
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
	templateFile.open(outPrefix + "W.nii.gz", Nifti::Mode::Write);
	Wv.writeTo(templateFile);
	templateFile.close();
	templateFile.open(outPrefix + "F.nii.gz", Nifti::Mode::Write);
	Fv.writeTo(templateFile);
	templateFile.close();
	templateFile.open(outPrefix + "A.nii.gz", Nifti::Mode::Write);
	Av.writeTo(templateFile);
	templateFile.close();
	cout << "All done." << endl;
	exit(EXIT_SUCCESS);
}

