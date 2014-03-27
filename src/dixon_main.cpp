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
#include <Eigen/Dense>

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
"Usage is: dixon [options] I0 I1 I2 \n\
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

	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvo:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
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
	if ((argc - optind) != 1) {
		cout << "Requires 1 complex input file with 3 echos." << endl << usage << endl;
		exit(EXIT_FAILURE);
	}
	Series<complex<float>> All;
	
	cout << "Opening input file: " << argv[optind] << endl;
	Nifti inputFile;
	inputFile.open(argv[optind++], Nifti::Mode::Read);
	Nifti templateFile(inputFile, 1);
	All.readFrom(inputFile);
	inputFile.close();

	Volume<float> Wv(templateFile.dims().head(3)), Fv(templateFile.dims().head(3));
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	ThreadPool pool;
	for (size_t k = 0; k < templateFile.dim(3); k++) {
		clock_t loopStart;
		// Read in data
		if (verbose)
			cout << "Starting slice " << k << "..." << flush;
		loopStart = clock();
		atomic<int> voxCount{0};
		
		auto I0s = All.viewSlice(0).viewSlice(k),
		     I1s = All.viewSlice(1).viewSlice(k),
			 I2s = All.viewSlice(2).viewSlice(k);
		auto Ws = Wv.viewSlice(k), Fs = Fv.viewSlice(k);
		
		for (size_t i = 0; i < I0s.dims().prod(); i++) {
			function<void (const size_t)> processVox = [&] (const size_t i) {
				complex<float> I0 = I0s[i], I1 = I1s[i], I2 = I2s[i];
				//float I0d = abs(I0);
				float phi_0 = arg(I0);
				
				complex<float> I1d = I1 * polar<float>(1., -phi_0);
				complex<float> I2d = I2 * polar<float>(1., -phi_0);
				
				float phi = arg(I2d) / 2.;
				float pc = (I1d * polar<float>(1., -phi)).real() / abs(I1d);
				
				Ws[i] = (sqrt(abs(I0)*abs(I2)) + pc * abs(I1)) / 2;
				Fs[i] = (sqrt(abs(I0)*abs(I2)) - pc * abs(I1)) / 2;
			};
			
			pool.for_loop(processVox, templateFile.dim(1));
		}
		
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
	templateFile.open(outPrefix + "water.nii.gz", Nifti::Mode::Write);
	Wv.writeTo(templateFile);
	templateFile.close();
	templateFile.open(outPrefix + "fat.nii.gz", Nifti::Mode::Write);
	Fv.writeTo(templateFile);
	templateFile.close();
	cout << "All done." << endl;
	exit(EXIT_SUCCESS);
}

