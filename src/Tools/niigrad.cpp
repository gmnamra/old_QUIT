/*
 *  niigrad.cpp
 *
 *  Created by Tobias Wood on 23/04/2014.
 *  Copyright (c) 2014 Tobias Wood.
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

#include "Nifti/Nifti.h"
#include "QUIT/Volume.h"
#include "QUIT/VolumeFunctions.h"
#include "QUIT/ThreadPool.h"

using namespace std;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: niigrad [options] input\n\
\n\
Default mode is to calculate the gradient (dI/dx + dI/dy + dI/dz) for each voxel\n\
in every volume. If you want to output the specific directional derivatives then\n\
specify the relevant options.\n\
\n\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--grad, -g        : Output summed _grad file (default)\n\
	--deriv, -d       : Output separate _dx, _dy, _dz files\n"
};

static bool verbose = false, outputGrad = false, outputDeriv = false;
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"grad", required_argument, 0, 'g'},
	{"deriv", required_argument, 0, 'd'},
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
	while ((c = getopt_long(argc, argv, "hvd:g:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'g': outputGrad = true; break;
			case 'd': outputDeriv = true; break;
			case 'h':
			case '?': // getopt will print an error message
				exit(EXIT_FAILURE);
		}
	}
	// Make sure we output something
	if (!(outputGrad || outputDeriv))
		outputGrad = true;

	if ((argc - optind) != 1) {
		cerr << "Incorrect number of files to process." << endl;
		exit(EXIT_FAILURE);
	}

	cout << "Opening input file: " << argv[optind] << endl;
	Nifti inFile(argv[optind], Nifti::Mode::Read);
	std::string basename = inFile.basePath();
	Nifti outFile(inFile); outFile.close(); outFile.open(basename + "_grad.nii.gz", Nifti::Mode::Write);
	cout << "Allocating working memory." << endl;
	auto d = inFile.dims().head(3);
	Volume<float> data(d, inFile.transform());
	Volume<float> grad(d, inFile.transform());
	Volume<Eigen::Vector3f> deriv(d, inFile.transform());
	cout << "Processing." << endl;
	for (size_t vol = 0; vol < inFile.dim(4); vol++) {
		inFile.readVolumes(data.data().begin(), data.data().end(), vol, 1);
		cout << "Calculating gradient for volume " << vol << endl;
		VolumeDerivative(data, grad, deriv);
		outFile.writeVolumes(grad.data().begin(), grad.data().end(), vol, 1);
	}
	inFile.close();
	outFile.close();
	cout << "Finished." << endl;
	exit(EXIT_SUCCESS);
}
