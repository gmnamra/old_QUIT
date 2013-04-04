/*
 *  afi_main.cpp
 *
 *  Created by Tobias Wood on 03/08/2012.
 *  Copyright (c) 2012-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <getopt.h>

#include "NiftiImage.h"
#include "procpar.h"

using namespace std;
#ifdef HAVE_NRECON
	#include "procpar.h"
	using namespace Recon;
#endif

const string usage("Usage is: afi [options] input output \n\
\
Options:\n\
	--mask, -m file  : Mask input with specified file.\n");
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	static struct option long_options[] =
	{
		{"mask", required_argument, 0, 'm'},
		{0, 0, 0, 0}
	};
	
	int indexptr = 0, c;
	string procPath, outPrefix;
	double n, nomFlip;
	double *tr1, *tr2, *flip, *B1, *mask = NULL;
	NiftiImage inFile;
	while ((c = getopt_long(argc, argv, "m:", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 'm':
				cout << "Reading mask." << endl;
				if (!inFile.open(optarg, NiftiImage::NIFTI_READ)) {
					exit(EXIT_FAILURE);
				}
				mask = inFile.readVolume<double>(0);
				inFile.close();
				break;
		}
	}
	if ((argc - optind) != 2) {
		cout << usage << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Opening input file " << argv[optind] << endl;
	if (!inFile.open(argv[optind], NiftiImage::NIFTI_READ)) {
		exit(EXIT_FAILURE);
	}
	#ifdef HAVE_NRECON
	ParameterList pars;
	if (ReadProcpar(inFile.basename() + ".procpar", pars)) {
		// From Sam Hurley. The sequence is implemented by waiting afi_dummy
		// periods after the first afi_tr.
		n = RealValue(pars, "afi_dummy") + 1;
		nomFlip = RealValue(pars, "flip1");
		cout << "Read TR2/TR1 ratio of " << n << " and flip-angle " << nomFlip << " degrees from procpar." << endl;
	} else
	#endif
	{
		cout << "Enter TR2/TR1 (ratio) and flip-angle (degrees): ";
		cin >> n >> nomFlip;
	}
	nomFlip = nomFlip * M_PI / 180.;
	tr1 = inFile.readVolume<double>(0);
	tr2 = inFile.readVolume<double>(1);
	inFile.close();
	outPrefix = string(argv[++optind]);
	flip = (double *)malloc(inFile.voxelsPerVolume() * sizeof(double));
	B1   = (double *)malloc(inFile.voxelsPerVolume() * sizeof(double));
	cout << "Allocated output memory." << endl;
	cout << "Processing..." << endl;
	for (size_t vox = 0; vox < inFile.voxelsPerVolume(); vox++) {
		if (!mask || mask[vox] > 0.) {
			double r = tr2[vox] / tr1[vox];
			double temp = (r*n - 1.) / (n - r);
			if (temp > 1.)
				temp = 1.;
			if (temp < -1.)
				temp = -1.;
			double alpha = acos(temp);
			flip[vox] = alpha * 180. / M_PI;
			B1[vox]   = alpha / nomFlip;
		}
		else
			B1[vox] = 1.; // So smoothing doesn't get messed up
	}
	NiftiImage outFile = inFile; // Could re-use infile, this is marginally clearer
	string outPath = outPrefix + "_flip.nii.gz";
	cout << "Writing actual flip angle to " << outPath << "..." << endl;
	outFile.setDim(4, 1);
	outFile.setDatatype(NIFTI_TYPE_FLOAT32);
	outFile.open(outPath, NiftiImage::NIFTI_WRITE);
	outFile.writeVolume(0, flip);
	outFile.close();
	
	outPath = outPrefix + "_B1.nii.gz";
	cout << "Writing B1 ratio to " << outPath << "..." << endl;
	outFile.open(outPath, NiftiImage::NIFTI_WRITE);
	outFile.writeVolume(0, B1);
	outFile.close();
	
	cout << "Finished." << endl;
    return EXIT_SUCCESS;
}

