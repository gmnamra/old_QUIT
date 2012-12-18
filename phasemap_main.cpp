/*
 *  phasemap_main.cpp
 *
 *  Created by Tobias Wood on 20/06/2012.
 *  Copyright (c) 2012 Tobias Wood. All rights reserved.
 */

#include <string>
#include <getopt.h>

#include "NiftiImage.h"
#include "procpar.h"

using namespace std;
#define USE_PROCPAR
#ifdef USE_PROCPAR
	#include "procpar.h"
	using namespace Recon;
#endif

const string usage
{
"Usage is: phasemap input_1 input_2 outprefix\n\
\n\
Echo times will be read from procpar if present.\n\
Options:\n\
	--mask, -m mask_file : Mask input with specified file\n\
	--phasetime T        : Calculate the phase accumulated in time T\n\
	--smooth             : Smooth output with a gaussian.\n"
};

int main(int argc, char** argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	static int smooth = false;
	static struct option long_options[] =
	{
		{"mask", required_argument, 0, 'm'},
		{"phasetime", required_argument, 0, 'p'},
		{"smooth", no_argument, &smooth, true},
		{0, 0, 0, 0}
	};
	
	int indexptr = 0, c;
	double TE1, TE2, deltaTE, phasetime = 0.;
	double *data1, *data2, *B0, *mask = NULL;
	NiftiImage inFile, outFile;
	while ((c = getopt_long(argc, argv, "m:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'm':
				cout << "Reading mask from " << optarg << endl;
				inFile.open(optarg, NiftiImage::NIFTI_READ);
				mask = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case 'p':
				phasetime = atof(optarg);
				break;
		}
	}
	if ((argc - optind) == 2) {
		cout << "Opening input file " << argv[optind] << "." << endl;
		inFile.open(argv[optind], NiftiImage::NIFTI_READ);
		
		#ifdef USE_PROCPAR
		ParameterList pars;
		if (ReadProcpar(inFile.basename() + ".procpar", pars)) {
			TE1 = RealValue(pars, "te", 0);
			TE2 = RealValue(pars, "te", 1);
		} else
		#endif
		{
			cout << "Enter TE2 & TE2 (seconds): ";
			cin >> TE1 >> TE2;
		}
		data1 = inFile.readVolume<double>(0);
		data2 = inFile.readVolume<double>(1);
		inFile.close();
	} else if ((argc - optind) == 3) {
		cout << "Opening input file 1" << argv[optind] << "." << endl;
		inFile.open(argv[optind], NiftiImage::NIFTI_READ);
		#ifdef USE_PROCPAR
		ParameterList pars;
		if (ReadProcpar(inFile.basename() + ".procpar", pars)) {
			TE1 = RealValue(pars, "te", 0);
		} else
		#endif
		{
			cout << "Enter TE1 (seconds): ";
			cin >> TE1;
		}
		data1 = inFile.readVolume<double>(0);
		inFile.close();
		cout << "Opening input file 2" << argv[++optind] << "." << endl;
		inFile.open(argv[optind], NiftiImage::NIFTI_READ);
		#ifdef USE_PROCPAR
		if (ReadProcpar(inFile.basename() + ".procpar", pars)){
			TE2 = RealValue(pars, "te", 0);
		} else
		#endif
		{
			cout << "Enter TE2 (seconds): ";
			cin >> TE2;
		}
		data2 = inFile.readVolume<double>(1);
		inFile.close();
	} else {
		cerr << usage << endl;
		exit(EXIT_FAILURE);
	}
	string outPrefix(argv[++optind]);
	if (TE2 < TE1) {	// Swap them
		fprintf(stdout, "TE2 < TE1, swapping.\n");
		double *tmp = data1;
		data1 = data2;
		data2 = tmp;
		double tmpTE = TE2;
		TE2 = TE1;
		TE1 = tmpTE;
	}
	deltaTE = TE2 - TE1;
	cout << "Delta TE = " << deltaTE << endl;
	B0 = (double *)malloc(inFile.voxelsPerVolume() * sizeof(double));
	cout << "Processing..." << endl;
	for (size_t vox = 0; vox < inFile.voxelsPerVolume(); vox++) {
		if (!mask || mask[vox] > 0.) {
			double deltaPhase = data2[vox] - data1[vox];
			B0[vox] = deltaPhase / (2 * M_PI * deltaTE);
			if (phasetime > 0.) {
				double ph = fmod(B0[vox] * 2 * M_PI * phasetime, 2 * M_PI);
				if (ph > M_PI) ph -= (2 * M_PI);
				if (ph < -M_PI) ph += (2 * M_PI);
				B0[vox] = ph;
			}
		}
	}
	cout << "Writing B0 map." << endl;
	string outPath = outPrefix + "_B0.nii.gz";
	outFile = inFile;
	outFile.setnt(1);
	outFile.setDatatype(NIFTI_TYPE_FLOAT32);
	outFile.open(outPath, NiftiImage::NIFTI_WRITE);
	outFile.writeVolume(0, B0);
	outFile.close();
	cout << "Finished." << endl;
    return EXIT_SUCCESS;
}

