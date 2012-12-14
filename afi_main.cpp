/*
 *  afi_main.cpp
 *
 *  Created by Tobias Wood on 03/08/2012.
 *  Copyright (c) 2012 Tobias Wood. All rights reserved.
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <getopt.h>

#include "NiftiImage.h"
#include "procpar.h"

const std::string usage("Usage is: afi [options] input output \n\
\
Options:\n\
	--mask, -m file  : Mask input with specified file.\n\
	--smooth         : Smooth output with a gaussian.\n");
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	static int smooth = false;
	static struct option long_options[] =
	{
		{"mask", required_argument, 0, 'm'},
		{"smooth", no_argument, &smooth, true},
		{0, 0, 0, 0}
	};
	
	int indexptr = 0, c;
	std::string procPath, outPrefix;
	par_t *pars;
	double n, nomFlip;
	double *tr1, *tr2, *flip, *B1, *mask = NULL;
	NiftiImage inFile;
	while ((c = getopt_long(argc, argv, "m:", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 'm':
				std::cout << "Reading mask." << std::endl;
				inFile.open(optarg, NiftiImage::NIFTI_READ);
				mask = inFile.readVolume<double>(0);
				inFile.close();
				break;
		}
	}
	if ((argc - optind) != 2)
	{
		std::cout << usage << std::endl;
		exit(EXIT_FAILURE);
	}
	std::cout << "Opening input file " << argv[optind] << std::endl;
	inFile.open(argv[optind], NiftiImage::NIFTI_READ);
	procPath = inFile.basename() + ".procpar";
	if ((pars = readProcpar(procPath.c_str())))
	{
		// From Sam Hurley. The sequence is implemented by waiting afi_dummy
		// periods after the first afi_tr.
		n = realVal(pars, "afi_dummy", 0) + 1;
		nomFlip = realVal(pars, "flip1", 0);
		fprintf(stdout, "Read TR2/TR1 ratio of %f and flip-angle %f degrees from procpar.\n", n, nomFlip);
	}
	else
	{
		fprintf(stdout, "Enter TR2/TR1 (ratio) and flip-angle (degrees): ");
		fscanf(stdin, "%lf %lf", &n, &nomFlip);
	}
	nomFlip = nomFlip * M_PI / 180.;
	tr1 = inFile.readVolume<double>(0);
	tr2 = inFile.readVolume<double>(1);
	inFile.close();
	outPrefix = std::string(argv[++optind]);
	flip = (double *)malloc(inFile.voxelsPerVolume() * sizeof(double));
	B1   = (double *)malloc(inFile.voxelsPerVolume() * sizeof(double));
	std::cout << "Allocated output memory." << std::endl;
	std::cout << "Processing..." << std::endl;
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
	std::string outPath = outPrefix + "_flip.nii.gz";
	std::cout << "Writing actual flip angle to " << outPath << "..." << std::endl;
	outFile.setnt(1);
	outFile.setDatatype(NIFTI_TYPE_FLOAT32);
	outFile.open(outPath, NiftiImage::NIFTI_WRITE);
	outFile.writeVolume(0, flip);
	outFile.close();
	
	outPath = outPrefix + "_B1.nii.gz";
	std::cout << "Writing B1 ratio to " << outPath << "..." << std::endl;
	outFile.open(outPath, NiftiImage::NIFTI_WRITE);
	outFile.writeVolume(0, B1);
	outFile.close();
	
	std::cout << "Finished." << std::endl;
    return EXIT_SUCCESS;
}

