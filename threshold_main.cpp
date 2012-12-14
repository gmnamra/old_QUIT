/*
 *  threshold_main.cpp
 *
 *  Created by Tobias Wood on 08/05/2012.
 *  Copyright (c) 2012 Tobias Wood. All rights reserved.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "NiftiImage.h"

const char *usage = "Usage is: threshold [options] input_file threshold output_file\n\
\
Options:\n\
	-v N   : Use volume N from input file.\n\
	-x/y/z L H : Only mask between planes L and H (Use -1 for end).\n";
	
int main(int argc, const char * argv[])
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	if (argc < 4)
	{
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	
	int thisArg = 1;
	unsigned int volume = 0, xl =  0, yl =  0, zl =  0,
					 xh = -1, yh = -1, zh = -1;
	while ((thisArg < argc) && (argv[thisArg][0] =='-'))
	{
		if (strcmp(argv[thisArg], "-v") == 0) {
			volume = atoi(argv[++thisArg]);
		} else if (strcmp(argv[thisArg], "-x") == 0) {
			xl = atoi(argv[++thisArg]);
			xh = atoi(argv[++thisArg]);
		} else if (strcmp(argv[thisArg], "-y") == 0) {
			yl = atoi(argv[++thisArg]);
			yh = atoi(argv[++thisArg]);
		} else if (strcmp(argv[thisArg], "-z") == 0) {
			zl = atoi(argv[++thisArg]);
			zh = atoi(argv[++thisArg]);
		} else {
			fprintf(stderr, "Undefined command line option\n%s", usage);
			exit(EXIT_FAILURE);
		}
		++thisArg;
	}
	
	NiftiImage image;
	image.open(argv[thisArg], 'r');
	fprintf(stdout, "Opened file to mask %s.\n", argv[thisArg]);
	if (volume >= image.nt())
		volume = image.nt() - 1;
	if (xl >= image.nx()) xl = image.nx();
	if (xh >= image.nx()) xh = image.nx();
	if (yl >= image.ny()) yl = image.ny();
	if (yh >= image.ny()) yh = image.ny();
	if (zl >= image.nz()) zl = image.nz();
	if (zh >= image.nz()) zh = image.nz();
	fprintf(stdout, "x %d %d y %d %d z %d %d\n", xl, xh,
	                                             yl, yh,
												 zl, zh);
	double *data = image.readVolume<double>(volume);
	image.close();
	float *mask = (float *)calloc(sizeof(float), image.voxelsPerVolume());
	double thresh = atof(argv[thisArg + 1]);
	fprintf(stdout, "Threshold is %f.\n", thresh);
	for (int i = 0; i < image.voxelsPerVolume(); i++)
	{
		if (data[i] >= thresh)
			mask[i] = 1;
	}
	
	image.setDims(image.nx(), image.ny(), image.nz(), 1);
	image.setDatatype(DT_FLOAT);
	image.open(argv[thisArg + 2], 'w');
	image.writeVolume<float>(0, mask);
	image.close();
	free(data);
	free(mask);
    return EXIT_SUCCESS;
}

