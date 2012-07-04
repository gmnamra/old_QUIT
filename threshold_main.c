//
//  main.c
//  threshold
//
//  Created by Tobias Wood on 08/05/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "fslio.h"

char *usage = "Usage is: threshold [options] input_file threshold output_file\n\
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
	
	FSLIO *inHdr = FslOpen(argv[thisArg], "r");
	size_t ndims;
	short nx, ny, nz, nvol;
	int nvox;
	FslGetDimensionality(inHdr, &ndims);
	FslGetDim(inHdr, &nx, &ny, &nz, &nvol);
	nvox = nx * ny * nz;
	
	fprintf(stdout, "Opened file to mask %s.\n", argv[thisArg]);
	fprintf(stdout, "%zu-D file, nx=%d, ny=%d, nz=%d, nvols=%d, nvox=%d.\n",
					ndims, nx, ny, nz, nvol, nvox);
	if (volume >= nvol)
		volume = nvol - 1;
	if (xl >= nx) xl = nx;
	if (xh >= nx) xh = nx;
	if (yl >= ny) yl = ny;
	if (yh >= ny) yh = ny;
	if (zl >= nz) zl = nz;
	if (zh >= nz) zh = nz;
	fprintf(stdout, "x %d %d y %d %d z %d %d\n", xl, xh,
	                                             yl, yh,
												 zl, zh);
	double ***data = FslGetVolumeAsScaledDouble(inHdr, volume);
	float *mask = calloc(sizeof(float), nvox);
	double thresh = atof(argv[thisArg + 1]);
	fprintf(stdout, "Threshold is %f.\n", thresh);
	
	for (size_t z = zl; z < zh; z++)
	{	for (size_t y = yl; y < yh; y++)
		{	for (size_t x = xl; x < xh; x++)
			{
				if (data[z][y][x] >= thresh)
					mask[(z*ny + y)*nx + x] = 1;
			}
		}
	}
	
	FSLIO *outHdr = FslOpen(argv[thisArg + 2], "wb");// FSL_TYPE_NIFTI_GZ);
	FslCloneHeader(outHdr, inHdr);
	FslSetDim(outHdr, nx, ny, nz, 1);
	FslSetDataType(outHdr, NIFTI_TYPE_FLOAT32);
	FslWriteHeader(outHdr);
	FslSeekVolume(outHdr, 0);
	FslWriteVolumes(outHdr, mask, 1);
	FslClose(inHdr);
	FslClose(outHdr);
    return EXIT_SUCCESS;
}

