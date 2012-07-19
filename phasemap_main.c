//
//  phasemap_main.c
//  DESPOT
//
//  Created by Tobias Wood on 20/06/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#include "fslio.h"
#include "procpar.h"

char *usage = "Usage is: phasemap input_1 input_2 output_file\n\
\n\
Echo times will be read from procpar if present.\n\
Options:\n\
	--mask mask_file : Mask input with specified file\n\
	--phasetime T    : Calculate the phase accumulated in time T\n";
	
int main(int argc, char** argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	if (argc < 5)
	{
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	
	static struct option long_options[] =
	{
		{"mask", required_argument, 0, 'm'},
		{"phasetime", required_argument, 0, 'p'},
		{0, 0, 0, 0}
	};
	
	int indexptr = 0, c;
	FSLIO *maskHdr = NULL;
	char procpar[MAXSTR];
	par_t *pars;
	double TE1, TE2, deltaTE, phasetime = 0.;
	
	while ((c = getopt_long(argc, argv, "", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 'm':
				maskHdr = FslOpen(optarg, "rb");
				break;
			case 'p':
				phasetime = atof(optarg);
				break;
		}
	}
	
	fprintf(stdout, "Opening input file 1 %s...\n", argv[optind]);
	FSLIO *in1 = FslOpen(argv[optind], "rb");
	strncpy(procpar, argv[optind], MAXSTR);
	strcat(procpar, ".procpar");
	if ((pars = readProcpar(procpar)))
		TE1 = realVal(pars, "te", 0);
	else {
		fprintf(stdout, "Enter TE1 (seconds): ");
		fscanf(stdin, "%lf", &TE1);
	}
	fprintf(stdout, "Opening input file 2 %s...\n", argv[++optind]);
	FSLIO *in2 = FslOpen(argv[optind], "rb");
	strncpy(procpar, argv[optind], MAXSTR);
	strcat(procpar, ".procpar");
	if ((pars = readProcpar(procpar)))
		TE2 = realVal(pars, "te", 0);
	else {
		fprintf(stdout, "Enter TE2 (seconds): ");
		fscanf(stdin, "%lf", &TE2);
	}

	FSLIO *out = FslOpen(argv[++optind], "wb");
	if (optind != argc - 1)
		fprintf(stderr, "Error during argument processing.\n");
	if (TE2 < TE1)
	{	// Swap them
		fprintf(stdout, "TE2 < TE1, swapping.\n");
		FSLIO *tmp = in2;
		in2 = in1;
		in1 = tmp;
		double tmpTE = TE2;
		TE2 = TE1;
		TE1 = tmpTE;
	}
	deltaTE = TE2 - TE1;
	fprintf(stdout, "Delta TE = %f s\n", deltaTE);
	short nx, ny, nz, nvol;
	FslGetDim(in1, &nx, &ny, &nz, &nvol);
	int nvox = nx * ny * nz;
	FslGetDim(in2, &nx, &ny, &nz, &nvol);
	if (nvox != nx * ny * nz)
	{
		fprintf(stderr, "File dimensions do not match.\n");
		exit(EXIT_FAILURE);
	}
	fprintf(stdout, "Image dimensions: %d %d %d\n", nx, ny, nz);
	
	double ***data1 = FslGetVolumeAsScaledDouble(in1, 0);
	fprintf(stdout, "Read image 1.\n");
	double ***data2 = FslGetVolumeAsScaledDouble(in2, 0);
	fprintf(stdout, "Read image 2.\n");
	double ***B0    = d3matrix(nz - 1, ny - 1, nx - 1);
	fprintf(stdout, "Allocated output memory.\n");
	double ***mask = NULL;
	if (maskHdr)
	{
		fprintf(stdout, "Reading mask.\n");
		mask = FslGetVolumeAsScaledDouble(maskHdr, 0);
	}
	
	fprintf(stdout, "Processing...");
	for (size_t z = 0; z < nz; z++)
	{	for (size_t y = 0; y < ny; y++)
		{	for (size_t x = 0; x < nx; x++)
			{
				if (!mask || mask[z][y][x] > 0.)
				{
					double deltaPhase = data2[z][y][x] - data1[z][y][x];
					B0[z][y][x] = deltaPhase / (2 * M_PI * deltaTE);
					if (phasetime > 0.)
					{
						double ph = fmod(B0[z][y][x] * 2 * M_PI * phasetime, 2 * M_PI);
						if (ph > M_PI) ph -= (2 * M_PI);
						if (ph < -M_PI) ph += (2 * M_PI);
						B0[z][y][x] = ph;
					}
				}
			}
		}
	}
	fprintf(stdout, "done.\n");
	FslCloneHeader(out, in1);
	FslSetDim(out, nx, ny, nz, 1);
	FslSetDataType(out, NIFTI_TYPE_FLOAT32);
	FslWriteHeader(out);
	FslWriteVolumeFromDouble(out, B0, 0);
	FslClose(in1);
	FslClose(in2);
	FslClose(out);
	if (maskHdr)
		FslClose(maskHdr);
	fprintf(stdout, "Wrote B0 map succesfully.\n");
    return EXIT_SUCCESS;
}

