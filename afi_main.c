/*
 *  afi_main.c
 *  DESPOT
 *
 *  Created by Tobias Wood on 03/08/2012.
 *  Copyright (c) 2012 Tobias Wood. All rights reserved.
*/

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <getopt.h>
#include <dispatch/dispatch.h>
#ifdef __APPLE__
	#include <libkern/OSAtomic.h>
	#define AtomicAdd OSAtomicAdd32
#else
	#define AtomicAdd(x, y) (*y) += x
#endif
#include "fslio.h"
#include "procpar.h"
#include "mathsUtil.h"
#include "mathsOps.h"

char *usage = "Usage is: afi [options] input output \n\
\
Options:\n\
	-m file  : Mask input with specified file.\n";
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
	char procpar[MAXSTR], *prefix;
	par_t *pars;
	double n, nomFlip;
	double ***tr1, ***tr2, ***flip, ***B1, 
	       ***smoothB1, ***mask = NULL;
	short nx, ny, nz, nvol;
	FSLIO *in = NULL, *out = NULL, *maskHdr = NULL;
	while ((c = getopt_long(argc, argv, "m:", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 'm':
				maskHdr = FslOpen(optarg, "rb");
				fprintf(stdout, "Reading mask.\n");
				mask = FslGetVolumeAsScaledDouble(maskHdr, 0);
				FslClose(maskHdr);
				break;
		}
	}
	if ((argc - optind) != 2)
	{
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	fprintf(stdout, "Opening input file %s...\n", argv[optind]);
	in = FslOpen(argv[optind], "rb");
	strncpy(procpar, argv[optind], MAXSTR);
	strcat(procpar, ".procpar");
	if ((pars = readProcpar(procpar)))
	{
		n = realVal(pars, "afi_n", 0);
		nomFlip = realVal(pars, "flip1", 0);
		fprintf(stdout, "Read TR2/TR1 ratio of %f and flip-angle %f degrees from procpar.\n", n, nomFlip);
	}
	else
	{
		fprintf(stdout, "Enter TR2/TR1 (ratio) and flip-angle (degrees): ");
		fscanf(stdin, "%lf %lf", &n, &nomFlip);
	}
	nomFlip = radians(nomFlip);
	FslGetDim(in, &nx, &ny, &nz, &nvol);
	tr1 = FslGetVolumeAsScaledDouble(in, 0);
	tr2 = FslGetVolumeAsScaledDouble(in, 1);
	prefix = argv[++optind];
	fprintf(stdout, "Image dimensions: %d %d %d\n", nx, ny, nz);
	flip  = d3matrix(nz - 1, ny - 1, nx - 1);
	B1    = d3matrix(nz - 1, ny - 1, nx - 1);
	smoothB1    = d3matrix(nz - 1, ny - 1, nx - 1);
	fprintf(stdout, "Allocated output memory.\n");
	fprintf(stdout, "Processing...");
	for (size_t z = 0; z < nz; z++)
	{	for (size_t y = 0; y < ny; y++)
		{	for (size_t x = 0; x < nx; x++)
			{
				if (!mask || mask[z][y][x] > 0.)
				{
					double r = tr2[z][y][x] / tr1[z][y][x];
					double temp = (r*n - 1.) / (n - r);
					if (temp > 1.)
						temp = 1.;
					if (temp < -1.)
						temp = -1.;
					double alpha = acos(temp);
					flip[z][y][x] = degrees(alpha);
					B1[z][y][x]   = alpha / nomFlip;
				}
			}
		}
	}
	fprintf(stdout, "done.\n");
	char outfile[1024];
	snprintf(outfile, 1024, "%s_flip.nii.gz", prefix);
	fprintf(stdout, "Writing to %s...\n", outfile);
	out = FslOpen(outfile, "wb");
	FslCloneHeader(out, in);
	FslSetDim(out, nx, ny, nz, 1);
	FslSetDataType(out, NIFTI_TYPE_FLOAT32);
	FslWriteHeader(out);
	FslWriteVolumeFromDouble(out, flip, 0);
	fprintf(stdout, "Wrote flip angle.\n");
	FslClose(out);
	snprintf(outfile, 1024, "%s_B1.nii.gz", prefix);
	out = FslOpen(outfile, "wb");
	FslCloneHeader(out, in);
	FslSetDim(out, nx, ny, nz, 1);
	FslSetDataType(out, NIFTI_TYPE_FLOAT32);
	FslWriteHeader(out);
	FslWriteVolumeFromDouble(out, B1, 0);
	fprintf(stdout, "Wrote B1 ratio.\n");
	FslClose(out);
	if (smooth)
	{
		fprintf(stdout, "Smoothing...");
		double *gaussKernel = gaussian3D(7, 7, 7, 1.5, 1.5, 1.5);
		convolve3D(smoothB1[0][0], B1[0][0], nx, ny, nz,
		           gaussKernel, 7, 7, 7);
		fprintf(stdout, "done.\n");
		snprintf(outfile, 1024, "%s_smooth_B1.nii.gz", prefix);
		out = FslOpen(outfile, "wb");
		FslCloneHeader(out, in);
		FslSetDim(out, nx, ny, nz, 1);
		FslSetDataType(out, NIFTI_TYPE_FLOAT32);
		FslWriteHeader(out);
		FslWriteVolumeFromDouble(out, smoothB1, 0);
		fprintf(stdout, "Wrote smoothed B1 ratio.\n");
		FslClose(out);
	}
	FslClose(in);
	fprintf(stdout, "Success.\n");
    return EXIT_SUCCESS;
}

