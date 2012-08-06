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
#include "mathsOps.h"
#include "mathsArray.h"

char *usage = "Usage is: phasemap input_1 input_2 outprefix\n\
\n\
Echo times will be read from procpar if present.\n\
Options:\n\
	--mask, -m mask_file : Mask input with specified file\n\
	--phasetime T        : Calculate the phase accumulated in time T\n\
	--smooth             : Smooth output with a gaussian.\n";

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
	char procpar[MAXSTR];
	par_t *pars;
	double TE1, TE2, deltaTE, phasetime = 0.;
	double *data1, *data2, *B0, *mask = NULL;
	short nx, ny, nz, nvol;
	FSLIO *in1 = NULL, *in2 = NULL, *out = NULL, *maskHdr = NULL;
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
			case 'p':
				phasetime = atof(optarg);
				break;
		}
	}
	if ((argc - optind) == 2)
	{
		fprintf(stdout, "Opening input file %s...\n", argv[optind]);
		in1 = FslOpen(argv[optind], "rb");
		strncpy(procpar, argv[optind], MAXSTR);
		strcat(procpar, ".procpar");
		if ((pars = readProcpar(procpar)))
		{
			TE1 = realVal(pars, "te", 0);
			TE2 = realVal(pars, "te", 1);
		}
		else
		{
			fprintf(stdout, "Enter TE2 & TE2 (seconds): ");
			fscanf(stdin, "%lf %lf", &TE1, &TE2);
		}
		FslGetDim(in1, &nx, &ny, &nz, &nvol);
		data1 = FslGetVolumeAsScaledDouble(in1, 0);
		data2 = FslGetVolumeAsScaledDouble(in1, 1);
	}
	else if ((argc - optind) == 3)
	{
		fprintf(stdout, "Opening input file 1 %s...\n", argv[optind]);
		in1 = FslOpen(argv[optind], "rb");
		strncpy(procpar, argv[optind], MAXSTR);
		strcat(procpar, ".procpar");
		if ((pars = readProcpar(procpar)))
			TE1 = realVal(pars, "te", 0);
		else {
			fprintf(stdout, "Enter TE1 (seconds): ");
			fscanf(stdin, "%lf", &TE1);
		}
		fprintf(stdout, "Opening input file 2 %s...\n", argv[++optind]);
		in2 = FslOpen(argv[optind], "rb");
		strncpy(procpar, argv[optind], MAXSTR);
		strcat(procpar, ".procpar");
		if ((pars = readProcpar(procpar)))
			TE2 = realVal(pars, "te", 0);
		else
		{
			fprintf(stdout, "Enter TE2 (seconds): ");
			fscanf(stdin, "%lf", &TE2);
		}
		FslGetDim(in1, &nx, &ny, &nz, &nvol);
		int nvox = nx * ny * nz;
		FslGetDim(in2, &nx, &ny, &nz, &nvol);
		if (nvox != nx * ny * nz)
		{
			fprintf(stderr, "File dimensions do not match.\n");
			exit(EXIT_FAILURE);
		}
		data1 = FslGetVolumeAsScaledDouble(in1, 0);
		data2 = FslGetVolumeAsScaledDouble(in2, 0);
	}
	else
	{
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	char *prefix = argv[++optind];
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
	fprintf(stdout, "Image dimensions: %d %d %d\n", nx, ny, nz);
	B0    = malloc(nz*ny*nx * sizeof(double));
	fprintf(stdout, "Allocated output memory.\n");
	fprintf(stdout, "Processing...");
	for (size_t vox = 0; vox < nz*ny*nx; vox++)
	{	
		if (!mask || mask[vox] > 0.)
		{
			double deltaPhase = data2[vox] - data1[vox];
			B0[vox] = deltaPhase / (2 * M_PI * deltaTE);
			if (phasetime > 0.)
			{
				double ph = fmod(B0[vox] * 2 * M_PI * phasetime, 2 * M_PI);
				if (ph > M_PI) ph -= (2 * M_PI);
				if (ph < -M_PI) ph += (2 * M_PI);
				B0[vox] = ph;
			}
		}
	}
	fprintf(stdout, "done.\n");
	char outfile[1024];
	snprintf(outfile, 1024, "%s_B0.nii.gz", prefix);
	out = FslOpen(outfile, "wb");
	FslCloneHeader(out, in1);
	FslSetDim(out, nx, ny, nz, 1);
	FslSetDataType(out, NIFTI_TYPE_FLOAT32);
	FslWriteHeader(out);
	FslWriteVolumeFromDouble(out, B0, 0);
	fprintf(stdout, "Wrote B0 map.\n");
	FslClose(out);
	if (smooth)
	{
		fprintf(stdout, "Smoothing...");
		array3d_t *B0_3d = array3d_from_buffer(B0, nz, ny, nx);
		array3d_t *gauss = gaussian3D(7, 7, 7, 1.5, 1.5, 1.5);
		array3d_t *smoothed = array3d_alloc(nz, ny, nx);
		convolve3D(smoothed, B0_3d, gauss);
		arrayMul(smoothed->array->data, smoothed->array->data, mask, nx * ny * nz);
		fprintf(stdout, "done.\n");
		snprintf(outfile, 1024, "%s_B0_smooth.nii.gz", prefix);
		out = FslOpen(outfile, "wb");
		FslCloneHeader(out, in1);
		FslSetDim(out, nx, ny, nz, 1);
		FslSetDataType(out, NIFTI_TYPE_FLOAT32);
		FslWriteHeader(out);
		FslWriteVolumeFromDouble(out, smoothed->array->data, 0);
		fprintf(stdout, "Wrote smoothed B0.\n");
		FslClose(out);	
	}
	FslClose(in1);
	if (in2)
		FslClose(in2);
	fprintf(stdout, "Success.\n");
    return EXIT_SUCCESS;
}

