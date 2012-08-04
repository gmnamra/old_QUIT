/*
 *  despot1_main.c
 *  MacRI
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright 2011 Tobias Wood. All rights reserved.
 *
 */

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
#include "DESPOT.h"
#include "fslio.h"
#include "procpar.h"
#include "mathsOps.h"
char *usage = "Usage is: despot1 [options] spgr_input output_prefix \n\
\
Options:\n\
	-m, --mask file : Mask input with specified file.\n\
	--B1 file       : Correct flip angles with specified B1 ratio.\n";
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//tests();
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	int nSPGR = 0;
	short nx, ny, nz, nv;
	char *outPrefix = NULL, *outExt = ".nii.gz", procpar[MAXSTR];
	double spgrTR = 0., *spgrAngles = NULL;
	FSLIO *spgrFile = NULL, *B1File = NULL, *maskFile = NULL;
	par_t *pars;
	
	static struct option long_options[] =
	{
		{"B1", required_argument, 0, '1'},
		{"mask", required_argument, 0, 'm'},
		{0, 0, 0, 0}
	};
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "m:n:z", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case '1':
				fprintf(stdout, "Opening B1 file: %s\n", optarg);
				B1File = FslOpen(optarg, "rb");
				break;
			case 'm':
				fprintf(stdout, "Opening mask file: %s\n", optarg);
				maskFile = FslOpen(optarg, "rb");
				break;
		}
	}
	if ((argc - optind) != 2)
	{
		fprintf(stderr, "Incorrect number of arguments.\n%s", usage);
		exit(EXIT_FAILURE);
	}
	fprintf(stdout, "Opening SPGR file: %s\n", argv[optind]);
	spgrFile = FslOpen(argv[optind], "rb");
	FslReadAllVolumes(spgrFile, argv[optind]);
	FslGetDim(spgrFile, &nx, &ny, &nz, &nv);
	nSPGR = nv;
	spgrAngles = malloc(nSPGR * sizeof(double));
	snprintf(procpar, MAXSTR, "%s.procpar", argv[optind]);
	pars = readProcpar(procpar);
	if (pars)
	{
		spgrTR = realVal(pars, "tr", 0);
		arrayCopy(spgrAngles, realVals(pars, "flip1", NULL), nSPGR);
		freeProcpar(pars);
	}
	else
	{
		fprintf(stdout, "Enter SPGR TR (s):");
		fscanf(stdin, "%lf", &spgrTR);
		fprintf(stdout, "Enter SPGR Flip Angles (degrees):");
		fgetArray(stdin, 'd', nSPGR, spgrAngles);
	}
	fprintf(stdout, "SPGR TR=%f s.\n", spgrTR);
	ARR_D(spgrAngles, nSPGR);
	arrayApply(spgrAngles, spgrAngles, radians, nSPGR);
		
	fprintf(stdout, "Ouput prefix will be: %s\n", argv[++optind]);
	outPrefix = argv[optind];
	//**************************************************************************	
	// Allocate memory for slices
	//**************************************************************************	
	int voxelsPerSlice = nx * ny;
	__block int voxCount;
	
	fprintf(stdout, "Reading SPGR data...\n");
	double ****SPGR = FslGetBufferAsScaledDouble(spgrFile);
	fprintf(stdout, "done.\n");
	double ***B1Data = NULL, ***maskData = NULL;
	if (B1File)
		B1Data = FslGetVolumeAsScaledDouble(B1File, 0);
	if (maskFile)
		maskData = FslGetVolumeAsScaledDouble(maskFile, 0);
	//**************************************************************************
	// Create results headers
	//**************************************************************************
	#define NR 3
	FSLIO **resultsHeaders = malloc(NR * sizeof(FSLIO *));
	double ****resultsData = d4matrix(NR, nz - 1, ny - 1, nx - 1);
	char *names[NR] = { "_M0", "_T1", "_despot1_res" };
	char outName[strlen(outPrefix) + 64];
	for (int r = 0; r < NR; r++)
	{
		strcpy(outName, outPrefix); strcat(outName, names[r]); strcat(outName, outExt);
		fprintf(stdout, "Writing result header:%s.\n", outName);
		resultsHeaders[r] = FslOpen(outName, "wb");
		FslCloneHeader(resultsHeaders[r], spgrFile);
		FslSetDim(resultsHeaders[r], nx, ny, nz, 1);
		FslSetDataType(resultsHeaders[r], DTYPE_FLOAT);
	}
		
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	for (int slice = 0; slice < nz; slice++)
	{
		clock_t loopStart, loopEnd;
		// Read in data
		fprintf(stdout, "Processing slice %d...\n", slice);
		loopStart = clock();
		voxCount = 0;
		//for (size_t vox = 0; vox < voxelsPerSlice; vox++)
		void (^processVoxel)(size_t vox) = ^(size_t vox)
		{
			double T1 = 0., M0 = 0., B1 = 1., res = 0.; // Place to restore per-voxel return values, assume B1 field is uniform for classic DESPOT
			if (B1File)
				B1 = B1Data[slice][0][vox];
			if ((!maskFile) || (maskData[slice][0][vox] > 0.))
			{
				AtomicAdd(1, &voxCount);
				double spgrs[nSPGR];
				for (int img = 0; img < nSPGR; img++)
					spgrs[img] = SPGR[img][slice][0][vox];
				//ARR_D(spgrs, nSPGR);
				res = calcDESPOT1(spgrAngles, spgrs, nSPGR, spgrTR, B1, &M0, &T1);
				
				// Sanity check
				M0 = clamp(M0, 0., 1.e7);
				T1 = clamp(T1, 0., 3.);
				//fprintf(stdout, " B1 %f M0 %f T1 %f res %f\n", B1, M0, T1, res);
			}
			
			resultsData[0][slice][0][vox] = M0;
			resultsData[1][slice][0][vox] = T1;
			resultsData[2][slice][0][vox] = res;
		};
		dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
		dispatch_apply(voxelsPerSlice, global_queue, processVoxel);
		
        loopEnd = clock();
        fprintf(stdout, "Finished slice %d", slice);
		if (voxCount)
		{
			fprintf(stdout, ", had %d unmasked voxels, CPU time per voxel was %f s.\n", 
			        voxCount, (loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC));
		}
		else
			fprintf(stdout, ", no unmasked voxels.\n");
	}
		
	for (size_t r = 0; r < NR; r++)
	{
		FslWriteHeader(resultsHeaders[r]);
		FslWriteVolumeFromDouble(resultsHeaders[r], resultsData[r], 0);
		FslClose(resultsHeaders[r]);
	}
	// Clean up memory
	free(SPGR);
	free(B1Data);
	free(maskData);
	fprintf(stdout, "All done.\n");
	exit(EXIT_SUCCESS);
}
