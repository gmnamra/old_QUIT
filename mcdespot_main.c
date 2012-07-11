/*
 *  mcdespot_main.c
 *  Fitting
 *
 *  Created by Tobias Wood on 14/02/2012.
 *  Copyright 2012 Tobias Wood. All rights reserved.
 *
 */

#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <dispatch/dispatch.h>
#include <getopt.h>
#include <signal.h>

#ifdef __APPLE__
	#include <libkern/OSAtomic.h>
	#define AtomicAdd OSAtomicAdd32
#else
	#define AtomicAdd(x, y) (*y)++
#endif

#include "DESPOT.h"
#include "fslio.h"
#include "mathsArray.h"

//******************************************************************************
// Constants
//******************************************************************************
char *usage = "Usage is: mcdespot [options] output_prefix spgr_file ssfp_file1 (ssfp_fileN)\n\
\
Options:\n\
	-m file   : Mask input with specified file.\n\
	-z        : Output .nii.gz files.\n\
	--M0 file : M0 Map file.\n\
	--B0 file : B0 Map file.\n\
	--B1 file : B1 Map file.\n\
	-b 3      : Boundaries suitable for 3T (default)\n\
	   7      : Boundaries suitable for 7T\n\
	   u      : User specified boundaries from stdin.\n";

const double lo3Bounds[6] = { 0.200, 0.500, 0.002, 0.040, 0.0, 0.050 },
		     hi3Bounds[6] = { 0.700, 2.500, 0.040, 0.200, 0.4, 2.000 },
		     lo7Bounds[6] = { 0.200, 0.500, 0.002, 0.040, 0.,   0.050 },
		     hi7Bounds[6] = { 1.500, 5.000, 0.050, 0.400, 0.5,  2.000 };
//******************************************************************************
// SIGTERM interrupt handler - for ensuring data gets saved even on a ctrl-c
//******************************************************************************
#define NR 7
FSLIO *resultsHeaders[NR] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL };
float *resultsData[NR];
const char *names[NR] = { "_T1_myel", "_T1_free", "_T2_myel", "_T2_free", "_frac_myel", "_tau_myel", "_res" };
void int_handler(int sig);
void int_handler(int sig)
{
	fprintf(stdout, "Processing terminated. Writing currently processed data.\n");
	for (size_t r = 0; r < NR; r++)
	{
		FslWriteHeader(resultsHeaders[r]);
		FslWriteVolumes(resultsHeaders[r], resultsData[r], 1);
		FslClose(resultsHeaders[r]);
	}
	exit(EXIT_FAILURE);
}


//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	if (argc < 4)
	{
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	
	char *outPrefix = NULL, *outExt = ".nii";
	size_t nSPGR, nPhases, nSSFP;
	double spgrTR, *spgrAngles,
		   ssfpTR, loBounds[NR - 1], hiBounds[NR - 1];
	double **bounds = malloc(2 * sizeof(double *));
	bounds[0] = loBounds; bounds[1] = hiBounds;
	bool loConstraint[7] = { true, true, true, true, true, true };
	bool hiConstraint[7] = { true, true, true, true, true, true };
	bool **constraints = malloc(2 * sizeof(double *));
	constraints[0] = loConstraint; constraints[1] = hiConstraint;
	FSLIO *maskFile = NULL, *M0File = NULL, *B0File = NULL, *B1File = NULL,
	      *spgrFile = NULL;
	
	static struct option long_options[] =
	{
		{"B0", required_argument, 0, '0'},
		{"B1", required_argument, 0, '1'},
		{"M0", required_argument, 0, 'M'},
		{0, 0, 0, 0}
	};
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "b:m:z", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 'm':
				if (!(maskFile = FslOpen(optarg, "rb")))
					exit(EXIT_FAILURE);
				break;
			case '0':
				if (!(B0File = FslOpen(optarg, "rb")))
					exit(EXIT_FAILURE);
				break;
			case '1':
				if (!(B1File = FslOpen(optarg, "rb")))
					exit(EXIT_FAILURE);
				break;
			case 'M':
				if (!(M0File = FslOpen(optarg, "rb")))
					exit(EXIT_FAILURE);
				break;
			case 'z':
				outExt = ".nii.gz";
				break;
			case 'b':
				switch (*optarg)
				{
					case '3':
						fprintf(stdout, "Using 3T boundaries.\n");
						arrayCopy(loBounds, lo3Bounds, NR - 1);
						arrayCopy(hiBounds, hi3Bounds, NR - 1);
						break;
					case '7':
						fprintf(stdout, "Using 7T boundaries.\n");
						arrayCopy(loBounds, lo7Bounds, NR - 1);
						arrayCopy(hiBounds, hi7Bounds, NR - 1);
						break;
					case 'u':
						fprintf(stdout, "Enter low boundaries (%d values):", NR - 1);
						fgetArray(stdin, 'd', NR - 1, loBounds);
						fprintf(stdout, "Enter high boundaries (%d values):", NR - 1);
						fgetArray(stdin, 'd', NR - 1, hiBounds);
						break;
					default:
						fprintf(stdout, "Unknown boundaries type '%s'.\n", optarg);
						exit(EXIT_FAILURE);
						break;
				}
				break;
		}
	}
	outPrefix = argv[optind++];
	fprintf(stdout, "Output prefix will be: %s\n", outPrefix);
	//**************************************************************************
	// Gather SPGR Data
	//**************************************************************************
	fprintf(stdout, "Reading SPGR header from %s.\n", argv[optind]);
	if (!(spgrFile = FslOpen(argv[optind++], "rb")))
		exit(EXIT_FAILURE);
	nSPGR = spgrFile->niftiptr->nt;
	spgrAngles = malloc(nSPGR * sizeof(double));
	fprintf(stdout, "Enter SPGR TR (seconds):");
	fscanf(stdin, "%f", &spgrTR);
	fprintf(stdout, "\nEnter %zu SPGR Flip Angles (degrees):", nSPGR);
	fgetArray(stdin, 'd', nSPGR, spgrAngles);
	arrayApply(spgrAngles, spgrAngles, radians, nSPGR);

	//**************************************************************************
	// Gather SSFP Data
	//**************************************************************************
	nPhases    = argc - optind;
	double *ssfpPhases = malloc(nPhases * sizeof(double));
	FSLIO *ssfpFiles[nPhases];
	fprintf(stdout, "Reading first SSFP header %s\n", argv[optind]);
	ssfpFiles[0] = FslOpen(argv[optind++], "rb");
	if (!FslCheckDims(spgrFile, ssfpFiles[0]))
	{
		fprintf(stderr, "Differing number of voxels in SPGR and SSFP data.\n");
		exit(EXIT_FAILURE);
	}
	fprintf(stdout, "Specified %zu phase cycling patterns.\n", nPhases);
	fprintf(stdout, "Enter SSFP TR (seconds):");
	fscanf(stdin, "%f", &ssfpTR);
	fprintf(stdout, "\nEnter %zu Phase-Cycling Patterns (degrees):", nPhases);
	fgetArray(stdin, 'd', nPhases, ssfpPhases); fprintf(stdout, "\n");
	arrayApply(ssfpPhases, ssfpPhases, radians, nPhases);
	
	nSSFP = ssfpFiles[0]->niftiptr->nt;
	double ssfpAngles[nSSFP];
	fprintf(stdout, "Enter %zu SSFP Flip Angles (degrees):", nSSFP);
	fgetArray(stdin, 'd', nSSFP, ssfpAngles); fprintf(stdout, "\n");
	arrayApply(ssfpAngles, ssfpAngles, radians, nSSFP);
	for (int p = 1; p < nPhases; p++)
	{
		fprintf(stdout, "Reading %f SSFP header from %s.\n", degrees(ssfpPhases[p]), argv[optind]);
		ssfpFiles[p] = FslOpen(argv[optind++], "rb");
		if (!FslCheckDims(spgrFile, ssfpFiles[p]))
		{
			fprintf(stderr, "Differing number of voxels in phase %d and %d headers.\n", 0, p);
			exit(EXIT_FAILURE);
		}
		if (ssfpFiles[p]->niftiptr->nt != nSSFP)
		{
			fprintf(stderr, "Wrong number of flip angles in SSFP file %s.\n", argv[optind - 1]);
			exit(EXIT_FAILURE);
		}
	}
	fprintf(stdout, "Read all SSFP files.\n");
	//**************************************************************************	
	// Allocate memory for slices
	//**************************************************************************
	int nx = spgrFile->niftiptr->nx;
	int ny = spgrFile->niftiptr->ny;
	int nz = spgrFile->niftiptr->nz;
	int voxelsPerSlice = nx * ny;
	int totalVoxels = voxelsPerSlice * nz;
	float *spgrData, **ssfpData, *maskData, *M0Data, *B0Data, *B1Data;
	spgrData = malloc(nSPGR * voxelsPerSlice * sizeof(float));
	ssfpData = malloc(nPhases * sizeof(float *));
	for (int p = 0; p < nPhases; p++)
		ssfpData[p] = malloc(nSSFP * voxelsPerSlice * sizeof(float));
	if (maskFile)
		maskData = malloc(voxelsPerSlice * sizeof(float));
	if (M0File)
		M0Data = malloc(voxelsPerSlice * sizeof(float));
	if (B0File)
		B0Data = malloc(voxelsPerSlice * sizeof(float));
	if (B1File)
		B1Data = malloc(voxelsPerSlice * sizeof(float));
	fprintf(stdout, "Allocated memory.\n");
	//**************************************************************************
	// Create results files
	// T1_m, T1_m, T2_f, T2_f,	f_m, tau_m, residue
	//**************************************************************************
	char outName[strlen(outPrefix) + 32];
	for (int r = 0; r < NR; r++)
	{
		strcpy(outName, outPrefix); strcat(outName, names[r]); strcat(outName, outExt);
		fprintf(stdout, "Creating result header %s\n", outName);
		resultsHeaders[r] = FslOpen(outName, "wb");
		FslCloneHeader(resultsHeaders[r], spgrFile);
		FslSetDim(resultsHeaders[r], nx, ny, nz, 1);
		FslSetDataType(resultsHeaders[r], DTYPE_FLOAT);
		resultsData[r] = malloc(totalVoxels * sizeof(float));
	}
	signal(SIGINT, int_handler);
	//**************************************************************************
	// Set up functions and angles together for regionContraction
	//**************************************************************************
	size_t *nD = malloc((1 + nPhases) * sizeof(size_t));
	eval_array_type **f = malloc((1 + nPhases) * sizeof(eval_array_type *));
	double **angles = malloc((1 + nPhases) * sizeof(double *));
	
	nD[0] = nSPGR;
	f[0]      = a2cSPGR;
	angles[0] = spgrAngles;
	for (int p = 0; p < nPhases; p++)
	{
		nD[p + 1]     = nSSFP;
		f[p + 1]      = a2cSSFP;
		angles[p + 1] = ssfpAngles;
	}
	
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	time_t allStart = time(NULL);
	struct tm *localStart = localtime(&allStart);
	char theTime[1024];
	strftime(theTime, 1024, "%H:%M:%S", localStart);
	fprintf(stdout, "Started processing at %s.\n", theTime);
	for (size_t slice = 0; slice < nz; slice++)
	{
		// Read in data
		fprintf(stdout, "Starting slice %zu...\n", slice);
		
		FslReadSliceSeries(spgrFile, (void *)spgrData, slice, nSPGR);
		if (B0File)
			FslReadSliceSeries(B0File, (void *)B0Data, slice, 1);
		if (B1File)
			FslReadSliceSeries(B1File, (void *)B1Data, slice, 1);
		if (M0File)
			FslReadSliceSeries(M0File, (void *)M0Data, slice, 1);
		if (maskFile)
			FslReadSliceSeries(maskFile, (void *)maskData, slice, 1);
		for (int p = 0; p < nPhases; p++)
			FslReadSliceSeries(ssfpFiles[p], (void *)ssfpData[p], slice, nSSFP);
		fprintf(stdout, "read data.\n");
		__block int voxCount = 0;
		int sliceOffset = slice * voxelsPerSlice;
		clock_t loopStart = clock();
		void (^processVoxel)(size_t vox) = ^(size_t vox)
		{
			double params[NR];
			arraySet(params, 0., NR);
			if (!maskFile || (maskData[vox] > 0.))
			{
				AtomicAdd(1, &voxCount);
				
				double M0 = 1, B0 = 0, B1 = 1.;
				if (M0File)
					M0 = (double)M0Data[vox];
				if (B0File)
					B0 = (double)B0Data[vox];
				if (B1File)
					B1 = (double)B1Data[vox];
				
				// Constants need to be set up here because B1 changes per-voxels
				double *signals[1 + nPhases], *consts[1 + nPhases];
				
				signals[0] = malloc(nSPGR * sizeof(double));
				for (int img = 0; img < nSPGR; img++)
					signals[0][img] = (double)spgrData[voxelsPerSlice * img + vox];
				consts[0] = malloc(2 * sizeof(double));
				consts[0][0] = spgrTR; consts[0][1] = B1;
				
				for (int p = 0; p < nPhases; p++)
				{
					signals[p + 1] = malloc(nSSFP * sizeof(double));
					for (int img = 0; img < nSSFP; img++)
						signals[p + 1][img] = (double)ssfpData[p][voxelsPerSlice * img + vox];
					
					consts[p + 1] = malloc(5 * sizeof(double));
					consts[p + 1][0] = ssfpTR; consts[p + 1][1] = M0;
					consts[p + 1][2] = B0;     consts[p + 1][3] = B1;
					consts[p + 1][4] = ssfpPhases[p];
				}

				// Store residual in final parameter
				regionContraction(params, NR - 1, consts, 1 + nPhases, angles,
				                  signals, nD, false, f, bounds, constraints,
								  1000, 10, 20, 0.1, 0.025, &(params[NR - 1]));		
				// Clean up memory
				for (int p = 0; p < 1 + nPhases; p++)
				{
					free(signals[p]);
					free(consts[p]);
				}
			}
			for (int p = 0; p < NR; p++)
				resultsData[p][sliceOffset + vox]  = (float)params[p];
		};
		dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
		dispatch_apply(voxelsPerSlice, global_queue, processVoxel);
		
        clock_t loopEnd = clock();
        fprintf(stdout, "Finished slice %zu", slice);
		if (voxCount > 0)
			fprintf(stdout, ", had %d unmasked voxels, CPU time per voxel was %f s.",
				    voxCount, (loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC));
		fprintf(stdout, ".\n");
	}
    time_t allEnd = time(NULL);
    struct tm *localEnd = localtime(&allEnd);
    strftime(theTime, 1024, "%H:%M:%S", localEnd);
	fprintf(stdout, "Finished processing at %s. Run-time was %f s.\n", theTime, difftime(allEnd, allStart));

	for (size_t r = 0; r < NR; r++)
	{
		FslWriteHeader(resultsHeaders[r]);
		FslWriteVolumes(resultsHeaders[r], resultsData[r], 1);
		FslClose(resultsHeaders[r]);
	}
	// Clean up memory
	free(spgrData);
	for (int p = 0; p < nPhases; p++)
		free(ssfpData[p]);
	free(M0Data);
	free(B0Data);
	free(B1Data);
	free(maskData);
	exit(EXIT_SUCCESS);
}

