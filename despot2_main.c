/*
 *  despot1_main.c
 *  MacRI
 *
 *  Created by Tobias Wood on 23/01/2012.
 *  Copyright 2012 Tobias Wood. All rights reserved.
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
#include "mathsArray.h"
#include "procpar.h"

char *usage = "Usage is: despot2 [options] output_prefix T1_map ssfp_180_file [additional ssfp files] \n\
\
Options:\n\
	-m, --mask file   : Mask input with specified file.\n\
	--B0 file         : B0 Map File.\n\
	--B1 file         : B1 Map File.\n\
	--lm              : Use Levenberg-Marquardt instead of Region Contraction.\n";
//******************************************************************************
// SIGINT interrupt handler - for ensuring data gets saved even on a ctrl-c
//******************************************************************************
#define NR 4
FSLIO *resultsHeaders[NR];
double *resultsData[NR];
char *names[NR] = { "_d2_M0", "_T2", "_d2_B0", "_d2_res" };
void int_handler(int sig);
void int_handler(int sig)
{
	fprintf(stdout, "Processing terminated. Writing currently processed data.\n");
	for (size_t r = 0; r < NR; r++)
	{
		FslWriteHeader(resultsHeaders[r]);
		FslWriteVolumeFromDouble(resultsHeaders[r], resultsData[r], 0);
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

	char *outPrefix = NULL, *outExt = ".nii.gz";
	size_t nPhases;
	double ssfpTR, *ssfpPhases = NULL, *ssfpAngles = NULL;
	FSLIO **ssfpFiles = NULL, *inFile = NULL;
	double *maskData = NULL, *B0Data = NULL, *B1Data = NULL, *T1Data = NULL;
	char procpar[MAXSTR];
	par_t *pars;
	
	static int levMar = false;
	static struct option long_options[] =
	{
		{"B0", required_argument, 0, '0'},
		{"B1", required_argument, 0, '1'},
		{"mask", required_argument, 0, 'm'},
		{"lm", no_argument, &levMar, true},
		{0, 0, 0, 0}
	};
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "m:z", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 'm':
				inFile = FslOpen(optarg, "rb");
				maskData = FslGetVolumeAsScaledDouble(inFile, 0);
				FslClose(inFile);
				break;
			case '0':
				inFile = FslOpen(optarg, "rb");
				B0Data = FslGetVolumeAsScaledDouble(inFile, 0);
				FslClose(inFile);
				break;
			case '1':
				inFile = FslOpen(optarg, "rb");
				B1Data = FslGetVolumeAsScaledDouble(inFile, 0);
				FslClose(inFile);
				break;
			case 0:
				// Just a flag
				break;
		}
	}
	
	fprintf(stdout, "Output prefix will be: %s\n", argv[optind]);
	outPrefix = argv[optind++];
	fprintf(stdout, "Reading T1 Map from: %s\n", argv[optind]);
	inFile = FslOpen(argv[optind++], "rb");
	T1Data = FslGetVolumeAsScaledDouble(inFile, 0);
	FslClose(inFile);
	//**************************************************************************
	// Gather SSFP Data
	//**************************************************************************
	nPhases = argc - optind;
	if (nPhases < 1)
	{
		fprintf(stderr, "Must have at least the 180 degree phase-cycling pattern to process.\n");
		exit(EXIT_FAILURE);
	}
	ssfpFiles  = malloc(nPhases * sizeof(nifti_image *));
	ssfpPhases = malloc(nPhases * sizeof(double));
	ssfpPhases[0] = M_PI;
	ssfpFiles[0] = FslOpen(argv[optind], "rb");
	size_t nx = ssfpFiles[0]->niftiptr->nx;
	size_t ny = ssfpFiles[0]->niftiptr->ny;
	size_t nz = ssfpFiles[0]->niftiptr->nz;
	size_t nSSFP = ssfpFiles[0]->niftiptr->nt;
	ssfpAngles = malloc(nSSFP * sizeof(double));	
	snprintf(procpar, MAXSTR, "%s.procpar", argv[optind]);
	pars = readProcpar(procpar);
	if (pars)
	{
		int check;
		fprintf(stdout, "Reading SSFP 180 parameters from procpar.\n");
		arrayCopy(ssfpAngles, realVals(pars, "flip1", &check), nSSFP);		
		if (check != nSSFP)
		{
			fprintf(stderr, "flip1 and nvols do not match.\n");
			exit(EXIT_FAILURE);
		}
		ssfpTR = realVal(pars, "tr", 0);
		fprintf(stdout, "TR (s): %f, Angles (deg) = ", ssfpTR);
		arrayPrint(stdout, ssfpAngles, nSSFP); fprintf(stdout, "\n");
		arrayApply(ssfpAngles, ssfpAngles, radians, nSSFP);
		freeProcpar(pars);
	}
	else
	{
		fprintf(stdout, "Enter %zu SSFP flip angles (degrees) :", nSSFP);
		fgetArray(stdin, 'd', nSSFP, ssfpAngles); fprintf(stdout, "\n");
		arrayApply(ssfpAngles, ssfpAngles, radians, nSSFP);
		fprintf(stdout, "Enter SSFP TR (ms):");
		fscanf(stdin, "%lf", &ssfpTR);		
	}

	for (size_t p = 1; p < nPhases; p++)
	{
		fprintf(stdout, "Reading SSFP header from %s.\n", argv[++optind]);
		ssfpFiles[p] = FslOpen(argv[optind], "rb");
		if (!FslCheckDims(ssfpFiles[0], ssfpFiles[p]))
		{
			fprintf(stderr, "Image %s has differing dimensions.\n", argv[optind]);
			exit(EXIT_FAILURE);
		}
		if (ssfpFiles[p]->niftiptr->nt != nSSFP)
		{
			fprintf(stderr, "Image %s has wrong number of flip angles.\n", argv[optind]);
			exit(EXIT_FAILURE);
		}
		snprintf(procpar, MAXSTR, "%s.procpar", argv[optind]);
		pars = readProcpar(procpar);
		if (pars)
			ssfpPhases[p] = radians(realVal(pars, "rfphase", 0));
		else
		{
			fprintf(stdout, "Enter phase-cycling (degrees):");
			fscanf(stdin, "%lf", ssfpPhases + p);
			ssfpPhases[p] = radians(ssfpPhases[p]);	
		}
	}
	//**************************************************************************	
	// Get input data
	//**************************************************************************
	int voxelsPerSlice = nx * ny;
	int totalVoxels = voxelsPerSlice * nz;
	double **ssfpData = malloc(nPhases * sizeof(double *));
	for (int p = 0; p < nPhases; p++)
		ssfpData[p] = FslGetAllVolumesAsScaledDouble(ssfpFiles[p]);
	fprintf(stdout, "Read SSFP data.\n");
	//**************************************************************************
	// Create results files
	// T2, residue
	//**************************************************************************
	char outName[strlen(outPrefix) + 15];
	for (int r = 0; r < NR; r++)
	{
		strcpy(outName, outPrefix); strcat(outName, names[r]); strcat(outName, outExt);
		fprintf(stdout, "Writing result header:%s.\n", outName);
		resultsHeaders[r] = FslOpen(outName, "wb");
		FslCloneHeader(resultsHeaders[r], ssfpFiles[0]);
		FslSetDim(resultsHeaders[r], nx, ny, nz, 1);
		FslSetDataType(resultsHeaders[r], DTYPE_FLOAT);
		resultsData[r] = malloc(totalVoxels * sizeof(double));
	}
	// Now register the SIGINT handler so we can ctrl-c and still get data
	signal(SIGINT, int_handler);
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
		__block int voxCount = 0;
		clock_t loopStart = clock();
		int sliceOffset = slice * voxelsPerSlice;
		//for (size_t vox = 0; vox < voxelsPerSlice; vox++)
		void (^processVoxel)(size_t vox) = ^(size_t vox)
		{
			double params[NR];
			arraySet(params, 0., NR);
			if (!maskData || ((maskData[sliceOffset + vox] > 0.) && (T1Data[sliceOffset + vox] > 0.)))
			{	// Zero T1 causes zero-pivot error.
				AtomicAdd(1, &voxCount);
				
				// Gather signals. 180 Phase data is required to be the first file passed into program
				double *signals[nPhases];
				for (int p = 0; p < nPhases; p++)
				{
					signals[p] = malloc(nSSFP * sizeof(double));
					for (int img = 0; img < nSSFP; img++)
						signals[p][img] = (double)ssfpData[p][img * totalVoxels + sliceOffset + vox];
				}
				// Some constants set up here because they change per-voxel
				double T1 = (double)T1Data[sliceOffset + vox];
				double B0 = 0, B1 = 1.;
				if (B0Data)
					B0 = (double)B0Data[sliceOffset + vox];
				if (B1Data)
					B1 = (double)B1Data[sliceOffset + vox];
				// Run classic DESPOT2 on 180 phase data
				if ((nPhases == 1) || levMar)
				{
					params[3] = classicDESPOT2(ssfpAngles, signals[0], nSSFP, ssfpTR, T1, B1, params);
					params[0] = clamp(params[0], 0, 1.e7);
					params[1] = clamp(params[1], 0.001, 0.250);
					params[2] = 0.;
				}
                                // Need to still trigger this if doing levMar, hence if not else if
				if (nPhases > 1)
				{
					double *xData[nPhases], *consts[nPhases];
					size_t dSize[nPhases];
					eval_array_type *fs[nPhases];
					int nP;
					
					for (int p = 0; p < nPhases; p++)
					{
						xData[p] = ssfpAngles;
						dSize[p] = nSSFP;
						consts[p] = arrayAlloc(5);
						consts[p][0] = ssfpTR;
						consts[p][1] = T1;
						consts[p][2] = B0;
						params[2]    = B0;
						consts[p][3] = B1;
						consts[p][4] = ssfpPhases[p];						
						if (levMar)
						{
						    nP = 3;
						    fs[p] = &a1cSSFPB0;
						}
						else
						{
						    nP = 3;
						    fs[p] = &a1cSSFPB0;
						}
					}
					
					double loBounds[3] = { 1.e5, 0.005, -.5 / ssfpTR };
					double hiBounds[3] = { 1.e6, 0.100,  .5 / ssfpTR };
					if (B0Data)
					{
						loBounds[2] = B0 * 0.95;
						hiBounds[2] = B0 * 1.05;
					}
					double *bounds[2] = { loBounds, hiBounds };
					bool loC[3] = { FALSE, TRUE, FALSE }, hiC[3] = { FALSE, FALSE, FALSE };
					bool *constrained[2] = { loC, hiC };
					if (levMar)
						levenbergMarquardt(params, nP, consts, xData, signals, fs, NULL, dSize, nPhases, NULL, NULL, false, params + 3);
					else
						regionContraction(params, nP, consts, nPhases, xData, signals, dSize, false, fs, bounds, constrained, 10000, 20, 10, 0.05, 0.05, params + 3);
					params[3] = fmod(params[2] * 2 * M_PI, 2 * M_PI);
					if (params[3] > M_PI) params[3] -= (2 * M_PI);
					if (params[3] < -M_PI) params[3] += (2 * M_PI);
				}
				// Clean up memory
				for (int p = 0; p < nPhases; p++)
					free(signals[p]);
			}
			for (int p = 0; p < NR; p++)
				resultsData[p][sliceOffset + vox]  = params[p];
		};
		dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
		dispatch_apply(voxelsPerSlice, global_queue, processVoxel);
		
        clock_t loopEnd = clock();
        fprintf(stdout, "Finished slice %zu", slice);
		if (voxCount > 0)
		{
			fprintf(stdout, ", had %d unmasked voxels, CPU time per voxel was %f s", 
			        voxCount, (loopEnd - loopStart) / ((double)voxCount * CLOCKS_PER_SEC));
		}
		fprintf(stdout, ".\n");
	}
    time_t allEnd = time(NULL);
    struct tm *localEnd = localtime(&allEnd);
    strftime(theTime, 1024, "%H:%M:%S", localEnd);
	fprintf(stdout, "Finished processing at %s. Run-time was %f s.\n", theTime, difftime(allEnd, allStart));
	
	for (size_t r = 0; r < NR; r++)
	{
		FslWriteHeader(resultsHeaders[r]);
		FslWriteVolumeFromDouble(resultsHeaders[r], resultsData[r], 0);
		FslClose(resultsHeaders[r]);
	}
	// Clean up memory
	for (int p = 0; p < nPhases; p++)
		free(ssfpData[p]);
	if (B0Data)
		free(B0Data);
	if (B1Data)
		free(B1Data);
	if (maskData)
		free(maskData);
	exit(EXIT_SUCCESS);
}
