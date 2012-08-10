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
#include "procpar.h"

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

const double lo3Bounds[7] = { 0.200, 0.500, 0.002, 0.040, 0.0, 0.050, 0. },
		     hi3Bounds[7] = { 0.700, 2.500, 0.040, 0.200, 0.4, 2.000, 0. },
		     lo7Bounds[7] = { 0.200, 0.500, 0.002, 0.040, 0.,   0.050, 0. },
		     hi7Bounds[7] = { 1.500, 5.000, 0.050, 0.400, 0.5,  2.000, 0. };
//******************************************************************************
// SIGTERM interrupt handler - for ensuring data gets saved even on a ctrl-c
//******************************************************************************
#define NR 8
FSLIO *resultsHeaders[NR] = { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
double *resultsData[NR];
const char *names[NR] = { "_T1_short", "_T1_long", "_T2_short", "_T2_long", "_frac_short", "_tau_short", "_mc_B0", "_mc_res" };
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
	
	char *outPrefix = NULL, *outExt = ".nii.gz", procpar[MAXSTR];
	size_t nSPGR, nPhases, nSSFP;
	double spgrTR, *spgrAngles, ssfpTR, *ssfpPhases, *ssfpAngles,
	       *maskData = NULL, *M0Data = NULL, *B0Data = NULL, *B1Data = NULL,
	       loBounds[NR - 1], hiBounds[NR - 1], **bounds = malloc(2 * sizeof(double *));
	bounds[0] = loBounds; bounds[1] = hiBounds;
	bool loConstraint[NR - 1] = { true, true, true, true, true, true, false };
	bool hiConstraint[NR - 1] = { true, true, true, true, true, true, false };
	bool **constraints = malloc(2 * sizeof(double *));
	constraints[0] = loConstraint; constraints[1] = hiConstraint;
	FSLIO *inFile = NULL, *spgrFile = NULL, **ssfpFiles;
	short nx, ny, nz, nt;
	par_t *pars;
	static int verbose = false, start_slice = -1, end_slice = -1;
	static struct option long_options[] =
	{
		{"B0", required_argument, 0, '0'},
		{"B1", required_argument, 0, '1'},
		{"M0", required_argument, 0, 'M'},
		{"mask", required_argument, 0, 'm'},
		{"verbose", no_argument, 0, 'v'},
		{"start_slice", required_argument, 0, 'S'},
		{"end_slice", required_argument, 0, 'E'},
		{0, 0, 0, 0}
	};
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "b:m:vz", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 'm':
				fprintf(stdout, "Reading mask file %s.\n", optarg);
				inFile = FslOpen(optarg, "rb");
				maskData = FslGetVolumeAsScaledDouble(inFile, 0);
				FslClose(inFile);
				break;
			case '0':
				fprintf(stdout, "Reading B0 file %s.\n", optarg);
				inFile = FslOpen(optarg, "rb");
				B0Data = FslGetVolumeAsScaledDouble(inFile, 0);
				FslClose(inFile);
				break;
			case '1':
				fprintf(stdout, "Reading B1 file %s.\n", optarg);
				inFile = FslOpen(optarg, "rb");
				B1Data = FslGetVolumeAsScaledDouble(inFile, 0);
				FslClose(inFile);
				break;
			case 'M':
				fprintf(stdout, "Reading M0 file %s.\n", optarg);
				inFile = FslOpen(optarg, "rb");
				M0Data = FslGetVolumeAsScaledDouble(inFile, 0);
				FslClose(inFile);
				break;
			case 'v':
				verbose = true;
				break;
			case 'S':
				start_slice = atoi(optarg);
				break;
			case 'E':
				end_slice = atoi(optarg);
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
			case 0:
				// Just a flag
				break;
		}
	}
	if ((argc - optind) < 3)
	{
		fprintf(stderr, "Insufficient number of arguments.\n%s", usage);
		exit(EXIT_FAILURE);
	}
	outPrefix = argv[optind++];
	fprintf(stdout, "Output prefix will be: %s\n", outPrefix);
	//**************************************************************************
	// Gather SPGR Data
	//**************************************************************************
	fprintf(stdout, "Opening SPGR file: %s\n", argv[optind]);
	spgrFile = FslOpen(argv[optind], "rb");
	FslGetDim(spgrFile, &nx, &ny, &nz, &nt);
	nSPGR = nt;
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
	optind++;
	//**************************************************************************
	// Gather SSFP Data
	//**************************************************************************
	nPhases = argc - optind;
	ssfpFiles  = malloc(nPhases * sizeof(nifti_image *));
	ssfpPhases = malloc(nPhases * sizeof(double));
	ssfpPhases[0] = M_PI;
	ssfpFiles[0] = FslOpen(argv[optind], "rb");
	nSSFP = ssfpFiles[0]->niftiptr->nt;
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
		{
			ssfpPhases[p] = radians(realVal(pars, "rfphase", 0));
			freeProcpar(pars);
		}
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
	fprintf(stdout, "Reading SPGR data...\n");
	double *SPGR = FslGetAllVolumesAsScaledDouble(spgrFile);
	fprintf(stdout, "Reading SSFP data...\n");
	double **SSFP = malloc(nPhases * sizeof(double *));
	for (int p = 0; p < nPhases; p++)
		SSFP[p] = FslGetAllVolumesAsScaledDouble(ssfpFiles[p]);
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
		resultsData[r] = malloc(totalVoxels * sizeof(double));
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
	fprintf(stdout, "Starting processing.\n");
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
    time_t procStart = time(NULL);
	if ((start_slice < 0) || (start_slice >= nz))
		start_slice = 0;
	if ((end_slice < 0) || (end_slice > nz))
		end_slice = nz;
	for (size_t slice = start_slice; slice < end_slice; slice++)
	{
		if (verbose)
			fprintf(stdout, "Starting slice %zu...\n", slice);
		__block int voxCount = 0;
		int sliceOffset = slice * voxelsPerSlice;
		clock_t loopStart = clock();
		void (^processVoxel)(size_t vox) = ^(size_t vox)
		{
			double M0 = 1, B0 = 0, B1 = 1., residual = 0.;
			if (M0Data)
				M0 = (double)M0Data[sliceOffset + vox];
			if (B0Data)
				B0 = (double)B0Data[sliceOffset + vox];
			if (B1Data)
				B1 = (double)B1Data[sliceOffset + vox];		
			double params[NR - 1] = { 0., 0., 0., 0., 0., 0., 0. };
			if (!maskData || (maskData[sliceOffset + vox] > 0.))
			{
				AtomicAdd(1, &voxCount);
				double *signals[1 + nPhases];				
				void **consts = malloc((1 + nPhases) * sizeof(void *));
				signals[0] = malloc(nSPGR * sizeof(double));
				for (int img = 0; img < nSPGR; img++)
					signals[0][img] = SPGR[totalVoxels * img + sliceOffset +  + vox];
				consts[0] = malloc(sizeof(SPGR_constants));
				((SPGR_constants *)consts[0])->TR = spgrTR;
				((SPGR_constants *)consts[0])->M0 = M0;
				((SPGR_constants *)consts[0])->B1 = B1;				
				for (int p = 0; p < nPhases; p++)
				{
					signals[p + 1] = malloc(nSSFP * sizeof(double));
					for (int img = 0; img < nSSFP; img++)
						signals[p + 1][img] = SSFP[p][totalVoxels * img + sliceOffset +  + vox];
					consts[p + 1] = malloc(sizeof(SSFP_constants));
					((SSFP_constants *)consts[p + 1])->TR = ssfpTR;
					((SSFP_constants *)consts[p + 1])->M0 = M0;
					((SSFP_constants *)consts[p + 1])->B0 = B0;
					((SSFP_constants *)consts[p + 1])->B1 = B1;
					((SSFP_constants *)consts[p + 1])->rfPhase = ssfpPhases[p];
				}
				
				bounds[0][6] = -0.5 / ssfpTR;
				bounds[1][6] =  0.5 / ssfpTR;
				regionContraction(params, NR - 1, consts, 1 + nPhases, angles,
				                  signals, nD, false, f, bounds, constraints,
								  1000, 10, 20, 0.1, 0.025, &residual);		
				// Clean up memory
				for (int p = 0; p < 1 + nPhases; p++)
				{
					free(signals[p]);
					free(consts[p]);
				}
			}
			for (int p = 0; p < NR - 1; p++)
				resultsData[p][sliceOffset + vox]  = params[p];
			resultsData[NR - 1][sliceOffset + vox] = residual;
		};
		dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
		dispatch_apply(voxelsPerSlice, global_queue, processVoxel);
		
		if (verbose)
		{
			clock_t loopEnd = clock();
			fprintf(stdout, "Finished slice %zu", slice);
			if (voxCount > 0)
				fprintf(stdout, ", had %d unmasked voxels, CPU time per voxel was %f s.",
						voxCount, (loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC));
			fprintf(stdout, ".\n");
		}
	}
    time_t procEnd = time(NULL);
    struct tm *localEnd = localtime(&procEnd);
	char theTime[MAXSTR];
    strftime(theTime, MAXSTR, "%H:%M:%S", localEnd);
	fprintf(stdout, "Finished processing at %s. Run-time was %f s.\n", theTime, difftime(procEnd, procStart));

	for (size_t r = 0; r < NR; r++)
	{
		FslWriteHeader(resultsHeaders[r]);
		FslWriteVolumeFromDouble(resultsHeaders[r], resultsData[r], 0);
		FslClose(resultsHeaders[r]);
	}
	// Clean up memory
	free(SPGR);
	for (int p = 0; p < nPhases; p++)
		free(SSFP[p]);
	free(SSFP);
	free(M0Data);
	free(B0Data);
	free(B1Data);
	free(maskData);
	exit(EXIT_SUCCESS);
}

