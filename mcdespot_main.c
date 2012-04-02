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

#ifdef __APPLE__
	#include <libkern/OSAtomic.h>
	#define AtomicAdd OSAtomicAdd32
#else
	#define AtomicAdd(x, y) (*y)++
#endif

#include "DESPOT.h"
#include "nifti_tools.h"
#include "znzlib.h"
#include "mathArray.h"

char *usage = "Usage is: mcdespot [options] output_prefix spgr_file ssfp_file1 (ssfp_fileN)\n\
\
Options:\n\
	-m file  : Mask input with specified file.\n\
	-z       : Output .nii.gz files.\n\
	-B1 file : B1 Map file.\n\
	-b 3     : Boundaries suitable for 3T (default)\n\
	   7     : Boundaries suitable for 7T\n\
	   u     : User specified boundaries from stdin.\n";

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
	size_t nSPGR, nPhases, *nSSFP;
	float spgrTR, *spgrAngles,
		   ssfpTR, *ssfpPhases, **ssfpAngles, loBounds[7], hiBounds[7],
		   lo3Bounds[7] = { 200.,  500.,  2.,  40., 0.0,   50., 0. },
		   hi3Bounds[7] = { 700., 2500., 40., 200., 0.4, 2000., 0. },
		   lo7Bounds[7] = {  200.,  500.,   2.,  40., 0.,   50., 0. },
		   hi7Bounds[7] = { 1500., 5000., 100., 400., 1., 2000., 0. },
		   **bounds = malloc(2 * sizeof(float *));
	bounds[0] = loBounds; bounds[1] = hiBounds;
	bool loConstraint[7] = { true, true, true, true, true, true, false };
	bool hiConstraint[7] = { true, true, true, true, true, true, false };
	bool **constraints = malloc(2 * sizeof(bool *));
	constraints[0] = loConstraint; constraints[1] = hiConstraint;
	nifti_image *spgrFile = NULL, *maskFile = NULL, *B1File, **ssfpFiles = NULL;
		
	int thisArg = 1;
	while ((thisArg < argc) && (argv[thisArg][0] =='-'))
	{
		if (strcmp(argv[thisArg], "-m") == 0) {
			fprintf(stdout, "Reading Mask from: %s\n", argv[++thisArg]);
			maskFile = nifti_image_read(argv[thisArg], FALSE);
		} else if (strcmp(argv[thisArg], "-B1") == 0) {
			fprintf(stdout, "Reading B1 Map from: %s\n", argv[++thisArg]);
			B1File = nifti_image_read(argv[thisArg], FALSE);
		} else if (strcmp(argv[thisArg], "-z") == 0) {
			outExt = ".nii.gz";
		} else if (strcmp(argv[thisArg], "-b") == 0) {
			switch (*argv[++thisArg])
			{
				case '3':
					fprintf(stdout, "Using 3T boundaries.\n");
					arrayCopy(loBounds, lo3Bounds, 7);
					arrayCopy(hiBounds, hi3Bounds, 7);
					break;
				case '7':
					fprintf(stdout, "Using 7T boundaries.\n");
					arrayCopy(loBounds, lo7Bounds, 7);
					arrayCopy(hiBounds, hi7Bounds, 7);
					break;
				case 'u':
					fprintf(stdout, "Enter low boundaries (6 values):");
					fgetArray(stdin, 'f', 6, loBounds);
					fprintf(stdout, "Enter high boundaries (6 values):");
					fgetArray(stdin, 'f', 6, hiBounds);
					break;
				default:
					fprintf(stdout, "Unknown boundaries type '%c'.\n", argv[thisArg][0]);
					exit(EXIT_FAILURE);
					break;
			}
		} else {
			fprintf(stderr, "Undefined command line option\n%s", usage);
			exit(EXIT_FAILURE);
		}
		++thisArg;
	}
	outPrefix = argv[thisArg]; thisArg++;
	fprintf(stdout, "Output prefix will be: %s\n", outPrefix);
	//**************************************************************************
	// Gather SPGR Data
	//**************************************************************************
	fprintf(stdout, "Reading SPGR header from %s.\n", argv[thisArg]);
	spgrFile = nifti_image_read(argv[thisArg], FALSE); thisArg++;
	nSPGR = spgrFile->nt;
	spgrAngles = malloc(nSPGR * sizeof(float));
	fprintf(stdout, "Enter SPGR TR (ms):");
	fscanf(stdin, "%f", &spgrTR);
	fprintf(stdout, "Enter %zu SPGR Flip Angles (degrees):", nSPGR);
	fgetArray(stdin, 'f', nSPGR, spgrAngles);
	arrayApply(spgrAngles, spgrAngles, radians, nSPGR);

	//**************************************************************************
	// Gather SSFP Data
	//**************************************************************************
	nPhases    = argc - thisArg;
	nSSFP      = malloc(nPhases * sizeof(int));
	ssfpPhases = malloc(nPhases * sizeof(float));
	ssfpAngles = malloc(nPhases * sizeof(float *));
	ssfpFiles  = malloc(nPhases * sizeof(nifti_image *));
	fprintf(stdout, "Specified %zu phase cycling patterns.\n", nPhases);
	fprintf(stdout, "Enter SSFP TR (ms):");
	fscanf(stdin, "%f", &ssfpTR);
	fprintf(stdout, "Enter %zu Phase-Cycling Patterns (degrees):", nPhases);
	fgetArray(stdin, 'f', nPhases, ssfpPhases); fprintf(stdout, "\n");
	arrayApply(ssfpPhases, ssfpPhases, radians, nPhases);
	for (int p = 0; p < nPhases; p++)
	{
		fprintf(stdout, "Reading SSFP header from %s.\n", argv[thisArg]);
		ssfpFiles[p] = nifti_image_read(argv[thisArg], FALSE); thisArg++;
		if ((p > 0) && (ssfpFiles[p - 1]->nx * ssfpFiles[p - 1]->ny * ssfpFiles[p - 1]->nz) != 
		               (ssfpFiles[0]->nx * ssfpFiles[0]->ny * ssfpFiles[0]->nz))
		{
			fprintf(stderr, "Differing number of voxels in phase %d and %d headers.\n", p - 1, p);
			exit(EXIT_FAILURE);
		}
		nSSFP[p] = ssfpFiles[p]->nt;
		ssfpAngles[p] = malloc(nSSFP[p] * sizeof(float));
		fprintf(stdout, "Enter %zu Angles (degrees) for Pattern %d:", nSSFP[p], p);
		fgetArray(stdin, 'f', nSSFP[p], ssfpAngles[p]); fprintf(stdout, "\n");
		arrayApply(ssfpAngles[p], ssfpAngles[p], radians, nSSFP[p]);
	}
	if ((spgrFile->nx * spgrFile->ny * spgrFile->nz) !=
	    (ssfpFiles[0]->nx * ssfpFiles[0]->ny * ssfpFiles[0]->nz))
	{
		fprintf(stderr, "Differing number of voxels in SPGR and SSFP data.\n");
		exit(EXIT_FAILURE);
	}
	
	//**************************************************************************	
	// Allocate memory for slices
	//**************************************************************************
	int voxelsPerSlice = spgrFile->nx * spgrFile->ny;	
	int totalVoxels = voxelsPerSlice * spgrFile->nz;
	float *spgrData, *maskData, *B1Data, **ssfpData;
	spgrData = malloc(nSPGR * voxelsPerSlice * sizeof(float));
	ssfpData = malloc(nPhases * sizeof(float *));
	for (int p = 0; p < nPhases; p++)
		ssfpData[p] = malloc(nSSFP[p] * voxelsPerSlice * sizeof(float));
	if (maskFile)
		maskData = malloc(voxelsPerSlice * sizeof(float));
	if (B1File)
		B1Data = malloc(voxelsPerSlice * sizeof(float));
	
	//**************************************************************************
	// Create results files
	// T1_m, T1_m, T2_f, T2_f,	f_m, tau_m, dw, residue
	// Need to write a full file of zeros first otherwise per-plane writing
	// won't produce a complete image.
	//**************************************************************************
	#define NR 8
	nifti_image **resultsHeaders = malloc(NR * sizeof(nifti_image *));
	int outDims[8] = {3, spgrFile->nx, spgrFile->ny, spgrFile->nz, 1, 1, 1, 1};
	float *blank = calloc(totalVoxels, sizeof(float));
	float **resultsSlices = malloc(NR * sizeof(float*));
	char *names[NR] = { "_T1_myel", "_T1_free", "_T2_myel", "_T2_free", "_frac_myel", "_tau_myel", "_dw", "_res" };
	char outName[strlen(outPrefix) + 18];
	for (int r = 0; r < NR; r++)
	{
		strcpy(outName, outPrefix); strcat(outName, names[r]); strcat(outName, outExt);
		fprintf(stdout, "Writing blank result file:%s.\n", outName);
		resultsHeaders[r] = nifti_copy_orientation(spgrFile, outDims, DT_FLOAT);
		nifti_set_filenames(resultsHeaders[r], outName, FALSE, TRUE);
		resultsHeaders[r]->data = blank;
		nifti_image_write(resultsHeaders[r]);
		resultsSlices[r] = malloc(voxelsPerSlice * sizeof(float));
	}

	//**************************************************************************
	// Set up functions and angles together for regionContraction
	//**************************************************************************
	size_t *nD = malloc((1 + nPhases) * sizeof(size_t));
	eval_array_type **f = malloc((1 + nPhases) * sizeof(eval_array_type *));
	float **angles = malloc((1 + nPhases) * sizeof(float *));

	loBounds[6] = 0.; hiBounds[6] = 1. / ssfpTR;
	nD[0] = nSPGR;
	f[0]      = a2cSPGR;
	angles[0] = spgrAngles;
	for (int p = 0; p < nPhases; p++)
	{
		nD[p + 1] = nSSFP[p];
		f[p + 1]      = a2cSSFP;
		angles[p + 1] = ssfpAngles[p];
	}
	
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	time_t allStart = time(NULL);
	struct tm *localStart = localtime(&allStart);
	char theTime[1024];
	strftime(theTime, 1024, "%H:%M:%S", localStart);
	fprintf(stdout, "Started processing at %s.\n", theTime);
	for (size_t slice = 0; slice < spgrFile->nz; slice++)
	{
		// Read in data
		fprintf(stdout, "Starting slice %zu...\n", slice);
		int sliceStart[7] = {0, 0, (int)slice, 0, 0, 0, 0};
		int sliceDim[7] = {spgrFile->nx, spgrFile->ny, 1, 1, 1, 1, 1};
		
		if (B1File)
			nifti_read_subregion_image(B1File, sliceStart, sliceDim, (void**)&(B1Data));
		if (maskFile)
		{
			if (maskFile->datatype == DT_FLOAT)
				nifti_read_subregion_image(maskFile, sliceStart, sliceDim, (void**)&(maskData));
			else
			{
				nifti_read_subregion_image(maskFile, sliceStart, sliceDim, &(maskFile->data));
				arrayConvert(maskData, maskFile->data, DT_FLOAT, maskFile->datatype, voxelsPerSlice);
			}
		}
		sliceDim[3] = nSPGR;
		nifti_read_subregion_image(spgrFile, sliceStart, sliceDim, (void**)&(spgrData));
		
		for (int p = 0; p < nPhases; p++)
		{
			sliceDim[3] = nSSFP[p];
			nifti_read_subregion_image(ssfpFiles[p], sliceStart, sliceDim, (void**)&(ssfpData[p]));
		}
		
		__block int voxCount = 0;
		clock_t loopStart = clock();
		void (^processVoxel)(size_t vox) = ^(size_t vox)
		{
			float params[NR];
			arraySet(params, 0., NR);
			if (!maskFile || (maskData[vox] > 0.))
			{
				AtomicAdd(1, &voxCount);
				
				float B1 = 1.;
				if (B1File)
					B1 = (float)B1Data[vox];
				
				// Constants need to be set up here because B1 changes per-voxels
				float *signals[1 + nPhases], *consts[1 + nPhases];
				
				signals[0] = malloc(nSPGR * sizeof(float));
				for (int img = 0; img < nSPGR; img++)
					signals[0][img] = (float)spgrData[voxelsPerSlice * img + vox];
				arrayScale(signals[0], signals[0], 1. / arrayMean(signals[0], nSPGR), nSPGR);
				consts[0] = malloc(2 * sizeof(float));
				consts[0][0] = spgrTR; consts[0][1] = B1;
								
				for (int p = 0; p < nPhases; p++)
				{
					signals[p + 1] = malloc(nSSFP[p] * sizeof(float));
					for (int img = 0; img < nSSFP[p]; img++)
						signals[p + 1][img] = (float)ssfpData[p][voxelsPerSlice * img + vox];
					arrayScale(signals[p + 1], signals[p + 1], 1. / arrayMean(signals[p + 1], nSSFP[p]), nSSFP[p]);
					
					consts[p + 1] = malloc(3 * sizeof(float));
					consts[p + 1][0] = ssfpTR; consts[p + 1][1] = B1; consts[p + 1][2] = ssfpPhases[p];
				}

				// Store residual in final parameter
				regionContraction(params, 7, consts, 1 + nPhases, angles,
				                  signals, nD, true, f, bounds, constraints,
								  5000, 25, 20, 0.1, 0.025, &(params[7]));
				params[6] = fmod(params[6], 1./ssfpTR); // Bring B0 back to one cycle 
				if (params[6] < -.5/ssfpTR)
					params[6] += 1./ssfpTR;  // And correct for signed modulus in C
				if (params[6] > .5/ssfpTR)
					params[6] -= 1./ssfpTR;
				params[6] *= 1.e3;           // Finally convert to Hz
				
				// Clean up memory
				for (int p = 0; p < 1 + nPhases; p++)
				{
					free(signals[p]);
					free(consts[p]);
				}
			}
			for (int p = 0; p < NR; p++)
				resultsSlices[p][vox]  = (float)params[p];
		};
		dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
		dispatch_apply(voxelsPerSlice, global_queue, processVoxel);
		
        clock_t loopEnd = clock();
        fprintf(stdout, "Finished slice %zu", slice);
		if (voxCount > 0)
		{
			fprintf(stdout, ", had %d unmasked voxels, CPU time per voxel was %f s. Writing to results files...", 
			        voxCount, (loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC));
			for (int p = 0; p < NR; p++)
				nifti_write_subregion_image(resultsHeaders[p], sliceStart, sliceDim, (void **)&(resultsSlices[p]));
            fprintf(stdout, "done");
		}
		fprintf(stdout, ".\n");
	}
    time_t allEnd = time(NULL);
    struct tm *localEnd = localtime(&allEnd);
    strftime(theTime, 1024, "%H:%M:%S", localEnd);
	fprintf(stdout, "Finished processing at %s. Run-time was %f s.\n", theTime, difftime(allEnd, allStart));

	// Clean up memory
	free(spgrData);
	for (int p = 0; p < nPhases; p++)
		free(ssfpData[p]);
	free(B1Data);
	free(maskData);
	exit(EXIT_SUCCESS);
}

