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
#include <libkern/OSAtomic.h>
#include "DESPOT.h"
#include "nifti_tools.h"
#include "znzlib.h"
#include "mathArray.h"

#define __DEBUG__ FALSE
#define __DEBUG_THRESH__ 60000

char *usage = "Usage is: mcdespot [options] output_prefix spgr_file ssfp_file1 (ssfp_fileN)\n\
\
Options:\n\
	-m file  : Mask input with specified file.\n\
	-z       : Output .nii.gz files.\n\
	-b file  : B1 Map file.\n";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	char *outPrefix = NULL, *outExt = ".nii";
	size_t nSPGR, nPhases, *nSSFP;
	double spgrTR, *spgrAngles,
		   ssfpTR, *ssfpPhases, **ssfpAngles;
	nifti_image *spgrFile = NULL, *maskFile = NULL, *B1File, **ssfpFiles = NULL;
	
	if (argc < 3)
	{
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	
	int thisArg = 1;
	while ((thisArg < argc) && (argv[thisArg][0] =='-'))
	{
		switch (argv[thisArg][1])
		{
			case 'm':
				fprintf(stdout, "Reading Mask from: %s\n", argv[++thisArg]);
				maskFile = nifti_image_read(argv[thisArg], FALSE);
				break;
			case 'b':
				fprintf(stdout, "Reading B1 Map from: %s\n", argv[++thisArg]);
				B1File = nifti_image_read(argv[thisArg], FALSE);
				break;
			case 'z':
				outExt = ".nii.gz";
				break; 
			default:
				fprintf(stderr, "Undefined command line option\n%s", usage);
				exit(EXIT_FAILURE);
				break;
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
	spgrAngles = malloc(nSPGR * sizeof(double));
	fprintf(stdout, "Enter SPGR TR (ms):");
	fscanf(stdin, "%lf", &spgrTR);
	fprintf(stdout, "Enter %zu SPGR Flip Angles (degrees):", nSPGR);
	fgetArray(stdin, 'd', nSPGR, spgrAngles);
	
	fprintf(stdout, "SPGR TR=%f ms. ", spgrTR);
	ARR_D(spgrAngles, nSPGR);
	arrayApply(spgrAngles, spgrAngles, radians, nSPGR);

	//**************************************************************************
	// Gather SSFP Data
	//**************************************************************************
	nPhases    = argc - thisArg;
	nSSFP      = malloc(nPhases * sizeof(int));
	ssfpPhases = malloc(nPhases * sizeof(double));
	ssfpAngles = malloc(nPhases * sizeof(double *));
	ssfpFiles  = malloc(nPhases * sizeof(nifti_image *));
	fprintf(stdout, "Specified %zu phase cycling patterns.\n", nPhases);
	fprintf(stdout, "Enter SSFP TR (ms):");
	fscanf(stdin, "%lf", &ssfpTR);
	fprintf(stdout, "Enter %zu Phase-Cycling Patterns (degrees):", nPhases);
	fgetArray(stdin, 'd', nPhases, ssfpPhases);
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
		ssfpAngles[p] = malloc(nSSFP[p] * sizeof(double));
		fprintf(stdout, "Enter %zu Angles (degrees) for Pattern %d:", nSSFP[p], p);
		fgetArray(stdin, 'd', nSSFP[p], ssfpAngles[p]);
		arrayApply(ssfpAngles[p], ssfpAngles[p], radians, nSSFP[p]);
		fprintf(stdout, "Specified pattern %d with phase %f deg.\n", p, degrees(ssfpPhases[p]));
		ARR_D(ssfpAngles[p], nSSFP[p]);
	}
	if ((spgrFile->nx * spgrFile->ny * spgrFile->nz) !=
	    (ssfpFiles[0]->nx * ssfpFiles[0]->ny * ssfpFiles[0]->nz))
	{
		fprintf(stderr, "Differing number of voxels in SPGR and SSFP data.\n");
		exit(EXIT_FAILURE);
	}
	
	//**************************************************************************	
	// Allocate memory for slices and results
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
	const int NR = 8;
	nifti_image **resultsHeaders = malloc(NR * sizeof(nifti_image *));
	int outDims[8] = {3, spgrFile->nx, spgrFile->ny, spgrFile->nz, 1, 1, 1, 1};
	float *blank = calloc(totalVoxels, sizeof(float));
	float **resultsSlices = malloc(NR * sizeof(float*));
	char *names[NR] = { "_T1_myel", "_T1_free", "_T2_myel", "_T2_free", "_frac_myel", "_tau_myel", "_dw", "_res" };
	char outName[strlen(outPrefix) + 18];
	for (int p = 0; p < NR; p++)
	{
		strcpy(outName, outPrefix); strcat(outName, names[p]); strcat(outName, outExt);
		fprintf(stdout, "Writing blank result file:%s.\n", outName);
		resultsHeaders[p] = nifti_copy_orientation(spgrFile, outDims, DT_FLOAT);
		nifti_set_filenames(resultsHeaders[p], outName, FALSE, TRUE);
		resultsHeaders[p]->data = blank;
		nifti_image_write(resultsHeaders[p]);
		resultsSlices[p] = malloc(voxelsPerSlice * sizeof(float));
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
		int sliceStart[NR] = {0, 0, (int)slice, 0, 0, 0, 0};
		int sliceDim[NR] = {spgrFile->nx, spgrFile->ny, 1, 1, 1, 1, 1};
		
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
			double params[NR];
			arraySet(params, 0., NR);
			if (!maskFile || (maskData[vox] > 0.))
			{
				OSAtomicAdd32(1, &voxCount);
				
				double spgrSignal[nSPGR];
				for (int img = 0; img < nSPGR; img++)
					spgrSignal[img] = (double)spgrData[voxelsPerSlice * img + vox];
				
				double *ssfpSignal[nPhases];
				for (int p = 0; p < nPhases; p++)
				{
					ssfpSignal[p] = malloc(nSSFP[p] * sizeof(double));
					for (int img = 0; img < nSSFP[p]; img++)
						ssfpSignal[p][img] = (double)ssfpData[p][voxelsPerSlice * img + vox];
				}
				
				double B1 = 1.;
				if (B1File)
					B1 = (double)B1Data[vox];

				mcDESPOT(nSPGR, spgrAngles, spgrSignal, spgrTR,
						 nPhases, nSSFP, ssfpPhases, ssfpAngles,
						 ssfpSignal, ssfpTR, B1, params);
				
				for (int p = 0; p < nPhases; p++)
					free(ssfpSignal[p]);
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
	return EXIT_SUCCESS;
	
	// Clean up memory
	free(spgrData);
	for (int p = 0; p < nPhases; p++)
		free(ssfpData[p]);
	free(B1Data);
	free(maskData);
}

