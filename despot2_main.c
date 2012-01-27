/*
 *  despot1_main.c
 *  MacRI
 *
 *  Created by Tobias Wood on 23/01/2012.
 *  Copyright 2012 Tobias Wood. All rights reserved.
 *
 */

#include <string.h>
#ifdef __APPLE__
	#include <dispatch/dispatch.h>
#else
	#define __block 
#endif
#include "DESPOT.h"
#include "nifti_tools.h"
#include "znzlib.h"

#define FALSE 0
#define TRUE 1

#define __DEBUG__ FALSE
#define __DEBUG_THRESH__ 60000

char *usage = "Usage is: despot2 [Output Prefix] [SSFP input file] [SSFP TR] [T1 Map File] [B1 Map File] <[Mask File]>\n\
\n\
Input File format:\n\
1st line - single integer (N) specifying number of files.\n\
N lines  - Path to file, space, parameter. This is the flip angle in degrees\n";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//tests();
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	int nSSFP = 0;
	char **ssfpFilenames;
	double ssfpTR, *ssfpAngles;
	char *outPrefix;
	
	if (argc == 1)
	{
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	
	if ((argc != 6) && (argc != 7))
	{
		fprintf(stderr, "Incorrect number of arguments (= %d) specified, should be 5 for DESPOT2 (6 if masked).\n",
				argc - 1); // Subtract off filename
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	
	outPrefix = argv[1];
	nSSFP = readRecordFile(argv[2], "sd", &ssfpFilenames, &ssfpAngles);
	for (int i = 0; i < nSSFP; i++)
		ssfpAngles[i] = radians(ssfpAngles[i]);
	ssfpTR = atof(argv[3]);
	fprintf(stdout, "Specified %d SSFP files with TR=%f ms.\n", nSSFP, ssfpTR);
	
	nifti_image *T1Map = NULL, *B1Map = NULL, *mask = NULL;
	T1Map = nifti_image_read(argv[4], FALSE);
	B1Map = nifti_image_read(argv[5], FALSE);
	if (argc == 7)
		mask = nifti_image_read(argv[6], FALSE);
	//**************************************************************************	
	// Read in headers / Allocate memory for slices and results
	//**************************************************************************
	nifti_image **ssfpFiles = (nifti_image **)malloc(nSSFP * sizeof(nifti_image *));
	loadHeaders(ssfpFilenames, ssfpFiles, nSSFP);
	
	int voxelsPerSlice = ssfpFiles[0]->nx * ssfpFiles[0]->ny;	
	int totalVoxels = voxelsPerSlice * ssfpFiles[0]->nz;
	__block float *T2Data = (float *)malloc(totalVoxels * sizeof(float));
	__block float *M0Data = (float *)malloc(totalVoxels * sizeof(float));
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	fprintf(stdout, "Fitting DESPOT2.\n");
	
	for (int slice = 0; slice < ssfpFiles[0]->nz; slice++)
	//void (^processSlice)(size_t slice) = ^(size_t slice)
	{
		// Read in data
		fprintf(stdout, "Processing slice %ld...\n", slice);
		double T2 = 0., M0 = 0.; // Place to restore per-voxel return values, assume B1 field is uniform for classic DESPOT
		
		int sliceStart[7] = {0, 0, slice, 0, 0, 0, 0};
		int sliceDim[7] = {ssfpFiles[0]->nx, ssfpFiles[0]->ny, 1, 1, 1, 1, 1};
		float *ssfpData[nSSFP];
		for (int i = 0; i < nSSFP; i++)
			ssfpData[i] = malloc(voxelsPerSlice * sizeof(float));
		float *T1Data, *B1Data, *maskData;
		T1Data = malloc(voxelsPerSlice * sizeof(float));
		B1Data = malloc(voxelsPerSlice * sizeof(float));
		if (mask)
			maskData = malloc(voxelsPerSlice * sizeof(float));
		
		for (int img = 0; img < nSSFP; img++)
			nifti_read_subregion_image(ssfpFiles[img], sliceStart, sliceDim, (void**)&(ssfpData[img]));
		nifti_read_subregion_image(T1Map, sliceStart, sliceDim, (void**)&(T1Data));
		nifti_read_subregion_image(B1Map, sliceStart, sliceDim, (void**)&(B1Data));
		if (mask)
			nifti_read_subregion_image(mask, sliceStart, sliceDim, (void**)&(maskData));
		
		int sliceIndex = slice * voxelsPerSlice;	
		for (int vox = 0; vox < voxelsPerSlice; vox++)
		{
			if (!mask || (maskData[vox] > 0.))
			{
				double ssfpSignal[nSSFP], T1, B1;
				for (int img = 0; img < nSSFP; img++)
					ssfpSignal[img] = (double)ssfpData[img][vox];		
				T1 = (double)T1Data[vox];
				B1 = (double)B1Data[vox];
				calcDESPOT2(ssfpAngles, ssfpSignal, nSSFP, ssfpTR, T1, B1, &M0, &T2);
				// Sanity check
				M0 = clamp(M0, 0., 1.e8);
				T2 = clamp(T2, 0., 3.e3);
			}
			else
			{
				M0 = 0.; T2 = 0.;
			}

						
			T2Data[sliceIndex + vox] = (float)T2;
			M0Data[sliceIndex + vox] = (float)M0;
		}
		
		// Clean up memory
		for (int i = 0; i < nSSFP; i++)
			free(ssfpData[i]);
	};
	//dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
	//dispatch_apply(SPGRFiles[0]->nz, global_queue, processSlice);
	fprintf(stdout, "Finished fitting. Writing results files.\n");

	//**************************************************************************
	// Create results files
	//**************************************************************************
	nifti_image *out = nifti_copy_nim_info(ssfpFiles[0]);
	char outName[strlen(outPrefix) + 4]; // Space for "T1" plus null
	strcpy(outName, outPrefix); strcat(outName, "_M0");
	writeResult(out, outName, (void*)M0Data);
	strcpy(outName, outPrefix); strcat(outName, "_T2");
	writeResult(out, outName, (void*)T2Data);
	fprintf(stdout, "All done.\n");
	return EXIT_SUCCESS;
}
