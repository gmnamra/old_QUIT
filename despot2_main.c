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

enum DESPOT2_MODE
{
	MODE_CLASSIC = 0,
	MODE_SIMPLEX
};

char *usage = "Usage is: despot2 [Output Prefix] [Classic/Simplex Switch] [SSFP0 input file] [SSFP TR] [T1 Map File] [B1 Map File] <[Mask File]>\n\
\n\
For classic DESPOT2 set switch to 0, simplex DESPOT2 set to 1.\n\
For classic DESPOT2 the input file requires 2 columns - SSFP 180 file, flip angle.\n\
For simplex DESPOT2 the input file requires 3 columns - SSFP 0 Phase file, SSFP 180 Phase file, Flip Angle.";

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
	char **ssfp0Filenames, **ssfp180Filenames;
	double ssfpTR, *ssfpAngles;
	char *outPrefix;
	
	if (argc == 1)
	{
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	
	if ((argc != 7) && (argc != 8))
	{
		fprintf(stderr, "Incorrect number of arguments (= %d) specified, should be 6 (7 if masked).\n",
				argc - 1); // Subtract off filename
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	
	outPrefix = argv[1];
	int mode = atoi(argv[2]);
	if (mode == MODE_CLASSIC)
	{
		fprintf(stdout, "Specified classic DESPOT2 fitting.\n");
		nSSFP = readRecordFile(argv[3], "sd", &ssfp180Filenames, &ssfpAngles);
	}
	else if (mode == MODE_SIMPLEX)
	{
		fprintf(stdout, "Specified simplex full DESPOT2 fitting.\n");
		nSSFP = readRecordFile(argv[3], "ssd", &ssfp0Filenames, &ssfp180Filenames, &ssfpAngles);
	}
	else
	{
		fprintf(stderr, "Unknown mode %d specified.\n", mode);
		exit(EXIT_FAILURE);
	}


	for (int i = 0; i < nSSFP; i++)
		ssfpAngles[i] = radians(ssfpAngles[i]);
	ssfpTR = atof(argv[4]);
	fprintf(stdout, "Specified %d SSFP files with TR=%f ms.\n", nSSFP, ssfpTR);
	
	nifti_image *T1Map = NULL, *B1Map = NULL, *mask = NULL;
	T1Map = nifti_image_read(argv[5], FALSE);
	B1Map = nifti_image_read(argv[6], FALSE);
	if (argc == 8)
		mask = nifti_image_read(argv[7], FALSE);
	//**************************************************************************	
	// Read in headers / Allocate memory for slices and results
	//**************************************************************************
	nifti_image **ssfp0Files = (nifti_image **)malloc(nSSFP * sizeof(nifti_image *));
	nifti_image **ssfp180Files = (nifti_image **)malloc(nSSFP * sizeof(nifti_image *));
	loadHeaders(ssfp180Filenames, ssfp180Files, nSSFP);
	if (mode == MODE_SIMPLEX)
	{
		loadHeaders(ssfp0Filenames, ssfp0Files, nSSFP);
		if (ssfp0Files[0]->nvox != ssfp180Files[0]->nvox)
		{
			fprintf(stderr, "Differing number of voxels in phase 0 and phase 180 headers.\n");
			exit(EXIT_FAILURE);
		}
	}
	
	int voxelsPerSlice = ssfp180Files[0]->nx * ssfp180Files[0]->ny;	
	int totalVoxels = voxelsPerSlice * ssfp180Files[0]->nz;
	__block float *T2Data = malloc(totalVoxels * sizeof(float));
	__block float *M0Data = malloc(totalVoxels * sizeof(float));
	__block float *dOData = malloc(totalVoxels * sizeof(float));
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	fprintf(stdout, "Fitting DESPOT2.\n");
	
	for (int slice = 0; slice < ssfp180Files[0]->nz; slice++)
	//void (^processSlice)(size_t slice) = ^(size_t slice)
	{
		// Read in data
		fprintf(stdout, "Processing slice %ld...\n", slice);
		double T2 = 0., M0 = 0., dO = 0.;
		
		int sliceStart[7] = {0, 0, slice, 0, 0, 0, 0};
		int sliceDim[7] = {ssfp180Files[0]->nx, ssfp180Files[0]->ny, 1, 1, 1, 1, 1};
		float *ssfp180Data[nSSFP], *ssfp0Data[nSSFP];
		for (int i = 0; i < nSSFP; i++)
		{
			ssfp180Data[i] = malloc(voxelsPerSlice * sizeof(float));
			if (mode == MODE_SIMPLEX)
				ssfp0Data[i] = malloc(voxelsPerSlice * sizeof(float));
		}
		float *T1Data, *B1Data, *maskData;
		T1Data = malloc(voxelsPerSlice * sizeof(float));
		B1Data = malloc(voxelsPerSlice * sizeof(float));
		if (mask)
			maskData = malloc(voxelsPerSlice * sizeof(float));
		
		for (int img = 0; img < nSSFP; img++)
		{
			nifti_read_subregion_image(ssfp180Files[img], sliceStart, sliceDim, (void**)&(ssfp180Data[img]));
			if (mode == MODE_SIMPLEX)
				nifti_read_subregion_image(ssfp0Files[img], sliceStart, sliceDim, (void**)&(ssfp0Data[img]));
		}
		nifti_read_subregion_image(T1Map, sliceStart, sliceDim, (void**)&(T1Data));
		nifti_read_subregion_image(B1Map, sliceStart, sliceDim, (void**)&(B1Data));
		if (mask)
			nifti_read_subregion_image(mask, sliceStart, sliceDim, (void**)&(maskData));
		
		int sliceIndex = slice * voxelsPerSlice;	
		for (int vox = 0; vox < voxelsPerSlice; vox++)
		{
			if (!mask || (maskData[vox] > 0.))
			{
				double ssfp180Signal[nSSFP], ssfp0Signal[nSSFP], T1, B1;
				for (int img = 0; img < nSSFP; img++)
				{
					ssfp180Signal[img] = (double)ssfp180Data[img][vox];
					if (mode == MODE_SIMPLEX)
						ssfp0Signal[img] = (double)ssfp0Data[img][vox];
				}
				
				T1 = (double)T1Data[vox];
				B1 = (double)B1Data[vox];
				dO = 125. / ssfpTR; // Guess from Sean's code
				classicDESPOT2(ssfpAngles, ssfp180Signal, nSSFP, ssfpTR, T1, B1, &M0, &T2);
				if (mode == MODE_SIMPLEX)
					simplexDESPOT2(ssfpAngles, ssfp0Signal, ssfp180Signal, nSSFP, ssfpTR, T1, B1, &M0, &T2, &dO);
				// Sanity check
				M0 = clamp(M0, 0., 1.e8);
				T2 = clamp(T2, 0., 3.e3);
				dO = clamp(dO, 0., INFINITY);
			}
			else
			{
				M0 = 0.; T2 = 0.; dO = 0.;
			}

						
			T2Data[sliceIndex + vox] = (float)T2;
			M0Data[sliceIndex + vox] = (float)M0;
			dOData[sliceIndex + vox] = (float)dO;
		}
		
		// Clean up memory
		for (int i = 0; i < nSSFP; i++)
		{
			free(ssfp180Data[i]);
			if (mode == MODE_SIMPLEX)
				free(ssfp0Data[i]);
		}
	};
	//dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
	//dispatch_apply(ssfp180Files[0]->nz, global_queue, processSlice);
	fprintf(stdout, "Finished fitting. Writing results files.\n");

	//**************************************************************************
	// Create results files
	//**************************************************************************
	nifti_image *out = nifti_copy_nim_info(ssfp180Files[0]);
	char outName[strlen(outPrefix) + 4]; // Space for "T1" plus null
	strcpy(outName, outPrefix); strcat(outName, "_M0");
	writeResult(out, outName, (void*)M0Data);
	strcpy(outName, outPrefix); strcat(outName, "_T2");
	writeResult(out, outName, (void*)T2Data);
	if (mode == MODE_SIMPLEX)
	{
		strcpy(outName, outPrefix); strcat(outName, "_dO");
		writeResult(out, outName, (void*)dOData);
	}
	fprintf(stdout, "All done.\n");
	return EXIT_SUCCESS;
}
