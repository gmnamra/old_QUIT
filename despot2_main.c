/*
 *  despot1_main.c
 *  MacRI
 *
 *  Created by Tobias Wood on 23/01/2012.
 *  Copyright 2012 Tobias Wood. All rights reserved.
 *
 */

#include <string.h>
#include <dispatch/dispatch.h>
#include "DESPOT.h"
#include "nifti_tools.h"
#include "znzlib.h"

#define FALSE 0
#define TRUE 1

#define __DEBUG__ FALSE
#define __DEBUG_THRESH__ 60000

char *usage = "Usage is: despot2 [Output Prefix] [Number of phase cycling patterns] <[SSFP input file] [Phase Cycling]>  [SSFP TR] [Simplex/Region Contraction] [T1 Map File] [B1 Map File] <[Mask File]>\n\
\n\
For classic DESPOT2 use 1 phase cycle, simplex DESPOT2 use multiple.\n\
The input file requires 2 columns - SSFP 180 file, flip angle.\n\
Specify the phase cycling pattern in degrees after each input file.\n";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//tests();
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	if (argc == 1)
	{
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	
	if (argc < 9)
	{
		fprintf(stderr, "Incorrect number of arguments (= %d) specified, should be at least 8.\n",
				argc - 1); // Subtract off filename
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}

	char *outPrefix = argv[1];
	size_t nPhases = atoi(argv[2]);
	
	size_t *nSSFP = malloc(nPhases * sizeof(size_t));
	char **ssfpFilenames[nPhases];
	double ssfpTR, *ssfpPhases = malloc(nPhases * sizeof(double)),
	       **ssfpAngles = malloc(nPhases * sizeof(double *));

	fprintf(stdout, "Specified %ld phase cycling patterns.\n", nPhases);	
	for (size_t p = 0; p < nPhases; p++)
	{
		nSSFP[p] = readRecordFile(argv[3 + (p * 2)], "sd", &(ssfpFilenames[p]), &(ssfpAngles[p]));
		ssfpPhases[p] = radians(atof(argv[4 + (p * 2)]));
		for (int i = 0; i < nSSFP[p]; i++)
			ssfpAngles[p][i] = radians(ssfpAngles[p][i]);
		fprintf(stdout, "Specified pattern %ld with phase %f deg and %ld files.\n", p, degrees(ssfpPhases[p]), nSSFP[p]);
		//ARR_D( ssfpAngles[p], nSSFP[p] );
	}
	ssfpTR = atof(argv[3 + (nPhases * 2)]);
	int mode = atoi(argv[4 + (nPhases * 2)]);
	fprintf(stdout, "Specified TR of %f ms.\n", ssfpTR);
	nifti_image *T1Map = NULL, *B1Map = NULL, *mask = NULL;
	fprintf(stdout, "Reading T1 Map: %s\n", argv[5 + (nPhases * 2)]);
	T1Map = nifti_image_read(argv[5 + (nPhases * 2)], FALSE);
	fprintf(stdout, "Reading B1 Map: %s\n", argv[6 + (nPhases * 2)]);
	B1Map = nifti_image_read(argv[6 + (nPhases * 2)], FALSE);
	fprintf(stdout, "argc %d comp to %ld\n", argc, 8 + (nPhases * 2));
	if (argc == (8 + (nPhases * 2)))
	{
		fprintf(stdout, "Reading Mask: %s\n", argv[7 + (nPhases * 2)]);
		mask = nifti_image_read(argv[7 + (nPhases * 2)], FALSE);
	}
	//**************************************************************************	
	// Read in headers / Allocate memory for slices and results
	//**************************************************************************
	nifti_image ***ssfpFiles = malloc(nPhases * sizeof(nifti_image **));
	for (int p = 0; p < nPhases; p++)
	{
		ssfpFiles[p] = (nifti_image **)malloc(nSSFP[p] * sizeof(nifti_image *));
		loadHeaders(ssfpFilenames[p], ssfpFiles[p], nSSFP[p]);
		if ((p > 0) && ssfpFiles[p - 1][0]->nvox != ssfpFiles[p][0]->nvox)
		{
			fprintf(stderr, "Differing number of voxels in phase %d and %d headers.\n", p - 1, p);
			exit(EXIT_FAILURE);
		}
	}
	
	int voxelsPerSlice = ssfpFiles[0][0]->nx * ssfpFiles[0][0]->ny;	
	int totalVoxels = voxelsPerSlice * ssfpFiles[0][0]->nz;
	__block float *T2Data = malloc(totalVoxels * sizeof(float));
	__block float *M0Data = malloc(totalVoxels * sizeof(float));
	__block float *dOData = malloc(totalVoxels * sizeof(float));
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	fprintf(stdout, "Fitting DESPOT2.\n");
	
	for (int slice = 0; slice < ssfpFiles[0][0]->nz; slice++)
	//void (^processSlice)(size_t slice) = ^(size_t slice)
	{
		// Read in data
		fprintf(stdout, "Processing slice %ld...\n", slice);
		double T2 = 0., M0 = 0., dO = 0.;
		
		int sliceStart[7] = {0, 0, slice, 0, 0, 0, 0};
		int sliceDim[7] = {ssfpFiles[0][0]->nx, ssfpFiles[0][0]->ny, 1, 1, 1, 1, 1};
		float ***ssfpData = malloc(nPhases * sizeof(float **));
		double **ssfpSignal = malloc(nPhases * sizeof(double *));
		for (int p = 0; p < nPhases; p++)
		{
			ssfpData[p] = malloc(nSSFP[p] * sizeof(float *));
			ssfpSignal[p] = malloc(nSSFP[p] * sizeof(double));
			for (int i = 0; i < nSSFP[p]; i++)
			{
				ssfpData[p][i] = malloc(voxelsPerSlice * sizeof(float));
				nifti_read_subregion_image(ssfpFiles[p][i], sliceStart, sliceDim, (void**)&(ssfpData[p][i]));
			}
		}
		float *T1Data, *B1Data, *maskData;
		T1Data = malloc(voxelsPerSlice * sizeof(float));
		B1Data = malloc(voxelsPerSlice * sizeof(float));
		if (mask)
			maskData = malloc(voxelsPerSlice * sizeof(float));
		
		nifti_read_subregion_image(T1Map, sliceStart, sliceDim, (void**)&(T1Data));
		nifti_read_subregion_image(B1Map, sliceStart, sliceDim, (void**)&(B1Data));
		if (mask)
			nifti_read_subregion_image(mask, sliceStart, sliceDim, (void**)&(maskData));
		
		int sliceIndex = slice * voxelsPerSlice;	
		for (int vox = 0; vox < voxelsPerSlice; vox++)
		{
			if (!mask || (maskData[vox] > 0.))
			{
				double T1, B1; 
				for (int p = 0; p < nPhases; p++)
				{
					for (int img = 0; img < nSSFP[p]; img++)
						ssfpSignal[p][img] = (double)ssfpData[p][img][vox];
				}
				
				T1 = (double)T1Data[vox];
				B1 = (double)B1Data[vox];
				dO = 125. / ssfpTR; // Guess from Sean's code
				if (nPhases == 1)
					classicDESPOT2(ssfpAngles[0], ssfpSignal[0], nSSFP[0], ssfpTR, T1, B1, &M0, &T2);
				else
				{
					if (mode == 0)
					{	// Use phase cycle with highest intensity to bootstrap simplex
						double *pSigs = malloc(nPhases * sizeof(double));
						for (int p = 0; p < nPhases; p++)
							pSigs[p] = ssfpSignal[p][0];
						size_t *ind = malloc(nPhases * sizeof(double));
						arraySort(pSigs, nPhases, ind);
						size_t hi = ind[nPhases - 1];
						classicDESPOT2(ssfpAngles[hi], ssfpSignal[hi], nSSFP[hi], ssfpTR, T1, B1, &M0, &T2);
						simplexDESPOT2(nPhases, nSSFP, ssfpPhases, ssfpAngles, ssfpSignal, ssfpTR, T1, B1, &M0, &T2, &dO);
						
						free(ind);
						free(pSigs);
					}
					else
					{	// Use region contraction
						contractDESPOT2(nPhases, nSSFP, ssfpPhases, ssfpAngles, ssfpSignal, ssfpTR, T1, B1, &M0, &T2, &dO);
					}
				}
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
		for (int p = 0; p < nPhases; p++)
		{
			for (int img = 0; img < nSSFP[p]; img++)
				free(ssfpData[p][img]);
			free(ssfpData[p]);
			free(ssfpSignal[p]);
		}
		free(ssfpSignal);
		free(T1Data);
		free(B1Data);
	};
	//dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
	//dispatch_apply(ssfpFiles[0][0]->nz, global_queue, processSlice);
	fprintf(stdout, "Finished fitting. Writing results files.\n");

	//**************************************************************************
	// Create results files
	//**************************************************************************
	nifti_image *out = nifti_copy_nim_info(ssfpFiles[0][0]);
	char outName[strlen(outPrefix) + 4]; // Space for "T1" plus null
	strcpy(outName, outPrefix); strcat(outName, "_M0");
	writeResult(out, outName, (void*)M0Data);
	strcpy(outName, outPrefix); strcat(outName, "_T2");
	writeResult(out, outName, (void*)T2Data);
	if (nPhases > 1)
	{
		strcpy(outName, outPrefix); strcat(outName, "_dO");
		writeResult(out, outName, (void*)dOData);
	}
	fprintf(stdout, "All done.\n");
	return EXIT_SUCCESS;
}
