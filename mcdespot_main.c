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
#include <dispatch/dispatch.h>
#include "DESPOT.h"
#include "nifti_tools.h"
#include "znzlib.h"

#define __DEBUG__ FALSE
#define __DEBUG_THRESH__ 60000

char *usage = "Usage is: mcdespot2 [Output Prefix] [SPGR input file] [SPGR TR] [Number of phase cycling patterns] <[SSFP input file] [Phase Cycling]>  [SSFP TR] [Simplex/Region Contraction] [T1 Map File] [B1 Map File] <[Mask File]>\n\
\n\
The input file requires 2 columns - SSFP 180 file, flip angle.\n\
Specify the phase cycling pattern in degrees after each input file.\n";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	if (argc == 1)
	{
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	
	if (argc < 13)
	{
		fprintf(stderr, "Incorrect number of arguments (= %d) specified, should be at least 12.\n",
				argc - 1); // Subtract off filename
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}

	char *outPrefix = argv[1];
	char **spgrFilenames;
	double *spgrAngles;
	size_t nSPGR = readRecordFile(argv[2], "sd", &spgrFilenames, &spgrAngles);
	double spgrTR = atof(argv[3]);
	for (int i = 0; i < nSPGR; i++)
		spgrAngles[i] = radians(spgrAngles[i]);
	fprintf(stdout, "Specified %ld SPGR files with TR = %f ms.\n", nSPGR, spgrTR);
	
	size_t nPhases = atoi(argv[4]);
	size_t *nSSFP = malloc(nPhases * sizeof(size_t));
	char **ssfpFilenames[nPhases];
	double ssfpTR, *ssfpPhases = malloc(nPhases * sizeof(double)),
	       **ssfpAngles = malloc(nPhases * sizeof(double *));

	fprintf(stdout, "Specified %ld phase cycling patterns.\n", nPhases);	
	for (size_t p = 0; p < nPhases; p++)
	{
		nSSFP[p] = readRecordFile(argv[5 + (p * 2)], "sd", &(ssfpFilenames[p]), &(ssfpAngles[p]));
		ssfpPhases[p] = radians(atof(argv[6 + (p * 2)]));
		for (int i = 0; i < nSSFP[p]; i++)
			ssfpAngles[p][i] = radians(ssfpAngles[p][i]);
		fprintf(stdout, "Specified pattern %ld with phase %f deg and %ld files.\n", p, degrees(ssfpPhases[p]), nSSFP[p]);
		//ARR_D( ssfpAngles[p], nSSFP[p] );
	}
	ssfpTR = atof(argv[5 + (nPhases * 2)]);
	int mode = atoi(argv[6 + (nPhases * 2)]);
	fprintf(stdout, "Specified TR of %f ms.\n", ssfpTR);
	nifti_image *T1Map = NULL, *B1Map = NULL, *mask = NULL;
	fprintf(stdout, "Reading T1 Map: %s\n", argv[7 + (nPhases * 2)]);
	T1Map = nifti_image_read(argv[7 + (nPhases * 2)], FALSE);
	if (!T1Map)
		exit(EXIT_FAILURE);	
	fprintf(stdout, "Reading B1 Map: %s\n", argv[8 + (nPhases * 2)]);
	B1Map = nifti_image_read(argv[8 + (nPhases * 2)], FALSE);
	if (!B1Map)
		exit(EXIT_FAILURE);
	if (argc == (9 + (nPhases * 2)))
	{}
	else if (argc == (10 + (nPhases * 2)))
	{
		fprintf(stdout, "Reading Mask: %s\n", argv[9 + (nPhases * 2)]);
		mask = nifti_image_read(argv[9 + (nPhases * 2)], FALSE);
		if (!mask)
			exit(EXIT_FAILURE);
	}
	else
	{
		fprintf(stdout, "Argument count was %d, expected %ld or %ld\n", argc - 1, 9 + (nPhases * 2), 10 + (nPhases * 2));
		exit(EXIT_FAILURE);
	}
	//**************************************************************************	
	// Read in headers / Allocate memory for slices and results
	//**************************************************************************
	nifti_image **spgrHeaders = malloc(nSPGR * sizeof(nifti_image *));
	loadHeaders(spgrFilenames, spgrHeaders, nSPGR);
	nifti_image ***ssfpHeaders = malloc(nPhases * sizeof(nifti_image **));
	for (int p = 0; p < nPhases; p++)
	{
		ssfpHeaders[p] = (nifti_image **)malloc(nSSFP[p] * sizeof(nifti_image *));
		loadHeaders(ssfpFilenames[p], ssfpHeaders[p], nSSFP[p]);
		if ((p > 0) && ssfpHeaders[p - 1][0]->nvox != ssfpHeaders[p][0]->nvox)
		{
			fprintf(stderr, "Differing number of voxels in phase %d and %d headers.\n", p - 1, p);
			exit(EXIT_FAILURE);
		}
	}
	if (ssfpHeaders[0][0]->nvox != spgrHeaders[0]->nvox)
	{
		fprintf(stderr, "Differing number of voxels in SPGR and SSFP data.\n");
		exit(EXIT_FAILURE);
	}
	int voxelsPerSlice = spgrHeaders[0]->nx * spgrHeaders[0]->ny;	
	int totalVoxels = voxelsPerSlice * spgrHeaders[0]->nz;
		
	//**************************************************************************
	// Create results files
	// T1_m, T1_m, T2_f, T2_f,	f_m, tau_m, dw, residue
	// Need to write a full file of zeros first otherwise per-plane writing
	// won't produce a complete image.
	//**************************************************************************
	#define NR 8
	nifti_image **resultsHeaders = malloc(NR * sizeof(nifti_image *));
	char outName[strlen(outPrefix) + 12];
	float *blank = calloc(totalVoxels, sizeof(float));
	char *names[NR] = { "_T1_myel", "_T1_free", "_T2_myel", "_T2_free", "_frac_myel", "_tau_myel", "_dw", "_res" };
	float **resultsSlices = malloc(NR * sizeof(float*));
	for (int p = 0; p < NR; p++)
	{
		strcpy(outName, outPrefix); strcat(outName, names[p]);
		fprintf(stdout, "Writing blank result file:%s.\n", outName);
		resultsHeaders[p] = nifti_copy_nim_info(spgrHeaders[0]);
		nifti_set_filenames(resultsHeaders[p], outName, FALSE, TRUE);
		resultsHeaders[p]->data = blank;
		nifti_image_write(resultsHeaders[p]);
		resultsSlices[p] = malloc(voxelsPerSlice * sizeof(float));
	}

	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	fprintf(stdout, "Fitting mcDESPOT.\n");
	
	for (size_t slice = 0; slice < ssfpHeaders[0][0]->nz; slice++)
	{
		// Read in data
		fprintf(stdout, "Starting slice %ld...\n", slice);
		int sliceStart[NR] = {0, 0, (int)slice, 0, 0, 0, 0};
		int sliceDim[NR] = {spgrHeaders[0]->nx, spgrHeaders[0]->ny, 1, 1, 1, 1, 1};
		float **spgrData = malloc(nSPGR * sizeof(float *));
		for (int i = 0; i < nSPGR; i++)
		{
			spgrData[i] = malloc(voxelsPerSlice * sizeof(float));
			nifti_read_subregion_image(spgrHeaders[i], sliceStart, sliceDim, (void**)&(spgrData[i]));
		}
		
		float ***ssfpData = malloc(nPhases * sizeof(float **));
		for (int p = 0; p < nPhases; p++)
		{
			ssfpData[p] = malloc(nSSFP[p] * sizeof(float *));
			for (int i = 0; i < nSSFP[p]; i++)
			{
				ssfpData[p][i] = malloc(voxelsPerSlice * sizeof(float));
				nifti_read_subregion_image(ssfpHeaders[p][i], sliceStart, sliceDim, (void**)&(ssfpData[p][i]));
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
		
		__block bool hasVoxels = false;
		__block size_t voxCount = 0;
		time_t loopStart = time(NULL);
		//for (int vox = 0; vox < voxelsPerSlice; vox++)
		void (^processVoxel)(size_t vox) = ^(size_t vox)
		{
			double params[NR];
			arraySet(params, 0., NR);
			if (!mask || (maskData[vox] > 0.))
			{
				hasVoxels = true;
				voxCount++;
				double spgrSignal[nSPGR];
				for (int img = 0; img < nSPGR; img++)
					spgrSignal[img] = (double)spgrData[img][vox];
				
				double *ssfpSignal[nPhases];
				for (int p = 0; p < nPhases; p++)
				{
					ssfpSignal[p] = malloc(nSSFP[p] * sizeof(double));
					for (int img = 0; img < nSSFP[p]; img++)
						ssfpSignal[p][img] = (double)ssfpData[p][img][vox];
				}
				
				double T1 = (double)T1Data[vox];
				double B1 = (double)B1Data[vox];
				if (mode == 0)
				{
				}
				else
				{	// Use region contraction
					mcDESPOT(nSPGR, spgrAngles, spgrSignal, spgrTR,
					         nPhases, nSSFP, ssfpPhases, ssfpAngles, ssfpSignal, ssfpTR,
							 T1, B1, params);
				}
				for (size_t p = 0; p < nPhases; p++)
					free(ssfpSignal[p]);
			}
			for (int p = 0; p < NR; p++)
				resultsSlices[p][vox]  = (float)params[p];
		};
		dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
		dispatch_apply(voxelsPerSlice, global_queue, processVoxel);
		time_t loopEnd = time(NULL);
		if (hasVoxels)
		{
			fprintf(stdout, "Slice %ld had %ld voxels to process, average time per voxel was %f s. Writing to results files...", 
			        slice, voxCount, difftime(loopEnd, loopStart) / voxCount);
			for (int p = 0; p < NR; p++)
				nifti_write_subregion_image(resultsHeaders[p], sliceStart, sliceDim, (void **)&(resultsSlices[p]));
		}
		
		// Clean up memory
		for (int img = 0; img < nSPGR; img++)
			free(spgrData[img]);
		for (int p = 0; p < nPhases; p++)
		{
			for (int img = 0; img < nSSFP[p]; img++)
				free(ssfpData[p][img]);
			free(ssfpData[p]);
		}
		free(T1Data);
		free(B1Data);	
		fprintf(stdout, "Finished slice %ld.\n", slice);
	};

	fprintf(stdout, "All done.\n");
	return EXIT_SUCCESS;
}

