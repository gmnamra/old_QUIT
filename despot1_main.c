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
#include <libkern/OSAtomic.h>
#include "DESPOT.h"
#include "fslio.h"

char *usage = "Usage is: despot1 [options] spgr_input output_prefix \n\
\
Options:\n\
	-m file  : Mask input with specified file.\n\
	-h file  : Use HIFI mode with specified SPGR-IR data.\n\
	-z       : Output .nii.gz files.\n\
	-n N     : Number of readout pulses for inversion mode.\n\
	-i 0-3   : Specify the scanner Inversion mode:\n\
	           0 (Default) Use raw segment TR from input file\n\
			   1 = 1.5T scanner, readout pulses div 2 + 2\n\
	           2 = 3T scanner, scale TI times by 0.9, readout pulses div 2 + 2\n\
	           3 = 3T scanner, scale TI times by 0.84, readout pulses + 2)\n\
			   1,2 & 3 MUST be followed by the PE Readout Step Count\n";
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//tests();
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	if (argc < 3)
	{
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	
	int nSPGR = 0, nIR = 0, nReadout = 0, invMode = 0;
	short nx, ny, nz, nv;
	char *outPrefix = NULL, *outExt = ".nii";
	double spgrTR = 0., *spgrAngles = NULL;
	double irTR = 0., irAngle = 0., *irTI = NULL, TIScale = 1.;
	FSLIO *spgrFile = NULL, *irFile = NULL, *maskFile = NULL;
	
	
	static struct option long_options[] =
	{
		{"B0", required_argument, 0, '0'},
		{"B1", required_argument, 0, '1'},
		{0, 0, 0, 0}
	};
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "h:i:m:n:z", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 'h':
				fprintf(stdout, "Opening IR file: %s\n", optarg);
				irFile = FslOpen(optarg, "rb");
				break;
			case 'i':
				invMode = atoi(optarg);
				break;
			case 'n':
				nReadout = atoi(optarg);
				break;
			case 'm':
				fprintf(stdout, "Opening mask file: %s\n", optarg);
				maskFile = FslOpen(optarg, "rb");
				break;
			case 'z':
				outExt = ".nii.gz";
				break;
		}
	}
	if (invMode && !nReadout)
	{
		fprintf(stderr, "Inversion mode was specified but not readout pulse count.\n");
		exit(EXIT_FAILURE);
	}
	
	switch (invMode)
	{
		case 1:
			TIScale = 1.0;
			nReadout = (nReadout / 2) + 2;
			break;
		case 2:
			TIScale = 0.9; // From Sean's code
			nReadout = (nReadout / 2) + 2;
			break;
		case 3:
			TIScale = 0.84; // From Sean's code
			nReadout = nReadout + 2;
			break;
		case 0:
			TIScale = 1.0;
			nReadout = 0;
			break;
		default:
			fprintf(stderr, "Inversion mode must be 0 - 3\n");
			exit(EXIT_FAILURE);
			break;
	}
	fprintf(stdout, "Opening SPGR file: %s\n", argv[optind]);
	spgrFile = FslOpen(argv[optind++], "rb");
	fprintf(stdout, "Ouput prefix will be: %s\n", argv[optind]);
	outPrefix = argv[optind++];
	if (irFile && !FslCheckDims(spgrFile, irFile))
	{
		fprintf(stderr, "Dimensions of SPGR & IR-SPGR files do not match.\n");
		exit(EXIT_FAILURE);
	}
	FslGetDim(spgrFile, &nx, &ny, &nz, &nv);
	nSPGR = nv;
	spgrAngles = malloc(nSPGR * sizeof(double));
	fprintf(stdout, "Enter SPGR TR (s):");
	fscanf(stdin, "%lf", &spgrTR);
	fprintf(stdout, "Enter SPGR Flip Angles (degrees):");
	fgetArray(stdin, 'd', nSPGR, spgrAngles);
	
	fprintf(stdout, "SPGR TR=%f ms. ", spgrTR);
	ARR_D(spgrAngles, nSPGR);
	for (int i = 0; i < nSPGR; i++)
		spgrAngles[i] = radians(spgrAngles[i]);

	//**************************************************************************	
	// Gather IR-SPGR Data
	//**************************************************************************	
	if (irFile)
	{
		nIR = irFile->niftiptr->nt;
		irTI = malloc(nIR * sizeof(double));
		fprintf(stdout, "Enter IR-SPGR Flip Angle (degrees):");
		fscanf(stdin, "%lf", &irAngle);
		irAngle = radians(irAngle);
		if (nReadout > 0)
		{
			fprintf(stdout, "Enter IR-SPGR TR (s):");
			fscanf(stdin, "%lf", &irTR);
			irTR = irTR * nReadout;
			fprintf(stdout, "Enter IR-SPGR TI times (s):");
			fgetArray(stdin, 'd', nIR, irTI);
			for (int i = 0; i < nIR; i++)
				irTI[i] = irTI[i] * TIScale;
		}
		else
		{
			fprintf(stdout, "Enter IR-SPGR TI times (s):");
			fgetArray(stdin, 'd', nIR, irTI);
			fprintf(stdout, "Enter first scan Segment TR (s):");
			fscanf(stdin, "%lf", &irTR);
			irTR -= irTI[0]; // Subtract off TI to get 
		}

		fprintf(stdout, "Specified %d SPGR-IR images with flip angle: %f degrees, TR = %f (s) ", nIR, degrees(irAngle), irTR);
		ARR_D(irTI, nIR);
	}

	//**************************************************************************	
	// Allocate memory for slices
	//**************************************************************************	
	int voxelsPerSlice = nx * ny;
	int totalVoxels = voxelsPerSlice * nz;
	__block int voxCount;
	float *SPGRData = malloc(nSPGR * voxelsPerSlice * sizeof(float));
	float *irData = NULL, *maskData = NULL;
	if (irFile)
		irData = malloc(nIR * voxelsPerSlice * sizeof(float));
	if (maskFile)
		maskData = malloc(voxelsPerSlice * sizeof(float));	
	//**************************************************************************
	// Create results headers
	//**************************************************************************
	#define NR 4
	FSLIO **resultsHeaders = malloc(NR * sizeof(FSLIO *));
	float **resultsData = malloc(NR * sizeof(float*));
	char *names[NR] = { "_M0", "_T1", "_B1", "_despot1_res" };
	char outName[strlen(outPrefix) + 64];
	for (int r = 0; r < NR; r++)
	{
		strcpy(outName, outPrefix); strcat(outName, names[r]); strcat(outName, outExt);
		fprintf(stdout, "Writing result header:%s.\n", outName);
		resultsHeaders[r] = FslOpen(outName, "wb");
		FslCloneHeader(resultsHeaders[r], spgrFile);
		FslSetDim(resultsHeaders[r], nx, ny, nz, 1);
		FslSetDataType(resultsHeaders[r], DTYPE_FLOAT);
		resultsData[r] = malloc(totalVoxels * sizeof(float));
	}
		
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	if (nIR == 0)
		fprintf(stdout, "Fitting classic DESPOT1.\n");
	else
		fprintf(stdout, "Fitting DESPOT1-HIFI.\n");
	
	for (int slice = 0; slice < nz; slice++)
	{
		clock_t loopStart, loopEnd;
		// Read in data
		fprintf(stdout, "Processing slice %d...\n", slice);
		FslReadSliceSeries(spgrFile, SPGRData, slice, nSPGR);
		if (irFile)
			FslReadSliceSeries(irFile, irData, slice, nIR);
		if (maskFile)
			FslReadSliceSeries(maskFile, maskData, slice, 1);
				
		loopStart = clock();
		voxCount = 0;
		int sliceIndex = slice * voxelsPerSlice;
		//for (size_t vox = 0; vox < voxelsPerSlice; vox++)
		void (^processVoxel)(size_t vox) = ^(size_t vox)
		{
			double T1 = 0., M0 = 0., B1 = 1., res = 0.; // Place to restore per-voxel return values, assume B1 field is uniform for classic DESPOT
			if ((!maskFile) || (maskData[vox] > 0.))
			{
				OSAtomicAdd32(1, &voxCount);
				double spgrs[nSPGR];
				for (int img = 0; img < nSPGR; img++)
					spgrs[img] = (double)SPGRData[voxelsPerSlice * img + vox];
				res = calcDESPOT1(spgrAngles, spgrs, nSPGR, spgrTR, B1, &M0, &T1);
				if (nIR > 0)
				{
					double irs[nIR];
					for (int img = 0; img < nIR; img++)
						irs[img] = (double)irData[voxelsPerSlice * img + vox];
					B1 = 1.;
					res = calcHIFI(spgrAngles, spgrs, nSPGR, spgrTR,
							       irTI, irs, nIR, irAngle, irTR,
							       &M0, &T1, &B1);
				}
				
				// Sanity check
				M0 = clamp(M0, 0., 1.e7);
				T1 = clamp(T1, 0., 5.);
				B1 = clamp(B1, 0., 2.);
			}
			
			resultsData[0][sliceIndex + vox] = (float)M0;
			resultsData[1][sliceIndex + vox] = (float)T1;
			resultsData[2][sliceIndex + vox] = (float)B1;
			resultsData[3][sliceIndex + vox] = (float)res;
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
		FslWriteVolumes(resultsHeaders[r], resultsData[r], 1);
		FslClose(resultsHeaders[r]);
	}
	// Clean up memory
	free(SPGRData);
	free(irData);
	fprintf(stdout, "All done.\n");
	exit(EXIT_SUCCESS);
}
