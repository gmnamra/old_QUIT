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
#include <dispatch/dispatch.h>
#include <libkern/OSAtomic.h>
#include "DESPOT.h"
#include "nifti_tools.h"
#include "znzlib.h"

char *usage = "Usage is: despot1 [options] output_prefix spgr_input \n\
\
Options:\n\
	-m file  : Mask input with specified file.\n\
	-h file  : Use HIFI mode with specified SPGR-IR data.\n\
	-z       : Output .nii.gz files.\n\
	-i 0-3 N : Specify the scanner Inversion mode:\n\
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
	
	int nSPGR = 0, nIR = 0, nReadout = 0;
	char *irFilename = NULL, *outPrefix = NULL, *outExt = ".nii";
	float spgrTR = 0., *spgrAngles = NULL;
	float irTR = 0., irAngle = 0., *irTI = NULL, TIScale = 1.;
	nifti_image *spgrFile = NULL, *irFile = NULL, *maskFile = NULL;
	
	int thisArg = 1;
	while ((thisArg < argc) && (argv[thisArg][0] =='-'))
	{
		switch (argv[thisArg][1])
		{
			case 'm':
				maskFile = nifti_image_read(argv[++thisArg], FALSE);
				break;
			case 'h':
				irFilename = argv[++thisArg];
				break;
			case 'z':
				outExt = ".nii.gz";
				break; 
			case 'i':
				switch (atoi(argv[++thisArg]))
				{
					case 1:
						TIScale = 1.0;
						nReadout = (atoi(argv[++thisArg]) / 2) + 2;
						break;
					case 2:
						TIScale = 0.9; // From Sean's code
						nReadout = (atoi(argv[++thisArg]) / 2) + 2;
						break;
					case 3:
						TIScale = 0.84; // From Sean's code
						nReadout = atoi(argv[++thisArg]) + 2;
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
			default:
				fprintf(stderr, "Undefined command line option\n%s", usage);
				exit(EXIT_FAILURE);
				break;
		}
		++thisArg;
	}
	outPrefix = argv[thisArg]; thisArg++;
	fprintf(stdout, "Reading headers.\n");
	spgrFile = nifti_image_read(argv[thisArg], FALSE);
	nSPGR = spgrFile->nt;
	spgrAngles = malloc(nSPGR * sizeof(float));
	fprintf(stdout, "Enter SPGR TR (ms):");
	fscanf(stdin, "%f", &spgrTR);
	fprintf(stdout, "Enter SPGR Flip Angles (degrees):");
	fgetArray(stdin, 'f', nSPGR, spgrAngles);
	
	fprintf(stdout, "SPGR TR=%f ms. ", spgrTR);
	ARR_D(spgrAngles, nSPGR);
	for (int i = 0; i < nSPGR; i++)
		spgrAngles[i] = radians(spgrAngles[i]);

	//**************************************************************************	
	// Gather IR-SPGR Data
	//**************************************************************************	
	if (irFilename)
	{
		irFile = nifti_image_read(irFilename, FALSE);
		nIR = irFile->nt;
		irTI = malloc(nIR * sizeof(float));
		fprintf(stdout, "Enter IR-SPGR Flip Angle (degrees):");
		fscanf(stdin, "%f", &irAngle);
		irAngle = radians(irAngle);
		if (nReadout > 0)
		{
			fprintf(stdout, "Enter IR-SPGR TR (ms):");
			fscanf(stdin, "%f", &irTR);
			irTR = irTR * nReadout;
			fprintf(stdout, "Enter IR-SPGR TI times (ms):");
			fgetArray(stdin, 'f', nIR, irTI);
			for (int i = 0; i < nIR; i++)
				irTI[i] = irTI[i] * TIScale;
		}
		else
		{
			fprintf(stdout, "Enter IR-SPGR TI times (ms):");
			fgetArray(stdin, 'f', nIR, irTI);
			fprintf(stdout, "Enter first scan Segment TR (ms):");
			fscanf(stdin, "%f", &irTR);
			irTR -= irTI[0]; // Subtract off TI to get 
		}

		fprintf(stdout, "Specified %d SPGR-IR images with flip angle: %f degrees, TR = %f (ms) ", nIR, degrees(irAngle), irTR);
		ARR_D(irTI, nIR);
	}

	//**************************************************************************	
	// Allocate memory for slices and results
	//**************************************************************************	
	int voxelsPerSlice = spgrFile->nx * spgrFile->ny;	
	int totalVoxels = voxelsPerSlice * spgrFile->nz;
	__block int voxCount;
	float *SPGRData = malloc(nSPGR * voxelsPerSlice * sizeof(float));
	float *irData = NULL, *maskData = NULL;
	if (irFile)
		irData = malloc(nIR * voxelsPerSlice * sizeof(float));
	if (maskFile)
		maskData = malloc(voxelsPerSlice * sizeof(float));	
	__block float *T1Data = (float *)malloc(totalVoxels * sizeof(float));
	__block float *M0Data = (float *)malloc(totalVoxels * sizeof(float));
	__block float *B1Data = (float *)malloc(totalVoxels * sizeof(float));
	__block float *T1Smooth = (float *)malloc(totalVoxels * sizeof(float));
	__block float *M0Smooth = (float *)malloc(totalVoxels * sizeof(float));
	__block float *B1Smooth = (float *)malloc(totalVoxels * sizeof(float));
	int kernelSize = 15;
	float *gaussKernel = matrixGaussianf(kernelSize, kernelSize, 3., 3.);
	fprintf(stdout, "Smoothing kernel is a gaussian, %d pixels wide/high.\n", kernelSize);
	
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	if (nIR == 0)
		fprintf(stdout, "Fitting classic DESPOT1.\n");
	else
	{
		if (totalVoxels != (irFile->nx * irFile->ny * irFile->nz))
		{
			fprintf(stderr, "SPGR and IR-SPGR files have different number of voxels (%ld and %ld)\n", spgrFile->nvox, irFile->nvox);
			exit(EXIT_FAILURE);
		}
		fprintf(stdout, "Fitting DESPOT1-HIFI.\n");
	}

	for (int slice = 0; slice < spgrFile->nz; slice++)
	{
		clock_t loopStart, loopEnd;
		// Read in data
		fprintf(stdout, "Processing slice %d...\n", slice);
		int sliceStart[7] = {0, 0, slice, 0, 0, 0, 0};
		int sliceDim[7] = {spgrFile->nx, spgrFile->ny, 1, 1, 1, 1, 1};
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
		nifti_read_subregion_image(spgrFile, sliceStart, sliceDim, (void**)&(SPGRData));
		if (nIR > 0)
		{
			sliceDim[3] = nIR;
			nifti_read_subregion_image(irFile, sliceStart, sliceDim, (void**)&(irData));
		}
		
		loopStart = clock();
		voxCount = 0;
		int sliceIndex = slice * voxelsPerSlice;
		void (^processVoxel)(size_t vox) = ^(size_t vox)
		{
			float T1 = 0., M0 = 0., B1 = 1.; // Place to restore per-voxel return values, assume B1 field is uniform for classic DESPOT
			if ((!maskFile) || (maskData[vox] > 0.))
			{
				OSAtomicAdd32(1, &voxCount);
				float spgrs[nSPGR];
				for (int img = 0; img < nSPGR; img++)
					spgrs[img] = (float)SPGRData[voxelsPerSlice * img + vox];		
				calcDESPOT1(spgrAngles, spgrs, nSPGR, spgrTR, B1, &M0, &T1);
				if (nIR > 0)
				{
					float irs[nIR];
					for (int img = 0; img < nIR; img++)
						irs[img] = (float)irData[voxelsPerSlice * img + vox];
					B1 = 1.;
					calcHIFI(spgrAngles, spgrs, nSPGR, spgrTR,
							 irTI, irs, nIR, irAngle, irTR,
							 &M0, &T1, &B1);
				}
				
				// Sanity check
				M0 = clamp(M0, 0., 1.e8);
				T1 = clamp(T1, 0., 3.e3);
				B1 = clamp(B1, 0., 2.);
			}
			T1Data[sliceIndex + vox] = (float)T1;
			M0Data[sliceIndex + vox] = (float)M0;
			B1Data[sliceIndex + vox] = (float)B1;
		};
		dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
		dispatch_apply(voxelsPerSlice, global_queue, processVoxel);
		
        loopEnd = clock();
        fprintf(stdout, "Finished first pass through slice %d", slice);
		if (voxCount)
		{
			fprintf(stdout, ", had %d unmasked voxels, CPU time per voxel was %f s.\n", 
			        voxCount, (loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC));
		}
		else
			fprintf(stdout, ", no unmasked voxels.\n");
		
		if (nIR > 0)
		{
			fprintf(stdout, "Smoothing...");
			matrixConvolvef(&(B1Smooth[sliceIndex]), &(B1Data[sliceIndex]),
						   spgrFile->nx, spgrFile->ny,
						   gaussKernel, kernelSize, kernelSize);
			fprintf(stdout, "done. Re-fitting with smoothed B1 parameter.\n");
			loopStart = clock();
			voxCount = 0;
			void (^processVoxel)(size_t vox) = ^(size_t vox)
			{
				float T1 = 0., M0 = 0., B1 = 1.; // Place to restore per-voxel return values, assume B1 field is uniform for classic DESPOT
				if (!maskFile || (maskData[vox] > 0))
				{
					OSAtomicAdd32(1, &voxCount);
					float spgrs[nSPGR];
					for (int img = 0; img < nSPGR; img++)
						spgrs[img] = (float)SPGRData[voxelsPerSlice * img + vox];
					B1 = (float)B1Smooth[sliceIndex + vox];
					calcDESPOT1(spgrAngles, spgrs, nSPGR, spgrTR, B1, &M0, &T1);
					// Sanity check
					M0 = clamp(M0, 0, 1.e8);
					T1 = clamp(T1, 0, 3.e3);
				}
				else
				{
					M0 = 0.; T1 = 0.; B1 = 0.;
				}

				T1Smooth[sliceIndex + vox] = T1;
				M0Smooth[sliceIndex + vox] = M0;
			};
			dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
			dispatch_apply(voxelsPerSlice, global_queue, processVoxel);
			
			loopEnd = clock();
	        fprintf(stdout, "Finished second pass through slice %d", slice);
			if (voxCount)
			{
				fprintf(stdout, ", had %d unmasked voxels, CPU time per voxel was %f s.\n", 
			    	    voxCount, (loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC));
			}
			else
				fprintf(stdout, ", no unmasked voxels.\n");
			}
	}
	// Clean up memory
	free(SPGRData);
	free(irData);

	//**************************************************************************
	// Create results files
	//**************************************************************************
	fprintf(stdout, "Finished fitting. Writing results files.\n");
	int outDims[8] = {3, spgrFile->nx, spgrFile->ny, spgrFile->nz, 1, 1, 1, 1};
	nifti_image *out = nifti_copy_orientation(spgrFile, outDims, DT_FLOAT);
	char outName[strlen(outPrefix) + 4]; // Space for "T1" plus null
	strcpy(outName, outPrefix); strcat(outName, "_M0"); strcat(outName, outExt);
	writeResult(out, outName, (void*)M0Data);
	strcpy(outName, outPrefix); strcat(outName, "_T1"); strcat(outName, outExt);
	writeResult(out, outName, (void*)T1Data);
	if (nIR > 0)
	{
		strcpy(outName, outPrefix); strcat(outName, "_B1"); strcat(outName, outExt);
		writeResult(out, outName, (void*)B1Data);
		strcpy(outName, outPrefix); strcat(outName, "_SmoothB1"); strcat(outName, outExt);
		writeResult(out, outName, (void*)B1Smooth);		
		strcpy(outName, outPrefix); strcat(outName, "_SmoothM0"); strcat(outName, outExt);
		writeResult(out, outName, (void*)M0Smooth);
		strcpy(outName, outPrefix); strcat(outName, "_SmoothT1"); strcat(outName, outExt);
		writeResult(out, outName, (void*)T1Smooth);
	}
	fprintf(stdout, "All done.\n");
	exit(EXIT_SUCCESS);
}
