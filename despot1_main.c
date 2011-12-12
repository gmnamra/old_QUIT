/*
 *  despot1_main.c
 *  MacRI
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright 2011 Tobias Wood. All rights reserved.
 *
 */

#include <string.h>
#include <util.h>
#include <dispatch/dispatch.h>
#include "DESPOT1.h"
#include "nifti_tools.h"
#include "znzlib.h"

#define FALSE 0
#define TRUE 1

#define __DEBUG__ FALSE
#define __DEBUG_THRESH__ 60000

char *usage = "Usage is: despot1 [Output Prefix] [SPGR input file] [SPGR TR] <[SPGR-IR input file] [TR] [Flip Angle] [PE Readout Step count] [Inversion Mode]>\n\
\n\
Input File format:\n\
1st line - single integer (N) specifying number of files.\n\
N lines  - Path to file, space, parameter. For the SPGR files\n\
		    this is the flip angle in degrees, for the SPGR-IR files this\n\
		    is the TI in ms.\n\
\n\
Inversion Modes: 0 = 1.5T scanner (readout pulses div 2 + 2)\n\
                 1 = 3.0T scanner 1 (scale TI times by 0.9)\n\
                 2 = 3.0T scanner 2 (scale TI times by 0.84, readout pulses + 2)\n";

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//tests();
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	int nSPGR = 0, nIR = 0;
	double spgrTR; char **spgrFilenames;
	double *spgrAngles;
	double irTR, irAngle, *irTI, TIScale = 1.; int nReadout; char **irFilenames;
	char *outPrefix;
	
	if (argc == 1)
	{
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	
	if (!((argc - 1) == 3 || (argc - 1) == 8))
	{
		fprintf(stderr, "Incorrect number of arguments (= %d) specified, should be 3 for DESPOT1 and 8 for DESPOT1-HIFI.\n",
				argc - 1); // Subtract off filename
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	
	outPrefix = argv[1];
	nSPGR = readRecordFile(argv[2], "sf", &spgrFilenames, &spgrAngles);
	for (int i = 0; i < nSPGR; i++)
		spgrAngles[i] = radians(spgrAngles[i]);
	spgrTR = atof(argv[3]);
	fprintf(stdout, "Specified %d SPGR files with TR=%f ms.\n", nSPGR, spgrTR);
	
	if ((argc - 1) == 8)
	{
		nIR = readRecordFile(argv[4], "sf", &irFilenames, &irTI);
		irTR = atof(argv[5]);
		irAngle = radians(atof(argv[6]));
		nReadout = atoi(argv[7]);
		int invMode = atoi(argv[8]);
		switch (invMode)
		{
			case 0:
				TIScale = 1.0;
				nReadout = (nReadout / 2) + 2;
				break;
			case 1:
				TIScale = 0.9; // From Sean's code
				nReadout = (nReadout / 2) + 2;
				break;
			case 2:
				TIScale = 0.84; // From Sean's code
				nReadout = nReadout + 2;
			default:
				fprintf(stderr, "Inversion mode must be 0, 1, or 2\n");
				exit(EXIT_FAILURE);
				break;
		}
		for (int i = 0; i < nIR; i++)
			irTI[i] = irTI[i] * TIScale;
		fprintf(stdout, "Specified %d SPGR-IR files with TR=%f ms, flip angle: %f degrees and (calculated) readout pulse count: %d\n", nIR, irTR, degrees(irAngle), nReadout);
	}
	//**************************************************************************	
	// Read in headers / Allocate memory for slices and results
	//**************************************************************************
	nifti_image **SPGRFiles = (nifti_image **)malloc(nSPGR * sizeof(nifti_image *));
	nifti_image **irFiles = (nifti_image **)malloc(nIR * sizeof(nifti_image *));
	loadHeaders(spgrFilenames, SPGRFiles, nSPGR);
	loadHeaders(irFilenames, irFiles, nIR);
	
	int voxelsPerSlice = SPGRFiles[0]->nx * SPGRFiles[0]->ny;	
	int totalVoxels = voxelsPerSlice * SPGRFiles[0]->nz;
	__block float *T1Data = (float *)malloc(totalVoxels * sizeof(float));
	__block float *M0Data = (float *)malloc(totalVoxels * sizeof(float));
	__block float *B1Data = (float *)malloc(totalVoxels * sizeof(float));
	__block float *T1Smooth = (float *)malloc(totalVoxels * sizeof(float));
	__block float *M0Smooth = (float *)malloc(totalVoxels * sizeof(float));
	__block float *B1Smooth = (float *)malloc(totalVoxels * sizeof(float));
	int kernelSize = 15;
	__block float *gaussKernel = matrixGaussianf(kernelSize, kernelSize, 3., 3.);
	fprintf(stdout, "Smoothing kernel is a gaussian, %d pixels wide/high.\n", kernelSize);
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	if (nIR > 0)
		fprintf(stdout, "Fitting DESPOT1-HIFI.\n");
	else
		fprintf(stdout, "Fitting classic DESPOT1.\n");
	
	for (int slice = 0; slice < SPGRFiles[0]->nz; slice++)
	//void (^processSlice)(size_t slice) = ^(size_t slice)
	{
		// Read in data
		fprintf(stdout, "Processing slice %ld...\n", slice);
		//float T1 = 0., M0 = 0., B1 = 1.; // Place to restore per-voxel return values, assume B1 field is uniform for classic DESPOT
		double pars[3] = {0., 1., 0.};
		
		int sliceStart[7] = {0, 0, slice, 0, 0, 0, 0};
		int sliceDim[7] = {SPGRFiles[0]->nx, SPGRFiles[0]->ny, 1, 1, 1, 1, 1};
		float *SPGRData[nSPGR];
		for (int i = 0; i < nSPGR; i++)
			SPGRData[i] = (float *)malloc(voxelsPerSlice * sizeof(float));
		float *irData[nIR];
		for (int i = 0; i < nIR; i++)
			irData[i] = (float *)malloc(voxelsPerSlice * sizeof(float));
			
		for (int img = 0; img < nSPGR; img++)
			nifti_read_subregion_image(SPGRFiles[img], sliceStart, sliceDim, (void**)&(SPGRData[img]));
		for (int img = 0; img < nIR; img++)
			nifti_read_subregion_image(irFiles[img], sliceStart, sliceDim, (void**)&(irData[img]));

		int sliceIndex = slice * voxelsPerSlice;	
		for (int vox = 0; vox < voxelsPerSlice; vox++)
		{
			double spgrs[nSPGR];
			for (int img = 0; img < nSPGR; img++)
				spgrs[img] = (double)SPGRData[img][vox];		
			if (nIR > 0)
			{
				double irs[nIR];
				for (int img = 0; img < nIR; img++)
					irs[img] = (double)irData[img][vox];
				B1 = calcHIFI(spgrs, spgrAngles, nSPGR, spgrTR,
							  irs, irTI, nIR, irAngle, irTR, nReadout);
			}
			else
				calcDESPOT1(spgrs, spgrAngles, nSPGR, spgrTR, B1, &T1, &M0);
			
			//ARR_D(pars, 3);
			// Sanity check
			if (pars[1] < 0.) pars[1] = 0.; if (pars[1] > 2.) pars[1] = 2.;
			if (pars[0] < 0.) pars[0] = 0.; if (pars[0] > 99999999.) pars[0] = 99999999.;
			if (pars[2] < 0.) pars[2] = 0.; if (pars[2] > 9999.) pars[2] = 9999.;		 
			
			T1Data[sliceIndex + vox] = pars[2];
			M0Data[sliceIndex + vox] = pars[0];
			B1Data[sliceIndex + vox] = pars[1];
		}
		
		if (nIR > 0)
		{
			matrixConvolvef(&(B1Smooth[sliceIndex]), &(B1Data[sliceIndex]),
						   SPGRFiles[0]->nx, SPGRFiles[0]->ny,
						   gaussKernel, kernelSize, kernelSize);
			for (int vox = 0; vox < voxelsPerSlice; vox++)
			{
				double spgrs[nSPGR];
				for (int img = 0; img < nSPGR; img++)
					spgrs[img] = (double)SPGRData[img][vox];
				calcDESPOT1(spgrAngles, spgrs, nSPGR, spgrTR, &(pars[0]), &(pars[2]), pars[1]);
				
				// Sanity check
				if (pars[0] < 0.) pars[0] = 0.; if (pars[0] > 99999999.) pars[0] = 99999999.;
				if (pars[2] < 0.) pars[2] = 0.; if (pars[2] > 9999.) pars[2] = 9999.;				
				T1Data[sliceIndex + vox] = pars[2];
				M0Data[sliceIndex + vox] = pars[0];
			}
		};
		
		
		// Clean up memory
		for (int i = 0; i < nSPGR; i++)
			free(SPGRData[i]);
		for (int i = 0; i < nIR; i++)
			free(irData[i]);
	};
	//dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
	//dispatch_apply(SPGRFiles[0]->nz, global_queue, processSlice);
	fprintf(stdout, "Finished fitting. Writing results files.\n");

	//**************************************************************************
	// Create results files
	//**************************************************************************
	nifti_image *out = nifti_copy_nim_info(SPGRFiles[0]);
	char outName[strlen(outPrefix) + 4]; // Space for "T1" plus null
	strcpy(outName, outPrefix); strcat(outName, "_T1");
	fprintf(stdout, "Writing T1 Map: %s\n", outName);
	nifti_set_filenames(out, outName, FALSE, TRUE);	
	out->data = (void *)T1Data;
	nifti_image_write(out);
	strcpy(outName, outPrefix); strcat(outName, "_M0");
	fprintf(stdout, "Writing M0 Map: %s\n", outName);
	nifti_set_filenames(out, outName, FALSE, TRUE);	
	out->data = (void *)M0Data;
	nifti_image_write(out);
	if (nIR > 0)
	{
		strcpy(outName, outPrefix); strcat(outName, "_B1");
		fprintf(stdout, "Writing B1 Map: %s\n", outName);
		nifti_set_filenames(out, outName, FALSE, TRUE);
		out->data = (void *)B1Data;
		nifti_image_write(out);
		strcpy(outName, outPrefix); strcat(outName, "_SmoothT1");
		fprintf(stdout, "Writing Smoothed T1 Map: %s\n", outName);
		nifti_set_filenames(out, outName, FALSE, TRUE);	
		out->data = (void *)T1Smooth;
		nifti_image_write(out);
		strcpy(outName, outPrefix); strcat(outName, "_SmoothM0");
		fprintf(stdout, "Writing Smoothed M0 Map: %s\n", outName);
		nifti_set_filenames(out, outName, FALSE, TRUE);	
		out->data = (void *)M0Smooth;
		nifti_image_write(out);	
		strcpy(outName, outPrefix); strcat(outName, "_SmoothB1");
		fprintf(stdout, "Writing Smoothed B1 Map: %s\n", outName);
		nifti_set_filenames(out, outName, FALSE, TRUE);
		out->data = (void *)B1Smooth;
		nifti_image_write(out);
	}
	fprintf(stdout, "All done.\n");
	return EXIT_SUCCESS;
}
