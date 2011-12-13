/*
 *  inversion_main.c
 *  MacRI
 *
 *  Created by Tobias Wood on 27/10/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <string.h>
#include <util.h>
#include <dispatch/dispatch.h>
#include "nifti_tools.h"
#include "znzlib.h"
#include "recovery.h"
#include "mathUtil.h"

int main(int argc, char **argv)
{
	if (argc != 3)
	{
		fprintf(stderr, "Incorrect number of arguments.\n");
		exit(EXIT_FAILURE);
	}
	
	char *outPrefix = argv[1];
	char **imageNames;
	double *TRs, *scales;
	int nImgs = readRecordFile(argv[2], "sdd", &imageNames, &TRs, &scales);
	nifti_image **headers = (nifti_image **)malloc(nImgs * sizeof(nifti_image *));
	loadHeaders(imageNames, headers, nImgs);
	
	int voxelsPerSlice = headers[0]->nx * headers[0]->ny;	
	int totalVoxels = voxelsPerSlice * headers[0]->nz;
	__block float *M0Data = (float *)malloc(totalVoxels * sizeof(float));
	__block float *T1Data = (float *)malloc(totalVoxels * sizeof(float));
	__block float *alphaData = (float *)malloc(totalVoxels * sizeof(float));
	__block float *resData = (float *)malloc(totalVoxels * sizeof(float));
	// Do fitting
	//for (int slice = 0; slice < headers[0]->nz; slice++)
	void (^processSlice)(size_t slice) = ^(size_t slice)
	{
		int sliceStart[7] = {0, 0, slice, 0, 0, 0, 0};
		int sliceDim[7] = {headers[0]->nx, headers[0]->ny, 1, 1, 1, 1, 1};
		float **data = (float **)malloc(nImgs * sizeof(float *));
		for (int i = 0; i < nImgs; i++)
			data[i] = (float *)malloc(voxelsPerSlice * sizeof(float));
		for (int img = 0; img < nImgs; img++)
			nifti_read_subregion_image(headers[img], sliceStart, sliceDim, (void**)&(data[img]));
		
		int sliceIndex = slice * voxelsPerSlice;
		fprintf(stdout, "Processing slice %ld...\n", slice);
		for (int vox = 0; vox < voxelsPerSlice; vox++)
		{
			double vals[nImgs];
			double M0 = 0., T1 = 0., alpha = 0., residual = 0.;			
			for (int img = 0; img < nImgs; img++)
			{
				float d = data[img][vox];
				double s = scales[img];
				double r = d * s;
				vals[img] = r;
			}
			
			if (vals[nImgs - 1] > 3000.)
			{
				residual = calcRecovery(vals, TRs, nImgs, &M0, &T1, &alpha);
				if (M0 < 0.)
					M0 = 0.;
				if (M0 > 1000000.) M0 = 1000000.;
				if (T1 < 0.)
					T1 = 0.;
				if (T1 > 5000.) T1 = 5000.;
			}
			else
			{
				M0 = 0.;
				T1 = 0.;
			}
			
			M0Data[sliceIndex + vox] = (float)M0;
			T1Data[sliceIndex + vox] = (float)T1;
			alphaData[sliceIndex + vox] = (float)alpha;
			resData[sliceIndex + vox] = (float)residual;
		}
		for (int i = 0; i < nImgs; i++)
			free(data[i]);
		free(data);
	};
	dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
	dispatch_apply(headers[0]->nz, global_queue, processSlice);	
	fprintf(stdout, "Finished fitting. Writing results files.\n");
	nifti_image *out = nifti_copy_nim_info(headers[0]);
	char outName[strlen(outPrefix) + 4]; // Space for "T1" plus null
	strcpy(outName, outPrefix); strcat(outName, "_T1");
	writeResult(out, outName, (void *)T1Data);
	strcpy(outName, outPrefix); strcat(outName, "_M0");
	writeResult(out, outName, (void *)M0Data);
	strcpy(outName, outPrefix); strcat(outName, "_alpha");
	writeResult(out, outName, (void *)alphaData);
	strcpy(outName, outPrefix); strcat(outName, "_residual");
	writeResult(out, outName, (void *)resData);
	
	fprintf(stdout, "All done.");
	return EXIT_SUCCESS;
}


