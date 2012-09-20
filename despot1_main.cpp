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
#ifdef __APPLE__
	#include <libkern/OSAtomic.h>
	#define AtomicAdd OSAtomicAdd32
#else
	#define AtomicAdd(x, y) (*y) += x
#endif
#include "DESPOT.h"
#include "NiftiImage.h"
#include "procpar.h"
const char *usage = "Usage is: despot1 [options] spgr_input output_prefix \n\
\
Options:\n\
	-m, --mask file : Mask input with specified file.\n\
	--B1 file       : Correct flip angles with specified B1 ratio.\n\
	-v, --verbose   : Print out more messages.\n";
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//tests();
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	char procpar[MAXSTR];
	double spgrTR = 0.;
	int nSPGR;
	double *B1Data = NULL, *maskData = NULL;
	NiftiImage spgrFile, inFile;
	par_t *pars;
	static int verbose = false;
	static struct option long_options[] =
	{
		{"B1", required_argument, 0, '1'},
		{"mask", required_argument, 0, 'm'},
		{"verbose", no_argument, 0, 'v'},
		{0, 0, 0, 0}
	};
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "m:v", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case '1':
				fprintf(stdout, "Opening B1 file: %s\n", optarg);
				inFile.open(optarg, 'r');
				B1Data = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case 'm':
				fprintf(stdout, "Opening mask file: %s\n", optarg);
				inFile.open(optarg, 'r');
				maskData = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case 'v':
				verbose = true;
				break;
			case '?': // getopt will print an error message
				abort();
		}
	}
	if ((argc - optind) != 2)
	{
		fprintf(stderr, "Incorrect number of arguments.\n%s", usage);
		exit(EXIT_FAILURE);
	}
	fprintf(stdout, "Opening SPGR file: %s\n", argv[optind]);
	spgrFile.open(argv[optind], 'r');
	nSPGR = spgrFile.nt();
	snprintf(procpar, MAXSTR, "%s.procpar", spgrFile.basename().c_str());
	pars = readProcpar(procpar);
	VectorXd spgrAngles(nSPGR);
	if (pars)
	{
		spgrTR = realVal(pars, "tr", 0);
		for (int i = 0; i < nSPGR; i++) spgrAngles[i] = realVal(pars, "flip1", i);
		freeProcpar(pars);
	}
	else
	{
		fprintf(stdout, "Enter SPGR TR (s):");
		fscanf(stdin, "%lf", &spgrTR);
		fprintf(stdout, "Enter SPGR Flip Angles (degrees):");
		for (int i = 0; i < nSPGR; i++) std::cin >> spgrAngles[i];
	}
	spgrAngles *= M_PI / 180.;
	const std::string outPrefix(argv[++optind]);
	if (verbose)
	{
		std::cout << "SPGR TR=" << spgrTR
		          << " s. Flip-angles: " << spgrAngles.transpose() * 180. / M_PI << std::endl;
		std::cout << "Ouput prefix will be: " << outPrefix << std::endl;
	}
	//**************************************************************************
	// Allocate memory for slices
	//**************************************************************************	
	int voxelsPerSlice = spgrFile.nx() * spgrFile.ny();
	int totalVoxels = spgrFile.voxelsPerVolume();
	__block int voxCount;
	
	fprintf(stdout, "Reading SPGR data...\n");
	double *SPGR = spgrFile.readAllVolumes<double>();
	spgrFile.close();
	fprintf(stdout, "done.\n");
	//**************************************************************************
	// Create results data storage
	//**************************************************************************
	#define NR 3
	double **resultsData   = (double **)malloc(NR * sizeof(double *));
	for (int i = 0; i < NR; i++)
		resultsData[i] = (double *)malloc(totalVoxels * sizeof(double));
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
	dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
	for (int slice = 0; slice < spgrFile.nz(); slice++)
	{
		clock_t loopStart, loopEnd;
		// Read in data
		if (verbose)
			fprintf(stdout, "Processing slice %d...\n", slice);
		loopStart = clock();
		voxCount = 0;
		int sliceOffset = slice * voxelsPerSlice;
		//for (size_t vox = 0; vox < voxelsPerSlice; vox++)
		void (^processVoxel)(size_t vox) = ^(size_t vox)
		{
			double T1 = 0., M0 = 0., B1 = 1., res = 0.; // Place to restore per-voxel return values, assume B1 field is uniform for classic DESPOT
			if ((!maskData) || (maskData[sliceOffset + vox] > 0.))
			{
				AtomicAdd(1, &voxCount);
				if (B1Data)
					B1 = B1Data[sliceOffset + vox];
				ArrayXd spgrs(nSPGR);
				for (int img = 0; img < nSPGR; img++)
					spgrs[img] = SPGR[img * totalVoxels + sliceOffset + vox];
				res = classicDESPOT1(spgrAngles, spgrs, spgrTR, B1, &M0, &T1);
				
				// Sanity check
				M0 = clamp(M0, 0., 1.e7);
				T1 = clamp(T1, 0., 15.);
			}
			resultsData[0][sliceOffset + vox] = M0;
			resultsData[1][sliceOffset + vox] = T1;
			resultsData[2][sliceOffset + vox] = res;
		};
		dispatch_apply(voxelsPerSlice, global_queue, processVoxel);
		
		if (verbose)
		{
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
	}
	const char *names[NR] = { "_M0", "_T1", "_despot1_res" };
	NiftiImage outFile(spgrFile);
	spgrFile.setDatatype(DT_FLOAT32);
	spgrFile.setDims(spgrFile.nx(), spgrFile.ny(), spgrFile.nz(), 1);
	for (int r = 0; r < NR; r++)
	{
		std::string outName = outPrefix + names[r] + ".nii.gz";
		if (verbose)
			std::cout << "Writing result header: " << outName << std::endl;
		spgrFile.open(outName, 'w');
		spgrFile.writeVolume<double>(0, resultsData[r]);
		spgrFile.close();
		free(resultsData[r]);
	}
	// Clean up memory
	free(SPGR);
	free(B1Data);
	free(maskData);
	fprintf(stdout, "All done.\n");
	exit(EXIT_SUCCESS);
}
