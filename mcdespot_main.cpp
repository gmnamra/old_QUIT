/*
 *  mcdespot_main.c
 *  Fitting
 *
 *  Created by Tobias Wood on 14/02/2012.
 *  Copyright 2012 Tobias Wood. All rights reserved.
 *
 */

#include <time.h>
#include <stdbool.h>
#include <dispatch/dispatch.h>
#include <getopt.h>
#include <signal.h>

#ifdef __APPLE__
	#include <libkern/OSAtomic.h>
	#define AtomicAdd OSAtomicAdd32
#else
	#define AtomicAdd(x, y) (*y)++
#endif

#include "DESPOT_Functors.h"
#include "RegionContraction.h"
#include "fslio.h"
#include "procpar.h"

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const char *usage = "Usage is: mcdespot [options] output_prefix spgr_file ssfp_file1 (ssfp_fileN)\n\
\
Options:\n\
	-v, --verbose     : Print extra information.\n\
	-m, --mask file   : Mask input with specified file.\n\
	--2, --3          : Use a 2 or 3 component model (default 2).\n\
	--M0 file         : M0 Map file.\n\
	--B0 file         : B0 Map file.\n\
	--B1 file         : B1 Map file.\n\
	--start_slice n   : Only start processing at slice n.\n\
	--end_slice n     : Finish at slice n-1.\n\
	-s, --samples n   : Use n samples for region contraction (Default 5000).\n\
	-r, --retain  n   : Retain n samples for new boundary (Default 50).\n\
	-c, --contract n  : Contract a maximum of n times (Default 10).\n\
	-e, --expand n    : Re-expand boundary by percentage n (Default 0).\n\
	-b 3              : Boundaries suitable for 3T\n\
	   7              : Boundaries suitable for 7T (default)\n\
	   u              : User specified boundaries from stdin.\n";

static int verbose = false, start_slice = -1, end_slice = -1,
		   samples = 2500, retain = 25, contract = 10, components = 2, tesla = 7, nP = 0;
static double expand = 0.;
static std::string outPrefix;
static struct option long_options[] =
{
	{"B0", required_argument, 0, '0'},
	{"B1", required_argument, 0, '1'},
	{"M0", required_argument, 0, 'M'},
	{"mask", required_argument, 0, 'm'},
	{"verbose", no_argument, 0, 'v'},
	{"start_slice", required_argument, 0, 'S'},
	{"end_slice", required_argument, 0, 'E'},
	{"samples", required_argument, 0, 's'},
	{"retain", required_argument, 0, 'r'},
	{"contract", required_argument, 0, 'c'},
	{"expand", required_argument, 0, 'e'},
	{"2", no_argument, &components, 2},
	{"3", no_argument, &components, 3},
	{0, 0, 0, 0}
};
//******************************************************************************
// SIGTERM interrupt handler - for ensuring data gets saved even on a ctrl-c
//******************************************************************************
FSLIO *savedHeader;
double **paramsData;
double *residualData;

void write_results();
void write_results()
{
	std::string outPath;
	short nx, ny, nz, nvol;
	for (int p = 0; p < nP; p++)
	{
		outPath = outPrefix + "_";
		if (components == 2)
			outPath += TwoComponent::names[p];
		else
			outPath += ThreeComponent::names[p];
		outPath += ".nii.gz";
		std::cout << "Writing parameter file: " << outPath << std::endl;
		FSLIO *outFile = FslOpen(outPath.c_str(), "wb");
		FslCloneHeader(outFile, savedHeader);
		FslGetDim(outFile, &nx, &ny, &nz, &nvol);
		FslSetDim(outFile, nx, ny, nz, 1);
		FslSetDimensionality(outFile, 3);
		FslSetDataType(outFile, NIFTI_TYPE_FLOAT32);
		FslWriteHeader(outFile);
		FslWriteVolumeFromDouble(outFile, paramsData[p], 0);
		FslClose(outFile);
		free(paramsData[p]);
	}
	free(paramsData);
	
	outPath = outPrefix + "_residual.nii.gz";
	std::cout << "Writing residual file: " << outPath << std::endl;
	FSLIO *outFile = FslOpen(outPath.c_str(), "wb");
	FslCloneHeader(outFile, savedHeader);
	FslGetDim(outFile, &nx, &ny, &nz, &nvol);
	FslSetDim(outFile, nx, ny, nz, 1);
	FslSetDimensionality(outFile, 3);
	FslSetDataType(outFile, NIFTI_TYPE_FLOAT32);
	FslWriteHeader(outFile);
	FslWriteVolumeFromDouble(outFile, residualData, 0);
	FslClose(outFile);
	FslClose(savedHeader);
	free(residualData);
}

void int_handler(int sig);
void int_handler(int sig)
{
	fprintf(stdout, "Processing terminated. Writing currently processed data.\n");
	write_results();
	exit(EXIT_FAILURE);
}


//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	if (argc < 4)
	{
		fprintf(stderr, "%s", usage);
		exit(EXIT_FAILURE);
	}
	Eigen::initParallel();
	char procpar[MAXSTR];
	double spgrTR, ssfpTR,
	       *maskData = NULL, *M0Data = NULL, *B0Data = NULL, *B1Data = NULL;
	FSLIO *inFile = NULL;
	short nx, ny, nz, nSPGR, nPhases, nSSFP;
	par_t *pars;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "b:m:vzs:r:c:e:", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 'm':
				fprintf(stdout, "Reading mask file %s.\n", optarg);
				inFile = FslOpen(optarg, "rb");
				maskData = FslGetVolumeAsScaledDouble(inFile, 0);
				FslClose(inFile);
				break;
			case '0':
				fprintf(stdout, "Reading B0 file %s.\n", optarg);
				inFile = FslOpen(optarg, "rb");
				B0Data = FslGetVolumeAsScaledDouble(inFile, 0);
				FslClose(inFile);
				break;
			case '1':
				fprintf(stdout, "Reading B1 file %s.\n", optarg);
				inFile = FslOpen(optarg, "rb");
				B1Data = FslGetVolumeAsScaledDouble(inFile, 0);
				FslClose(inFile);
				break;
			case 'M':
				fprintf(stdout, "Reading M0 file %s.\n", optarg);
				inFile = FslOpen(optarg, "rb");
				M0Data = FslGetVolumeAsScaledDouble(inFile, 0);
				FslClose(inFile);
				break;
			case 'v': verbose = true; break;
			case 'S': start_slice = atoi(optarg); break;
			case 'E': end_slice = atoi(optarg); break;
			case 's': samples  = atoi(optarg); break;
			case 'r': retain   = atoi(optarg); break;
			case 'c': contract = atoi(optarg); break;
			case 'e': expand   = atof(optarg); break; 
			case 'b':
				switch (*optarg)
				{
					case '3':
						std::cout << "Using 3T boundaries.\n";
						tesla = 3;
						break;
					case '7':
						std::cout << "Using 7T boundaries.\n";
						tesla = 7;
						break;
					case 'u':
						tesla = -1;
						break;
					default:
						std::cout << "Unknown boundaries type " << optarg << std::endl;
						exit(EXIT_FAILURE);
						break;
				}
				break;
			case 0:
				// Just a flag
				break;
		}
	}
	outPrefix = argv[optind++];
	if (verbose)
		std::cout << "Output prefix will be: " << outPrefix << std::endl;
	//**************************************************************************
	// Gather SPGR Data
	//**************************************************************************
	std::cout << "Opening SPGR file: " << argv[optind] << std::endl;
	inFile = FslOpen(argv[optind], "rb");
	FslGetDim(inFile, &nx, &ny, &nz, &nSPGR);
	VectorXd spgrAngles(nSPGR);
	snprintf(procpar, MAXSTR, "%s.procpar", argv[optind]);
	pars = readProcpar(procpar);
	if (pars)
	{
		spgrTR = realVal(pars, "tr", 0);
		for (int i = 0; i < nSPGR; i++) spgrAngles[i] = realVal(pars, "flip1", i);
		freeProcpar(pars);
	}
	else
	{
		std::cout << "Enter SPGR TR (s): " << std::flush;
		std::cin >> spgrTR;
		std::cout << "Enter SPGR Flip Angles (degrees): " << std::flush;
		for (int i = 0; i < nSPGR; i++) std::cin >> spgrAngles[i];
	}
	spgrAngles *= M_PI / 180.;
	int voxelsPerSlice = nx * ny;
	int totalVoxels = voxelsPerSlice * nz;
	std::cout << "Reading SPGR data..." << std::endl;
	double *SPGR = FslGetAllVolumesAsScaledDouble(inFile);	
	// Save this header to output the results files.
	savedHeader = inFile;
	optind++;
	//**************************************************************************
	// Gather SSFP Data
	//**************************************************************************
	nPhases = argc - optind;
	VectorXd ssfpPhases(nPhases), ssfpAngles;
	double **SSFP = (double **)malloc(nPhases * sizeof(double *));
	for (size_t p = 0; p < nPhases; p++)
	{
		short inX, inY, inZ, inVols;
		std::cout << "Reading SSFP header from " << argv[optind] << std::endl;
		inFile = FslOpen(argv[optind], "rb");
		FslGetDim(inFile, &inX, &inY, &inZ, &inVols);
		eigen_assert((inX == nx) && (inY == ny) && (inZ == nz));
		snprintf(procpar, MAXSTR, "%s.procpar", argv[optind]);
		pars = readProcpar(procpar);		
		if (p == 0)
		{	// Read nSSFP, TR and flip angles from first file
			nSSFP = inVols;
			ssfpAngles.resize(nSSFP, 1);
			if (pars)
			{
				ssfpTR = realVal(pars, "tr", 0);
				for (int i = 0; i < nSSFP; i++)
					ssfpAngles[i] = realVal(pars, "flip1", i);
			}
			else
			{
				std::cout << "Enter SSFP TR (s): " << std::flush;
				std::cin >> ssfpTR;
				std::cout << "Enter " << nSSFP << " flip angles (degrees): " << std::flush;
				for (int i = 0; i < ssfpAngles.size(); i++)
					std::cin >> ssfpAngles[i];
			}
		}
		
		eigen_assert((inVols == nSSFP));
		if (pars)
		{
			ssfpPhases[p] = realVal(pars, "rfphase", 0);
			freeProcpar(pars);
		}
		else
		{
			std::cout << "Enter phase-cycling (degrees): " << std::flush;
			std::cin >> ssfpPhases[p];
		}
		std::cout << "Reading SSFP data..." << std::endl;
		SSFP[p] = FslGetAllVolumesAsScaledDouble(inFile);
		FslClose(inFile);
		optind++;
	}
	ssfpAngles *= M_PI / 180.;
	ssfpPhases *= M_PI / 180.;
	
	if (optind != argc)
	{
		std::cerr << "Unprocessed arguments supplied.\n" << usage;
		exit(EXIT_FAILURE);
	}	
	//**************************************************************************
	// Allocate results memory and set up boundaries
	//**************************************************************************
	if (components == 2)
	{
		std::cout << "Using 2 component model." << std::endl;
		nP = TwoComponent::nP;
	} else {
		std::cout << "Using 3 component model." << std::endl;
		nP = ThreeComponent::nP;
	}
		
	VectorXd loBounds(nP), hiBounds(nP);
	VectorXi loConstraints(nP), hiConstraints(nP);
	
	residualData = (double *)malloc(totalVoxels * sizeof(double));
	paramsData = (double **)malloc(nP * sizeof(double *));
	if (tesla < 0)
		std::cout << "Enter " << nP << " parameter pairs (low then high): " << std::flush;
	for (int i = 0; i < nP; i++)
	{
		paramsData[i] = (double *)malloc(totalVoxels * sizeof(double));
		loConstraints[i] = true; hiConstraints[i] = true;
		if (tesla == 3)
		{
			if (components == 2)
			{
				loBounds[i] = TwoComponent::lo3Bounds[i];
				hiBounds[i] = TwoComponent::hi3Bounds[i];
			} else {
				loBounds[i] = ThreeComponent::lo3Bounds[i];
				hiBounds[i] = ThreeComponent::hi3Bounds[i];
			}
		} else if (tesla == 7) {
			if (components == 2)
			{
				loBounds[i] = TwoComponent::lo7Bounds[i];
				hiBounds[i] = TwoComponent::hi7Bounds[i];
			} else {
				loBounds[i] = ThreeComponent::lo7Bounds[i];
				hiBounds[i] = ThreeComponent::hi7Bounds[i];
			}
		} else {
			std::cin >> loBounds[i] >> hiBounds[i];
		}
	}
	if (B0Data)
	{	// User is supplying a B0 estimate
		loBounds(nP - 1) = 0.;
		hiBounds(nP - 1) = 0.;
	}
	else
	{	// Make sure the B0 bounds are sensible
		loBounds(nP - 1) = (-0.5 / ssfpTR);
		hiBounds(nP - 1) = (0.5 / ssfpTR);
		loConstraints(nP - 1) = false;
		hiConstraints(nP - 1) = false;
	}
	if (verbose)
	{
		std::cout << "SPGR TR (s): " << spgrTR << " SPGR flip angles (deg): " << spgrAngles.transpose() * 180. / M_PI << std::endl;		
		std::cout << "SSFP TR (s): " << ssfpTR << " SSFP flip angles (deg): " << ssfpAngles.transpose() * 180. / M_PI << std::endl;	
		std::cout << "Phase Cycling Patterns (degrees): " << ssfpPhases.transpose() * 180. / M_PI << std::endl;
		std::cout << "Low bounds: " << loBounds.transpose() << std::endl;
		std::cout << "Hi bounds:  " << hiBounds.transpose() << std::endl;
	}	
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
    time_t procStart = time(NULL);
	if ((start_slice < 0) || (start_slice >= nz))
		start_slice = 0;
	if ((end_slice < 0) || (end_slice > nz))
		end_slice = nz;
	signal(SIGINT, int_handler);	// If we've got here there's actually allocated data to save
	std::cout << "Starting processing." << std::endl;
	dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);	
	for (size_t slice = start_slice; slice < end_slice; slice++)
	{
		if (verbose)
			std::cout << "Starting slice " << slice << "..." << std::flush;
		__block int voxCount = 0;
		const int sliceOffset = slice * voxelsPerSlice;
		clock_t loopStart = clock();
		//for (int vox = 0; vox < voxelsPerSlice; vox++)
		void (^processVoxel)(size_t vox) = ^(size_t vox)
		{
			double M0 = 1, B0 = INFINITY, B1 = 1., residual = 0.;
			if (M0Data) M0 = (double)M0Data[sliceOffset + vox];
			if (B0Data) B0 = (double)B0Data[sliceOffset + vox];
			if (B1Data) B1 = (double)B1Data[sliceOffset + vox];		
			VectorXd params(nP); params.setZero();
			if (!maskData || (maskData[sliceOffset + vox] > 0.))
			{
				AtomicAdd(1, &voxCount);
				VectorXd SPGR_signal(nSPGR);
				for (int vol = 0; vol < nSPGR; vol++)
					SPGR_signal[vol] = SPGR[totalVoxels*vol + sliceOffset + vox];
				std::vector<VectorXd> SSFP_signals;
				for (int p = 0; p < nPhases; p++)
				{
					VectorXd temp(nSSFP);
					for (int vol = 0; vol < nSSFP; vol++)
						temp[vol] = SSFP[p][totalVoxels*vol + sliceOffset + vox];
					SSFP_signals.push_back(temp);
				}
				if (components == 2)
				{
					TwoComponent tc(spgrAngles, SPGR_signal, ssfpAngles, ssfpPhases, SSFP_signals, spgrTR, ssfpTR, M0, B0, B1);
					residual = regionContraction<TwoComponent>(params, tc, loBounds, hiBounds,
															   loConstraints, hiConstraints,
															   samples, retain, contract, 0.05, expand, vox + time(NULL));
				} else {
					ThreeComponent tc(spgrAngles, SPGR_signal, ssfpAngles, ssfpPhases, SSFP_signals, spgrTR, ssfpTR, M0, B0, B1);
					residual = regionContraction<ThreeComponent>(params, tc, loBounds, hiBounds,
																 loConstraints, hiConstraints,
																 samples, retain, contract, 0.05, expand, vox + time(NULL));
				}
			}
			for (int p = 0; p < nP; p++)
				paramsData[p][sliceOffset + vox]  = params[p];
			residualData[sliceOffset + vox] = residual;
		};
		dispatch_apply(voxelsPerSlice, global_queue, processVoxel);
		
		if (verbose)
		{
			clock_t loopEnd = clock();
			if (voxCount > 0)
				std::cout << voxCount << " unmasked voxels, CPU time per voxel was "
				          << ((loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC)) << " s, ";
			std::cout << "finished." << std::endl;
		}
	}
    time_t procEnd = time(NULL);
    struct tm *localEnd = localtime(&procEnd);
	char theTime[MAXSTR];
    strftime(theTime, MAXSTR, "%H:%M:%S", localEnd);
	std::cout << "Finished processing at " << theTime << "Run-time was " 
	          << difftime(procEnd, procStart) << " s." << std::endl;

	write_results();
	// Clean up memory
	free(SPGR);
	for (int p = 0; p < nPhases; p++)
		free(SSFP[p]);
	free(SSFP);
	free(M0Data);
	free(B0Data);
	free(B1Data);
	free(maskData);
	exit(EXIT_SUCCESS);
}

