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

#include "procpar.h"
#include "DESPOT_Functors.h"
#include "RegionContraction.h"

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const std::string usage {
"Usage is: mcdespot [options] output_prefix spgr_file ssfp_file1 (ssfp_fileN)\n\
\
Options:\n\
	-v, --verbose     : Print extra information.\n\
	-m, --mask file   : Mask input with specified file.\n\
	--1, --2, --3     : Use 1, 2 or 3 component model (default 2).\n\
	--M0 file         : M0 Map file.\n\
	--B0 file         : B0 Map file.\n\
	--B1 file         : B1 Map file.\n\
	--start_slice n   : Only start processing at slice n.\n\
	--end_slice n     : Finish at slice n-1.\n\
	-d, --drop        : Drop certain flip-angles (Read from standard in).\n\
	-n, --normalise   : Normalise signals to maximum (Ignore M0).\n\
	-s, --samples n   : Use n samples for region contraction (Default 5000).\n\
	-r, --retain  n   : Retain n samples for new boundary (Default 50).\n\
	-c, --contract n  : Contract a maximum of n times (Default 10).\n\
	-e, --expand n    : Re-expand boundary by percentage n (Default 0).\n\
	-b 3              : Boundaries suitable for 3T\n\
	   7              : Boundaries suitable for 7T (default)\n\
	   u              : User specified boundaries from stdin.\n"
};

static int verbose = false, normalise = false, drop = false, start_slice = -1, end_slice = -1,
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
	{"normalise", no_argument, &normalise, true},
	{"samples", required_argument, 0, 's'},
	{"retain", required_argument, 0, 'r'},
	{"contract", required_argument, 0, 'c'},
	{"expand", required_argument, 0, 'e'},
	{"drop", no_argument, &drop, true},
	{"1", no_argument, &components, 1},
	{"2", no_argument, &components, 2},
	{"3", no_argument, &components, 3},
	{0, 0, 0, 0}
};
//******************************************************************************
// SIGTERM interrupt handler - for ensuring data gets saved even on a ctrl-c
//******************************************************************************
NiftiImage savedHeader;
double **paramsData;
double *residualData;

void int_handler(int sig);
void int_handler(int sig)
{
	fprintf(stdout, "Processing terminated. Writing currently processed data.\n");
	switch (components)
	{
		case 1: write_results<OneComponent>(outPrefix, paramsData, residualData, savedHeader); break;
		case 2: write_results<TwoComponent>(outPrefix, paramsData, residualData, savedHeader); break;
		case 3: write_results<ThreeComponent>(outPrefix, paramsData, residualData, savedHeader); break;
	}
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
		std::cerr << usage << std::endl;
		exit(EXIT_FAILURE);
	}
	Eigen::initParallel();
	__block double spgrTR, ssfpTR,
	               *maskData = NULL, *M0Data = NULL,
				   *B0Data = NULL, *B1Data = NULL;
	NiftiImage inFile;
	int nSPGR, nPhases, nSSFP;
	std::string procPath;
	par_t *pars;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "b:m:vznds:r:c:e:", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 'm':
				std::cout << "Reading mask file " << optarg << std::endl;
				inFile.open(optarg, NIFTI_READ);
				maskData = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case '0':
				std::cout << "Reading B0 file " << optarg << std::endl;
				inFile.open(optarg, NIFTI_READ);
				B0Data = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case '1':
				std::cout << "Reading B1 file " << optarg << std::endl;
				inFile.open(optarg, NIFTI_READ);
				B1Data = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case 'M':
				std::cout << "Reading M0 file " << optarg << std::endl;
				inFile.open(optarg, NIFTI_READ);
				M0Data = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case 'v': verbose = true; break;
			case 'S': start_slice = atoi(optarg); break;
			case 'E': end_slice = atoi(optarg); break;
			case 'n': normalise = true; break;
			case 'd': drop = true; break;
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
						abort();
						break;
				}
				break;
			case 0:
				// Just a flag
				break;
			case '?': // getopt will print an error message
				abort();
		}
	}
	outPrefix = argv[optind++];
	if (verbose)
		std::cout << "Output prefix will be: " << outPrefix << std::endl;
	//**************************************************************************
	// Gather SPGR Data
	//**************************************************************************
	std::cout << "Opening SPGR file: " << argv[optind] << std::endl;
	inFile.open(argv[optind], NIFTI_READ);
	nSPGR = inFile.nt();
	VectorXd spgrAngles(nSPGR);
	procPath = inFile.basename() + ".procpar";
	pars = readProcpar(procPath.c_str());
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
	int voxelsPerSlice = inFile.nx() * inFile.ny();
	int totalVoxels = inFile.voxelsPerVolume();
	std::cout << "Reading SPGR data..." << std::endl;
	__block double *SPGR = inFile.readAllVolumes<double>();
	// Save this header to output the results files.
	inFile.close();
	savedHeader = inFile;
	optind++;
	//**************************************************************************
	// Gather SSFP Data
	//**************************************************************************
	nPhases = argc - optind;
	VectorXd ssfpPhases(nPhases), ssfpAngles;
	__block double **SSFP = (double **)malloc(nPhases * sizeof(double *));
	for (size_t p = 0; p < nPhases; p++)
	{
		std::cout << "Reading SSFP header from " << argv[optind] << std::endl;
		inFile.open(argv[optind], NIFTI_READ);
		procPath = inFile.basename() + ".procpar";
		pars = readProcpar(procPath.c_str());
		if (p == 0)
		{	// Read nSSFP, TR and flip angles from first file
			nSSFP = inFile.nt();
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
		// Check that all files have consistent dims
		eigen_assert((inFile.nx() == savedHeader.nx()) &&
		             (inFile.ny() == savedHeader.ny()) &&
					 (inFile.nx() == savedHeader.nz()));
		eigen_assert((inFile.nt() == nSSFP));
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
		SSFP[p] = inFile.readAllVolumes<double>();
		inFile.close();
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
	std::cout << "Using " << components << " component model." << std::endl;
	switch (components)
	{
		case 1: nP = OneComponent::nP; break;
		case 2: nP = TwoComponent::nP; break;
		case 3: nP = ThreeComponent::nP; break;
	}
		
	__block VectorXd loBounds(nP), hiBounds(nP);
	__block VectorXi loConstraints(nP), hiConstraints(nP);
	
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
			switch (components)
			{
				case 1: loBounds[i] = OneComponent::lo3Bounds[i]; hiBounds[i] = OneComponent::hi3Bounds[i]; break;
				case 2: loBounds[i] = TwoComponent::lo3Bounds[i]; hiBounds[i] = TwoComponent::hi3Bounds[i]; break;
				case 3: loBounds[i] = ThreeComponent::lo3Bounds[i]; hiBounds[i] = ThreeComponent::hi3Bounds[i]; break;
			}
		} else if (tesla == 7) {
			switch (components)
			{
				case 1: loBounds[i] = OneComponent::lo7Bounds[i]; hiBounds[i] = OneComponent::hi7Bounds[i]; break;
				case 2: loBounds[i] = TwoComponent::lo7Bounds[i]; hiBounds[i] = TwoComponent::hi7Bounds[i]; break;
				case 3: loBounds[i] = ThreeComponent::lo7Bounds[i]; hiBounds[i] = ThreeComponent::hi7Bounds[i]; break;
			}
		} else {
			std::cin >> loBounds[i] >> hiBounds[i];
		}
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
	// Select which angles to use in the analysis
	//**************************************************************************	
	VectorXi spgrKeep(nSPGR), ssfpKeep(nSSFP);
	if (drop)
	{
		std::cout << "Choose SPGR angles to use (1 to keep, 0 to drop, " << nSPGR << " values): ";
		for (int i = 0; i < nSPGR; i++) std::cin >> spgrKeep[i];
		std::cout << "Using " << spgrKeep.sum() << " SPGR angles. " << spgrKeep.transpose() << std::endl;
		std::cout << "Choose SSFP angles to use (1 to keep, 0 to drop, " << nSSFP << " values): ";
		for (int i = 0; i < nSSFP; i++) std::cin >> ssfpKeep[i];
		std::cout << "Using " << ssfpKeep.sum() << " SSFP angles. " << ssfpKeep.transpose() << std::endl;
		VectorXd temp = spgrAngles;
		spgrAngles.resize(spgrKeep.sum());
		int angle = 0;
		for (int i = 0; i < nSPGR; i++)
			if (spgrKeep(i)) spgrAngles(angle++) = temp(i);
		temp = ssfpAngles;
		ssfpAngles.resize(ssfpKeep.sum());
		angle = 0;
		for (int i = 0; i < nSSFP; i++)
			if (ssfpKeep(i)) ssfpAngles(angle++) = temp(i);
	}
	else
	{
		spgrKeep.setOnes();
		ssfpKeep.setOnes();
	}
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
    time_t procStart = time(NULL);
	if ((start_slice < 0) || (start_slice >= savedHeader.nz()))
		start_slice = 0;
	if ((end_slice < 0) || (end_slice > savedHeader.nz()))
		end_slice = savedHeader.nz();
	signal(SIGINT, int_handler);	// If we've got here there's actually allocated data to save
	std::cout << "Starting processing." << std::endl;
	dispatch_queue_t global_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
	for (int slice = start_slice; slice < end_slice; slice++)
	{
		if (verbose)
			std::cout << "Starting slice " << slice << "..." << std::flush;
		__block int voxCount = 0;
		__block const int sliceOffset = slice * voxelsPerSlice;
		clock_t loopStart = clock();
		//for (int vox = 0; vox < voxelsPerSlice; vox++)
		void (^processVoxel)(size_t vox) = ^(size_t vox)
		{
			double B0 = 0., B1 = 1., residual = 0.;
			VectorXd params(nP); params.setZero();
			if ((!maskData || (maskData[sliceOffset + vox] > 0.)) &&
			    (!M0Data || (M0Data[sliceOffset + vox] > 0.)))
			{
				AtomicAdd(1, &voxCount);
				// Must take copies before altering inside a block
				VectorXd localLo = loBounds, localHi = hiBounds;
				if (M0Data)
				{
					localLo(0) = (double)M0Data[sliceOffset + vox];
					localHi(0) = (double)M0Data[sliceOffset + vox];
				}
				
				if (B0Data) B0 = (double)B0Data[sliceOffset + vox];
				if (B1Data) B1 = (double)B1Data[sliceOffset + vox];
				
				VectorXd SPGR_signal(spgrKeep.sum());
				int vol = 0;
				for (int i = 0; i < nSPGR; i++)
				{
					if (spgrKeep(i))
						SPGR_signal[vol++] = SPGR[totalVoxels*i + sliceOffset + vox];
				}
				if (normalise)
					SPGR_signal /= SPGR_signal.mean();
				std::vector<VectorXd> SSFP_signals;
				for (int p = 0; p < nPhases; p++)
				{
					VectorXd temp(ssfpKeep.sum());
					vol = 0;
					for (int i = 0; i < nSSFP; i++)
					{
						if (ssfpKeep(i))
							temp[vol++] = SSFP[p][totalVoxels*i + sliceOffset + vox];
					}
					if (normalise)
						temp /= temp.mean();
					SSFP_signals.push_back(temp);
				}
				
				int rSeed = static_cast<int>(time(NULL));
				switch (components)
				{
					case 1: {
						OneComponent tc(spgrAngles, SPGR_signal, ssfpAngles, ssfpPhases, SSFP_signals, spgrTR, ssfpTR, B0, B1, normalise);
						residual = regionContraction<OneComponent>(params, tc, localLo, localHi,
															       loConstraints, hiConstraints,
															       samples, retain, contract, 0.05, expand, rSeed);
						break; }
					case 2: {
						TwoComponent tc(spgrAngles, SPGR_signal, ssfpAngles, ssfpPhases, SSFP_signals, spgrTR, ssfpTR, B0, B1, normalise);
						residual = regionContraction<TwoComponent>(params, tc, localLo, localHi,
															       loConstraints, hiConstraints,
															       samples, retain, contract, 0.05, expand, rSeed);
						break; }
					case 3: {
						ThreeComponent tc(spgrAngles, SPGR_signal, ssfpAngles, ssfpPhases, SSFP_signals, spgrTR, ssfpTR, B0, B1, normalise);
						residual = regionContraction<ThreeComponent>(params, tc, localLo, localHi,
																 	loConstraints, hiConstraints,
																 	samples, retain, contract, 0.05, expand, rSeed);
						break; }
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
	switch (components)
	{
		case 1: write_results<OneComponent>(outPrefix, paramsData, residualData, savedHeader); break;
		case 2: write_results<TwoComponent>(outPrefix, paramsData, residualData, savedHeader); break;
		case 3: write_results<ThreeComponent>(outPrefix, paramsData, residualData, savedHeader); break;
	}
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

