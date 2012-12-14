/*
 *  despot1_main.c
 *  MacRI
 *
 *  Created by Tobias Wood on 23/01/2012.
 *  Copyright 2012 Tobias Wood. All rights reserved.
 *
 */

#include <time.h>
#include <getopt.h>
#include <signal.h>

#include <atomic>

#include <unsupported/Eigen/NonLinearOptimization>
#include "DESPOT.h"
#include "DESPOT_Functors.h"

#define USE_PROCPAR

#ifdef USE_PROCPAR
	#include "procpar.h"
	using namespace Recon;
#endif

const std::string usage("Usage is: despot2 [options] output_prefix T1_map ssfp_files\n\
\
Options:\n\
	--mask, -m file   : Mask input with specified file.\n\
	--B0 file         : B0 Map file.\n\
	--B1 file         : B1 Map file.\n\
	--M0 file         : Proton density file.\n\
	--lm              : Use Levenberg-Marquardt instead of Region Contraction.\n\
	--verbose, -v     : Print slice processing times.\n\
	--start_slice N   : Start processing from slice N.\n\
	--end_slice   N   : Finish processing at slice N.\n");

static int levMar = false, verbose = false, start_slice = -1, end_slice = -1;
static std::string outPrefix;
static struct option long_options[] =
{
	{"B0", required_argument, 0, '0'},
	{"B1", required_argument, 0, '1'},
	{"M0", required_argument, 0, 'M'},
	{"mask", required_argument, 0, 'm'},
	{"lm", no_argument, &levMar, true},
	{"verbose", no_argument, 0, 'v'},
	{"start_slice", required_argument, 0, 'S'},
	{"end_slice", required_argument, 0, 'E'},
	{0, 0, 0, 0}
};

//******************************************************************************
// SIGTERM interrupt handler - for ensuring data gets saved even on a ctrl-c
//******************************************************************************
NiftiImage savedHeader;
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
	int nSSFP, nPhases;
	double ssfpTR;
	NiftiImage inFile;
	double *maskData = NULL, *B0Data = NULL, *B1Data = NULL, *T1Data = NULL,
	       *M0Data = NULL;
	std::string procPath;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "m:vz", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 'm':
				fprintf(stdout, "Reading mask file %s.\n", optarg);
				inFile.open(optarg, 'r');
				maskData = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case '0':
				fprintf(stdout, "Reading B0 file %s.\n", optarg);
				inFile.open(optarg, 'r');
				B0Data = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case '1':
				fprintf(stdout, "Reading B1 file %s.\n", optarg);
				inFile.open(optarg, 'r');
				B1Data = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case 'M':
				fprintf(stdout, "Reading M0 file %s.\n", optarg);
				inFile.open(optarg, 'r');
				M0Data = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case 'v':
				verbose = true;
				break;
			case 'S':
				start_slice = atoi(optarg);
				break;
			case 'E':
				end_slice = atoi(optarg);
				break;
			case 0:
				// Just a flag
				break;
			case '?': // getopt will print an error message
				abort();
		}
	}
	
	std::cout << "Output prefix will be: " << argv[optind] << std::endl;
	outPrefix = argv[optind++];
	std::cout << "Reading T1 Map from: " << argv[optind] << std::endl;
	inFile.open(argv[optind++], 'r');
	T1Data = inFile.readVolume<double>(0);
	inFile.close();
	savedHeader = inFile;
	//**************************************************************************
	// Gather SSFP Data
	//**************************************************************************
	nPhases = argc - optind;
	VectorXd ssfpPhases(nPhases), ssfpAngles;
	int voxelsPerSlice, totalVoxels;
	double **ssfpData = (double **)malloc(nPhases * sizeof(double *));
	for (size_t p = 0; p < nPhases; p++)
	{
		std::cout << "Reading SSFP header from " << argv[optind] << std::endl;
		inFile.open(argv[optind], 'r');
		if (p == 0)
		{	// Read nSSFP, TR and flip angles from first file
			nSSFP = inFile.nt();
			voxelsPerSlice = inFile.nx() * inFile.ny();
			totalVoxels = voxelsPerSlice * inFile.nz();
			ssfpAngles.resize(nSSFP, 1);
			
			#ifdef USE_PROCPAR
			ParameterList pars;
			if (ReadProcpar(inFile.basename() + ".procpar", pars)) {
				ssfpTR = RealValue(pars, "tr");
				for (int i = 0; i < nSSFP; i++)
					ssfpAngles[i] = RealValue(pars, "flip1", i);
			} else
			#endif
			{
				std::cout << "Enter SSFP TR (s): " << std::flush;
				std::cin >> ssfpTR;
				std::cout << "Enter " << nSSFP << " flip angles (degrees): " << std::flush;
				for (int i = 0; i < ssfpAngles.size(); i++)
					std::cin >> ssfpAngles[i];
			}
		}
		#ifdef USE_PROCPAR
		ParameterList pars;
		if (ReadProcpar(inFile.basename() + ".procpar", pars)) {
			ssfpPhases[p] = RealValue(pars, "rfphase");
		} else
		#endif
		{
			std::cout << "Enter phase-cycling (degrees): " << std::flush;
			std::cin >> ssfpPhases[p];
		}
		std::cout << "Reading SSFP data..." << std::endl;
		ssfpData[p] = inFile.readAllVolumes<double>();
		// Don't close the first header because we've saved it to write the
		// results, and FSLIO gets fussy about cloning closed headers
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
	
	if (verbose)
	{
		std::cout << "SSFP Angles (deg): " << ssfpAngles.transpose() * 180 / M_PI << std::endl
		          << "SSFP Phases (deg): " << ssfpPhases.transpose() * 180 / M_PI << std::endl;
	}
	//**************************************************************************
	// Set up results data
	//**************************************************************************
	double *residualData = new double[totalVoxels];
	double *T2Data = new double[totalVoxels];
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
    time_t procStart = time(NULL);
	if ((start_slice < 0) || (start_slice >= inFile.nz()))
		start_slice = 0;
	if ((end_slice < 0) || (end_slice > inFile.nz()))
		end_slice = inFile.nz();
	for (int slice = start_slice; slice < end_slice; slice++) {
		// Read in data
		if (verbose)
			cout << "Starting slice " << slice << "..." << flush;
		
		atomic<int> voxCount{0};
		const int sliceOffset = slice * voxelsPerSlice;
		clock_t loopStart = clock();
		function<void (const int&)> processVox = [&] (const int &vox) {
			// Set up parameters and constants
			double M0 = 0., T1 = 0., T2 = 0., B0 = 0, B1 = 1., residual = 0.;
			if (!maskData || ((maskData[sliceOffset + vox] > 0.) && (T1Data[sliceOffset + vox] > 0.)))
			{	// Zero T1 causes zero-pivot error.
				voxCount++;
				T1 = T1Data[sliceOffset + vox];
				if (B0Data) B0 = B0Data[sliceOffset + vox];
				if (B1Data)	B1 = B1Data[sliceOffset + vox];
				// Gather signals.
				std::vector<VectorXd> signals;
				for (int p = 0; p < nPhases; p++)
				{
					VectorXd temp(nSSFP);
					for (int i = 0; i < nSSFP; i++)
						temp(i) = ssfpData[p][i*totalVoxels + sliceOffset + vox];
					signals.push_back(temp);
				}
				// Choose phase with accumulated phase closest to 180
				int index = 0;
				double bestPhase = DBL_MAX;
				for (int p = 0; p < nPhases; p++)
				{
					double thisPhase = (B0 * ssfpTR * 2 * M_PI) + ssfpPhases[p];
					if (fabs(fmod(thisPhase - M_PI, 2 * M_PI)) < bestPhase)
					{
						bestPhase = fabs(fmod(thisPhase - M_PI, 2 * M_PI));
						index = p;
					}
				}
				residual = classicDESPOT2(ssfpAngles, signals[index], ssfpTR, T1, B1, &M0, &T2);
			}
			T2Data[sliceOffset + vox] = clamp(T2, 0., 0.5);
			residualData[sliceOffset + vox] = residual;
		};
		apply_for(voxelsPerSlice, processVox);
		
		if (verbose) {
			clock_t loopEnd = clock();
			if (voxCount > 0)
				cout << voxCount << " unmasked voxels, CPU time per voxel was "
				          << ((loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC)) << " s, ";
			cout << "finished." << endl;
		}
	}
    time_t procEnd = time(NULL);
    struct tm *localEnd = localtime(&procEnd);
	char theTime[512];
    strftime(theTime, 512, "%H:%M:%S", localEnd);
	fprintf(stdout, "Finished processing at %s. Run-time was %f s.\n", theTime, difftime(procEnd, procStart));
	savedHeader.setnt(1);
	savedHeader.setDatatype(NIFTI_TYPE_FLOAT32);
	savedHeader.open(outPrefix + "_T2.nii.gz", NiftiImage::NIFTI_WRITE);
	savedHeader.writeVolume(0, T2Data);
	savedHeader.close();
	// Clean up memory
	for (int p = 0; p < nPhases; p++)
		free(ssfpData[p]);
	if (B0Data)
		free(B0Data);
	if (B1Data)
		free(B1Data);
	if (maskData)
		free(maskData);
	exit(EXIT_SUCCESS);
}
