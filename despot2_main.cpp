/*
 *  despot1_main.cpp
 *
 *  Created by Tobias Wood on 23/01/2012.
 *  Copyright 2012 Tobias Wood. All rights reserved.
 *
 */

#include <time.h>
#include <getopt.h>
#include <iostream>
#include <atomic>
#include <Eigen/Dense>

#include "NiftiImage.h"
#include "DESPOT.h"
#include "DESPOT_Functors.h"
#include "RegionContraction.h"

#define USE_PROCPAR
#ifdef USE_PROCPAR
	#include "procpar.h"
	using namespace Recon;
#endif

using namespace std;
using namespace Eigen;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string credit {
"despot2 - Written by tobias.wood@kcl.ac.uk, based on work by Sean Deoni. \n\
Acknowledgements greatfully received, grant discussions welcome."
};

const string usage {
"Usage is: despot2 [options] T1_map ssfp_files\n\
\
Options:\n\
	--help, -h        : Print this message.\n\
	--mask, -m file   : Mask input with specified file.\n\
	--out, -o path    : Add a prefix to the output filenames.\n\
	--B0 file         : B0 Map file.\n\
	--B1 file         : B1 Map file.\n\
	--M0 file         : Proton density file.\n\
	--verbose, -v     : Print slice processing times.\n\
	--start_slice N   : Start processing from slice N.\n\
	--end_slice   N   : Finish processing at slice N.\n\
	--tesla, -t 3     : Enables DESPOT-FM with boundaries suitable for 3T\n\
	            7     : Boundaries suitable for 7T (default)\n\
	            u     : User specified boundaries from stdin.\n"

};

// tesla == 0 means NO DESPOT-FM
static int tesla = 0, fitB0 = false, verbose = false, start_slice = -1, end_slice = -1;
static string outPrefix;
static struct option long_options[] =
{
	{"B0", required_argument, 0, '0'},
	{"B1", required_argument, 0, '1'},
	{"M0", required_argument, 0, 'M'},
	{"help", no_argument, 0, 'h'},
	{"mask", required_argument, 0, 'm'},
	{"tesla", required_argument, 0, 't'},
	{"verbose", no_argument, 0, 'v'},
	{"start_slice", required_argument, 0, 'S'},
	{"end_slice", required_argument, 0, 'E'},
	{0, 0, 0, 0}
};

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	cout << credit << endl;
	if (argc < 4) {
		cout << usage << endl;
		return EXIT_FAILURE;
	}
	Eigen::initParallel();
	NiftiImage inFile, savedHeader;
	double *maskData = NULL, *B0Data = NULL, *B1Data = NULL, *T1Data = NULL,
	       *M0Data = NULL;
	string procPath;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hm:o:vt:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				inFile.open(optarg, 'r');
				maskData = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case '0':
				cout << "Reading B0 file: " << optarg << endl;
				inFile.open(optarg, 'r');
				B0Data = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case '1':
				cout << "Reading B1 file: " << optarg << endl;
				inFile.open(optarg, 'r');
				B1Data = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case 'M':
				cout << "Reading M0 file " << optarg;
				inFile.open(optarg, 'r');
				M0Data = inFile.readVolume<double>(0);
				inFile.close();
				break;
			case 't':
				switch (*optarg) {
					case '3': tesla = 3; break;
					case '7': tesla = 7; break;
					case 'u': tesla = -1; break;
					default:
						cout << "Unknown boundaries type " << optarg << endl;
						abort();
						break;
				}
				cout << "Using " << tesla << "T boundaries." << endl;
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
			case 'h':
				cout << usage << endl;				
				return EXIT_FAILURE;
		}
	}
	if ((tesla != 0) && !B0Data)
		fitB0 = true;
	cout << "Reading T1 Map from: " << argv[optind] << endl;
	inFile.open(argv[optind++], 'r');
	T1Data = inFile.readVolume<double>(0);
	inFile.close();
	savedHeader = inFile;
	//**************************************************************************
	// Gather SSFP Data
	//**************************************************************************
	int nFlip, nPhases;
	nPhases = argc - optind;
	vector<DESPOTConstants> consts(nPhases);
	VectorXd ssfpAngles;
	int voxelsPerSlice, voxelsPerVolume;
	double **ssfpData = (double **)malloc(nPhases * sizeof(double *));
	for (size_t p = 0; p < nPhases; p++) {
		cout << "Reading SSFP header from " << argv[optind] << endl;
		inFile.open(argv[optind], 'r');
		if (p == 0)
		{	// Read nFlip, TR and flip angles from first file
			nFlip = inFile.dim(4);
			voxelsPerSlice = inFile.voxelsPerSlice();
			voxelsPerVolume = inFile.voxelsPerVolume();
			ssfpAngles.resize(nFlip, 1);
			
			#ifdef USE_PROCPAR
			ParameterList pars;
			if (ReadProcpar(inFile.basename() + ".procpar", pars)) {
				consts[0].TR = RealValue(pars, "tr");
				for (int i = 0; i < nFlip; i++)
					ssfpAngles[i] = RealValue(pars, "flip1", i);
			} else
			#endif
			{
				cout << "Enter SSFP TR (s): " << flush;
				cin >> consts[0].TR;
				cout << "Enter " << nFlip << " flip angles (degrees): " << flush;
				for (int i = 0; i < ssfpAngles.size(); i++)
					cin >> ssfpAngles[i];
			}
		} else {
			consts[p].TR = consts[0].TR;
		}
		#ifdef USE_PROCPAR
		ParameterList pars;
		if (ReadProcpar(inFile.basename() + ".procpar", pars)) {
			consts[p].phase = RealValue(pars, "rfphase") * M_PI / 180.;
		} else
		#endif
		{
			cout << "Enter phase-cycling (degrees): " << flush;
			cin >> consts[p].phase; consts[p].phase *= M_PI / 180.;
		}
		cout << "Reading SSFP data..." << endl;
		ssfpData[p] = inFile.readAllVolumes<double>();
		// Don't close the first header because we've saved it to write the
		// results, and FSLIO gets fussy about cloning closed headers
		inFile.close();
		optind++;
	}
	ssfpAngles *= M_PI / 180.;
	
	if (optind != argc) {
		cerr << "Unprocessed arguments supplied.\n" << usage;
		exit(EXIT_FAILURE);
	}
	
	// Set up boundaries for DESPOT-FM if needed
	ArrayXd loBounds, hiBounds;
	const long nP = DESPOT2FM::inputs();
	if (tesla != 0) {
		if (tesla > 0) {
			loBounds = DESPOT2FM::defaultLo(tesla);
			hiBounds = DESPOT2FM::defaultHi(tesla);
		} else if (tesla < 0) {
			cout << "Enter " << nP << " parameter pairs (low then high): " << flush;
			for (int i = 0; i < nP; i++) cin >> loBounds[i] >> hiBounds[i];
		}
		// If fitting, give a suitable range and allocate results memory
		if (fitB0) {
			loBounds[0] = -0.5 / consts[0].TR;
			hiBounds[0] =  0.5 / consts[0].TR;
			B0Data = new double[voxelsPerVolume];
		} else { // Otherwise fix and let functors pick up the specified value
			loBounds[0] = 0.;
			hiBounds[0] = 0.;
		}
	}
	
	if (verbose) {
		cout << "SSFP Angles (deg): " << ssfpAngles.transpose() * 180 / M_PI << endl;
		if (tesla != 0)
			cout << "Low bounds: " << loBounds.transpose() << endl
				 << "Hi bounds:  " << hiBounds.transpose() << endl;
	}
	//**************************************************************************
	// Set up results data
	//**************************************************************************
	vector<double *> paramsData(nP);
	for (int p = 0; p < nP; p++)
		paramsData[p] = new double[voxelsPerVolume];
	double *residuals = new double[voxelsPerVolume];
	//**************************************************************************
	// Do the fitting
	//**************************************************************************
    time_t procStart = time(NULL);
	if ((start_slice < 0) || (start_slice >= inFile.dim(3)))
		start_slice = 0;
	if ((end_slice < 0) || (end_slice > inFile.dim(3)))
		end_slice = inFile.dim(3);
	for (int slice = start_slice; slice < end_slice; slice++) {
		// Read in data
		if (verbose)
			cout << "Starting slice " << slice << "..." << flush;
		
		atomic<int> voxCount{0};
		const int sliceOffset = slice * voxelsPerSlice;
		clock_t loopStart = clock();
		function<void (const int&)> processVox = [&] (const int &vox) {
			// Set up parameters and constants
			double residual = 0., T1 = 0.;
			ArrayXd params(nP);
			if (!maskData || ((maskData[sliceOffset + vox] > 0.) && (T1Data[sliceOffset + vox] > 0.)))
			{	// Zero T1 causes zero-pivot error.
				voxCount++;
				T1 = T1Data[sliceOffset + vox];
				// Gather signals.
				vector<VectorXd> signals;
				for (int p = 0; p < nPhases; p++) {
					if (B0Data) consts[p].B0 = B0Data[sliceOffset + vox];
					if (B1Data)	consts[p].B1 = B1Data[sliceOffset + vox];
					VectorXd temp(nFlip);
					for (int i = 0; i < nFlip; i++)
						temp(i) = ssfpData[p][i*voxelsPerVolume + sliceOffset + vox];
					signals.push_back(temp);
				}
				
				if (tesla == 0) {
					// Choose phase with accumulated phase closest to 180 and then classic DESPOT2
					int index = 0;
					double bestPhase = DBL_MAX;
					for (int p = 0; p < nPhases; p++) {
						double thisPhase = (consts[p].B0 * consts[p].TR * 2 * M_PI) + consts[p].phase;
						if (fabs(fmod(thisPhase - M_PI, 2 * M_PI)) < bestPhase) {
							bestPhase = fabs(fmod(thisPhase - M_PI, 2 * M_PI));
							index = p;
						}
					}
					residual = classicDESPOT2(ssfpAngles, signals[index], consts[index].TR, T1, consts[index].B1, params[1], params[2]);
					params[0] = consts[index].B0;
				} else {
					// DESPOT2-FM
					DESPOT2FM tc(ssfpAngles, signals, consts, T1, false, fitB0);
					ArrayXd params(nP);
					residual = regionContraction<DESPOT2FM>(params, tc, loBounds, hiBounds);				
				}
			}
			for (int p = 0; p < nP; p++) {
				paramsData[p][sliceOffset + vox] = params[p];
			}
			residuals[sliceOffset + vox] = residual;
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
	cout << "Finished processing at " << theTime << ". Run-time was " 
	          << difftime(procEnd, procStart) << " s." << endl;
	savedHeader.setDim(4, 1);
	savedHeader.setDatatype(NIFTI_TYPE_FLOAT32);
	for (int p = 0; p < nP; p++) {
		savedHeader.open(outPrefix + DESPOT2FM::names()[p], NiftiImage::NIFTI_WRITE);
		savedHeader.writeVolume(0, paramsData[p]);
		savedHeader.close();
	}
	savedHeader.open(outPrefix + "D2_Residual.nii.gz", NiftiImage::NIFTI_WRITE);
	savedHeader.writeVolume(0, residuals);
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
	return EXIT_SUCCESS;
}
