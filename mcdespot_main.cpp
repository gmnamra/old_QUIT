/*
 *  mcdespot_main.c
 *  Fitting
 *
 *  Created by Tobias Wood on 14/02/2012.
 *  Copyright 2012 Tobias Wood. All rights reserved.
 *
 */

#include <time.h>
#include <getopt.h>
#include <signal.h>

#include <iostream>
#include <fstream>

#include <functional>
#include <thread>
#include <atomic>

using namespace std;

#include "procpar.h"
#include "DESPOT_Functors.h"
#include "RegionContraction.h"

void apply_for(const int max, const function<void(int)> f) {
	vector<thread> pool;
	
	atomic<int> complete{0};
	
	function<void()> worker = [&complete, &max, &f]() {
		cout << "Worker " << this_thread::get_id() << " launched" << endl;
		while (complete < max) {
			int local = atomic_fetch_add(&complete, 1);
			f(local);
		}
		cout << "Worker " << this_thread::get_id() << " finished" << endl;
	};
	
	cout << "Hardware concurrency: " << thread::hardware_concurrency() << endl;
	for (int t = 0; t < thread::hardware_concurrency(); t++)
		pool.push_back(thread(worker));
	for (int t = 0; t < pool.size(); t++)
		pool[t].join();

}


//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: mcdespot [options] input_file output_prefix\n\
\n\
The input file must consist of at least one line such as:\n\
SPGR path_to_spgr_file (path_to_B1_file) \n\
SSFP path_to_ssfp_file phase_cycling (path_to_B0) (path_to_B1) \n\
\n\
Phase-cycling is specified in degrees. If no B0/B1 correction is specified,\n\
default values of B0 = 0 Hz and B1 = 1 will be used.\n\
\n\
Options:\n\
	-v, --verbose     : Print extra information.\n\
	-m, --mask file   : Mask input with specified file.\n\
	--1, --2, --3     : Use 1, 2 or 3 component model (default 2).\n\
	--M0 file         : M0 Map file.\n\
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
static string outPrefix;
static struct option long_options[] =
{
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
		cerr << usage << endl;
		exit(EXIT_FAILURE);
	}
	Eigen::initParallel();
	double spgrTR, ssfpTR,
	               *maskData = NULL, *M0Data = NULL;
	vector<NiftiImage *> SPGR_files, SSFP_files, SPGR_B1_files, SSFP_B0_files, SSFP_B1_files;
	NiftiImage inHeader;
	int nSPGR = 0, nSSFP = 0;
	string procPath;
	par_t *pars;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "b:m:vznds:r:c:e:", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				inHeader.open(optarg, NiftiImage::NIFTI_READ);
				maskData = inHeader.readVolume<double>(0);
				inHeader.close();
				break;
			case 'M':
				cout << "Reading M0 file " << optarg << endl;
				inHeader.open(optarg, NiftiImage::NIFTI_READ);
				M0Data = inHeader.readVolume<double>(0);
				inHeader.close();
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
						cout << "Using 3T boundaries.\n";
						tesla = 3;
						break;
					case '7':
						cout << "Using 7T boundaries.\n";
						tesla = 7;
						break;
					case 'u':
						tesla = -1;
						break;
					default:
						cout << "Unknown boundaries type " << optarg << endl;
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
	if ((argc - optind) != 2)
	{
		cerr << "Incorrect number of arguments.";
		exit(EXIT_FAILURE);
	}

	//**************************************************************************
	// Read input file and set up corresponding SPGR & SSFP lists
	//**************************************************************************
	VectorXd spgrAngles, ssfpAngles;
	vector<double> ssfpPhases;
	cout << "Opening input file: " << argv[optind] << endl;
	ifstream inFile(argv[optind++]);
	if (!inFile)
	{
		cerr << "Could not open input file." << endl;
		exit(EXIT_FAILURE);
	}
	string nextLine;
	
	while (getline(inFile, nextLine))
	{
		stringstream thisLine(nextLine);
		string type, path;
		NiftiImage *inHdr;
		
		if (!(thisLine >> type >> path)) {
			cerr << "Could not read image type and path from input file line: " << nextLine << endl;
			exit(EXIT_FAILURE);
		}
		
		cout << "Opening image header: " << path << endl;
		inHdr = new NiftiImage(path, NiftiImage::NIFTI_READ);
		if ((SPGR_files.size() == 0) && (SSFP_files.size() == 0))
			savedHeader = *inHdr;
		
		if (!inHdr->volumesCompatible(savedHeader)) {
			cerr << "Physical space is incompatible with first image." << endl;
			exit(EXIT_FAILURE);
		}
		
		if (type == "SPGR") {
			if (SPGR_files.size() == 0) {
				nSPGR = inHdr->nt();
				spgrAngles.resize(nSPGR, 1);
				pars = readProcpar((inHdr->basename() + ".procpar").c_str());
				if (pars) {
					spgrTR = realVal(pars, "tr", 0);
					for (int i = 0; i < nSPGR; i++) spgrAngles[i] = realVal(pars, "flip1", i);
					freeProcpar(pars);
				} else {
					thisLine >> spgrTR;
					for (int i = 0; i < nSPGR; i++) thisLine >> spgrAngles[i];
				}
				spgrAngles *= M_PI / 180.;
			}
			SPGR_files.push_back(inHdr);
			if (thisLine >> path) // Read a path to B1 file
				inHdr = new NiftiImage(path, NiftiImage::NIFTI_READ);
			else
				inHdr = NULL;
			SPGR_B1_files.push_back(inHdr);
		} else if (type == "SSFP") {
			pars = readProcpar((inHdr->basename() + ".procpar").c_str());
			double phase;
			if (pars)
				phase = realVal(pars, "rfphase", 0);
			else
				inFile >> phase;
			ssfpPhases.push_back(phase * M_PI / 180.);
			if (SSFP_files.size() == 0) {
				nSSFP = inHdr->nt();
				ssfpAngles.resize(nSSFP, 1);
				if (pars) {
					ssfpTR = realVal(pars, "tr", 0);
					for (int i = 0; i < nSSFP; i++)
						ssfpAngles[i] = realVal(pars, "flip1", i);
				} else {
					inFile >> ssfpTR;
					for (int i = 0; i < ssfpAngles.size(); i++) inFile >> ssfpAngles[i];
				}
				ssfpAngles *= M_PI / 180.;
			}
			freeProcpar(pars);
			SSFP_files.push_back(inHdr);
			if (thisLine >> path)	// Read a path to B0 file
				inHdr = new NiftiImage(path, NiftiImage::NIFTI_READ);
			else
				inHdr = NULL;
			SSFP_B0_files.push_back(inHdr);
			if (thisLine >> path) // Read a path to B1 file
				inHdr = new NiftiImage(path, NiftiImage::NIFTI_READ);
			else
				inHdr = NULL;
			SSFP_B1_files.push_back(inHdr);
		} else {
			cerr << "Unknown scan type: " << type << ", must be SPGR or SSFP." << endl;
			exit(EXIT_FAILURE);
		}
	}
	if (SPGR_files.size() == 0 && SSFP_files.size() == 0)
	{
		cerr << "No input images specified." << endl;
		exit(EXIT_FAILURE);
	}
	//**************************************************************************
	// Select which angles to use in the analysis
	//**************************************************************************	
	VectorXi spgrKeep(nSPGR), ssfpKeep(nSSFP);
	spgrKeep.setOnes(); ssfpKeep.setOnes();
	if (drop) {
		cout << "Choose SPGR angles to use (1 to keep, 0 to drop, " << nSPGR << " values): ";
		for (int i = 0; i < nSPGR; i++) cin >> spgrKeep[i];
		cout << "Choose SSFP angles to use (1 to keep, 0 to drop, " << nSSFP << " values): ";
		for (int i = 0; i < nSSFP; i++) cin >> ssfpKeep[i];
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
	outPrefix = argv[optind++];
	if (verbose)
		cout << "Output prefix will be: " << outPrefix << endl;
	//**************************************************************************
	// Allocate results memory and set up boundaries
	//**************************************************************************
	int voxelsPerSlice = savedHeader.voxelsPerSlice();
	int voxelsPerVolume = savedHeader.voxelsPerVolume();	
	cout << "Using " << components << " component model." << endl;
	switch (components)
	{
		case 1: nP = OneComponent::nP; break;
		case 2: nP = TwoComponent::nP; break;
		case 3: nP = ThreeComponent::nP; break;
	}
		
	VectorXd loBounds(nP), hiBounds(nP);
	VectorXi loConstraints(nP), hiConstraints(nP);
	
	residualData = (double *)malloc(voxelsPerVolume * sizeof(double));
	paramsData = (double **)malloc(nP * sizeof(double *));
	if (tesla < 0)
		cout << "Enter " << nP << " parameter pairs (low then high): " << flush;
	for (int i = 0; i < nP; i++)
	{
		paramsData[i] = (double *)malloc(voxelsPerVolume * sizeof(double));
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
			cin >> loBounds[i] >> hiBounds[i];
		}
	}
	if (normalise)
	{
		loBounds[0] = 1.;
		hiBounds[0] = 1.;
	}
	if (verbose)
	{
		cout << "SPGR TR (s): " << spgrTR << " SPGR flip angles (deg): " << spgrAngles.transpose() * 180. / M_PI << endl;		
		cout << "SSFP TR (s): " << ssfpTR << " SSFP flip angles (deg): " << ssfpAngles.transpose() * 180. / M_PI << endl;	
		cout << "Low bounds: " << loBounds.transpose() << endl;
		cout << "Hi bounds:  " << hiBounds.transpose() << endl;
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
	cout << "Starting processing." << endl;
	for (int slice = start_slice; slice < end_slice; slice++)
	{
		if (verbose)
			cout << "Starting slice " << slice << "..." << flush;
		atomic<int> voxCount{0};
		const int sliceOffset = slice * voxelsPerSlice;
		
		// Read data for slice
		vector<double *> SPGR_vols, SPGR_B1s, SSFP_vols, SSFP_B0s, SSFP_B1s;
		for (int i = 0; i < SPGR_files.size(); i++)
		{
			SPGR_vols.push_back(SPGR_files[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1));
			if (SPGR_B1_files[i])
				SPGR_B1s.push_back(SPGR_B1_files[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1));
			else
				SPGR_B1s.push_back(NULL);
		}
		for (int i = 0; i < SSFP_files.size(); i++)
		{
			SSFP_vols.push_back(SSFP_files[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1));
			if (SSFP_B0_files[i])
				SSFP_B0s.push_back(SSFP_B0_files[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1));
			else
				SSFP_B0s.push_back(NULL);
			if (SSFP_B1_files[i])
				SSFP_B1s.push_back(SSFP_B1_files[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1));
			else
				SSFP_B1s.push_back(NULL);
		}
		
		clock_t loopStart = clock();
		//for (int vox = 0; vox < voxelsPerSlice; vox++)
		//void (^processVoxel)(size_t vox) = ^(size_t vox)
		function<void (const size_t&)> processVox = [&] (const size_t &vox)
		{
			double residual = 0.;
			VectorXd params(nP); params.setZero();
			if ((!maskData || (maskData[sliceOffset + vox] > 0.)) &&
			    (!M0Data || (M0Data[sliceOffset + vox] > 0.)))
			{
				//cout << "Voxel " << (vox % savedHeader.nx()) << " " << (vox / savedHeader.nx()) << endl;
				voxCount++;
				// Must take copies before altering inside a block
				VectorXd localLo = loBounds, localHi = hiBounds;
				if (M0Data)
				{
					localLo(0) = (double)M0Data[sliceOffset + vox];
					localHi(0) = (double)M0Data[sliceOffset + vox];
				}
								
				vector<VectorXd> SPGR_signals, SSFP_signals;
				VectorXd SPGR_B1(SPGR_vols.size()),
				         SSFP_B0(SSFP_vols.size()), SSFP_B1(SSFP_vols.size());
				
				for (int i = 0; i < SPGR_vols.size(); i++) {
					VectorXd temp(spgrKeep.sum());
					int vol = 0;
					for (int j = 0; j < nSPGR; j++) {
						if (spgrKeep(j))
							temp[vol++] = SPGR_vols[i][voxelsPerSlice*j + vox];
					}
					if (normalise)
						temp /= temp.mean();
					SPGR_signals.push_back(temp);
					SPGR_B1[i] = SPGR_B1s[i] ? SPGR_B1s[i][vox] : 1.;
				}
				
				for (int i = 0; i < SSFP_vols.size(); i++)
				{
					VectorXd temp(ssfpKeep.sum());
					int vol = 0;
					for (int j = 0; j < nSSFP; j++)
					{
						if (ssfpKeep(j))
							temp[vol++] = SSFP_vols[i][voxelsPerSlice*j + vox];
					}
					if (normalise)
						temp /= temp.mean();
					SSFP_signals.push_back(temp);
					SSFP_B0[i] = SSFP_B0s[i] ? SSFP_B0s[i][vox] : 0.;
					SSFP_B1[i] = SSFP_B1s[i] ? SSFP_B1s[i][vox] : 1.;
				}
				
				int rSeed = static_cast<int>(time(NULL));
				switch (components)
				{
					case 1: {
						OneComponent tc(spgrAngles, SPGR_signals, SPGR_B1, spgrTR,
						                ssfpAngles, ssfpPhases, SSFP_signals, SSFP_B0, SSFP_B1, ssfpTR,
										normalise);
						residual = regionContraction<OneComponent>(params, tc, localLo, localHi,
															       loConstraints, hiConstraints,
															       samples, retain, contract, 0.05, expand, rSeed);
						break; }
					case 2: {
						TwoComponent tc(spgrAngles, SPGR_signals, SPGR_B1, spgrTR,
						                ssfpAngles, ssfpPhases, SSFP_signals, SSFP_B0, SSFP_B1, ssfpTR,
										normalise);
						residual = regionContraction<TwoComponent>(params, tc, localLo, localHi,
															       loConstraints, hiConstraints,
															       samples, retain, contract, 0.05, expand, rSeed);
						break; }
					case 3: {
						ThreeComponent tc(spgrAngles, SPGR_signals, SPGR_B1, spgrTR,
						                ssfpAngles, ssfpPhases, SSFP_signals, SSFP_B0, SSFP_B1, ssfpTR,
										normalise);
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
		
		apply_for(voxelsPerSlice, processVox);
		
		
		// Clean up memory
		for (int i = 0; i < SPGR_vols.size(); i++)
		{
			free(SPGR_vols[i]);
			if (SPGR_B1s[i]) free(SPGR_B1s[i]);
		}
		for (int i = 0; i < SSFP_vols.size(); i++)
		{
			free(SSFP_vols[i]);
			if (SSFP_B0s[i]) free(SSFP_B0s[i]);
			if (SSFP_B1s[i]) free(SSFP_B1s[i]);
		}
		
		if (verbose)
		{
			clock_t loopEnd = clock();
			if (voxCount > 0)
				cout << voxCount << " unmasked voxels, CPU time per voxel was "
				          << ((loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC)) << " s, ";
			cout << "finished." << endl;
		}
	}
    time_t procEnd = time(NULL);
    struct tm *localEnd = localtime(&procEnd);
	char theTime[MAXSTR];
    strftime(theTime, MAXSTR, "%H:%M:%S", localEnd);
	cout << "Finished processing at " << theTime << "Run-time was " 
	          << difftime(procEnd, procStart) << " s." << endl;
	switch (components)
	{
		case 1: write_results<OneComponent>(outPrefix, paramsData, residualData, savedHeader); break;
		case 2: write_results<TwoComponent>(outPrefix, paramsData, residualData, savedHeader); break;
		case 3: write_results<ThreeComponent>(outPrefix, paramsData, residualData, savedHeader); break;
	}

	for (int i = 0; i < SPGR_files.size(); i++)
	{
		SPGR_files[i]->close();
		if (SPGR_B1_files[i]) SPGR_files[i]->close();
	}
	for (int i = 0; i < SSFP_files.size(); i++)
	{
		SSFP_files[i]->close();
		if (SSFP_B0_files[i]) SSFP_B0_files[i]->close();
		if (SSFP_B1_files[i]) SSFP_B1_files[i]->close();
	}
	free(M0Data);
	free(maskData);
	exit(EXIT_SUCCESS);
}

