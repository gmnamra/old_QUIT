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
#include <atomic>

using namespace std;

#include "procpar.h"
#include "DESPOT.h"
#include "DESPOT_Functors.h"
#include "RegionContraction.h"

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: mcdespot [options] input_file output_prefix\n\
\n\
The input file must consist of at least one line such as:\n\
SPGR path_to_spgr_file (path_to_B1_file) \n\
SSFP path_to_ssfp_file (path_to_B1) (path_to_B0) \n\
\n\
Phase-cycling is specified in degrees. If no B0/B1 correction is specified,\n\
default values of B0 = 0 Hz and B1 = 1 will be used. B0 value is only used if\n\
--B0 option specified, otherwise it is a free parameter.\n\
\n\
Options:\n\
	-v, --verbose     : Print extra information.\n\
	-m, --mask file   : Mask input with specified file.\n\
	--1, --2, --3     : Use 1, 2 or 3 component model (default 2).\n\
	--M0 file         : M0 Map file.\n\
	--B0              : Use the default or specified B0 value (don't fit).\n\
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

static int verbose = false, normalise = false, drop = false, fitB0 = true,
           start_slice = -1, end_slice = -1, slice = 0,
		   samples = 2500, retain = 25, contract = 10,
		   components = 2, tesla = 7, nP = 0;
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
	{"B0", no_argument, &fitB0, false},
	{0, 0, 0, 0}
};
//******************************************************************************
// SIGTERM interrupt handler - for ensuring data gets saved even on a ctrl-c
//******************************************************************************
NiftiImage savedHeader;
NiftiImage *paramsHdrs, residualHdr;
double **paramsData;
double *residualData;

void int_handler(int sig);
void int_handler(int sig)
{
	fprintf(stdout, "Processing terminated. Writing currently processed data.\n");
	for (int p = 0; p < nP; p++) {
		paramsHdrs[p].writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, paramsData[p]);
		paramsHdrs[p].close();
	}
	residualHdr.writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, residualData);
	residualHdr.close();
	exit(EXIT_FAILURE);
}


//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	#pragma mark Argument Processing
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
	#pragma mark  Read input file and set up corresponding SPGR & SSFP lists
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
		
		inHdr = new NiftiImage(path, NiftiImage::NIFTI_READ);
		if ((SPGR_files.size() == 0) && (SSFP_files.size() == 0))
			savedHeader = *inHdr;
		inHdr->checkCompatible(savedHeader);
		
		if (type == "SPGR") {
			cout << "Opened SPGR header: " << path << endl;
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
			if (thisLine >> path) { // Read a path to B1 file
				inHdr = new NiftiImage(path, NiftiImage::NIFTI_READ);
				inHdr->checkCompatible(savedHeader);
				cout << "Opened B1 correction header: " << path << endl;
			} else
				inHdr = NULL;
			SPGR_B1_files.push_back(inHdr);
		} else if (type == "SSFP") {
			cout << "Opened SSFP header: " << path << endl;
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
			if (thisLine >> path) { // Read a path to B1 file
				inHdr = new NiftiImage(path, NiftiImage::NIFTI_READ);
				inHdr->checkCompatible(savedHeader);
				cout << "Opened B1 correction header: " << path << endl;
			} else
				inHdr = NULL;
			SSFP_B1_files.push_back(inHdr);
			if (thisLine >> path) {	// Read a path to B0 file
				inHdr = new NiftiImage(path, NiftiImage::NIFTI_READ);
				inHdr->checkCompatible(savedHeader);
				cout << "Opened B0 correction header: " << path << endl;
			} else
				inHdr = NULL;
			SSFP_B0_files.push_back(inHdr);
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
	#pragma mark Select which angles to use in the analysis
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
	#pragma mark Allocate memory and set up boundaries.
	// Use NULL to indicate that default values should be used -
	// 0 for B0, 1 for B1
	//**************************************************************************
	int voxelsPerSlice = savedHeader.voxelsPerSlice();
	//int voxelsPerVolume = savedHeader.voxelsPerVolume();
	
	vector<double *> SPGR_vols(SPGR_files.size(), NULL), SPGR_B1s(SPGR_files.size(), NULL);
	vector<double *> SSFP_vols(SSFP_files.size(), NULL), SSFP_B0s(SSFP_files.size(), NULL), SSFP_B1s(SSFP_files.size(), NULL);
	for (int i = 0; i < SPGR_files.size(); i++) {
		SPGR_vols[i] = new double[voxelsPerSlice * nSPGR];
		if (SPGR_B1_files[i]) SPGR_B1s[i] = new double[voxelsPerSlice * nSPGR];
	}
	for (int i = 0; i < SSFP_files.size(); i++) {
		SSFP_vols[i] = new double[voxelsPerSlice * nSSFP];
		if (SSFP_B0_files[i]) SSFP_B0s[i] = new double[voxelsPerSlice * nSSFP];
		if (SSFP_B1_files[i]) SSFP_B1s[i] = new double[voxelsPerSlice * nSSFP];
	}
	
	cout << "Using " << components << " component model." << endl;
	switch (components)
	{
		case 1: nP = OneComponent::nP; break;
		case 2: nP = TwoComponent::nP; break;
		case 3: nP = ThreeComponent::nP; break;
	}
		
	VectorXd loBounds(nP), hiBounds(nP);
	VectorXi loConstraints(nP), hiConstraints(nP);
	
	residualData = new double[voxelsPerSlice];
	paramsData = new double *[nP];
	if (tesla < 0)
		cout << "Enter " << nP << " parameter pairs (low then high): " << flush;
	paramsHdrs = new NiftiImage[nP];
	residualHdr = savedHeader;
	residualHdr.setnt(1); residualHdr.setDatatype(NIFTI_TYPE_FLOAT32);
	residualHdr.open(outPrefix + "_residual.nii.gz", NiftiImage::NIFTI_WRITE);
	for (int i = 0; i < nP; i++) {
		paramsData[i] = new double[voxelsPerSlice];
		paramsHdrs[i] = savedHeader;
		paramsHdrs[i].setnt(1); paramsHdrs[i].setDatatype(NIFTI_TYPE_FLOAT32);
		switch (components) {
			case 1:	paramsHdrs[i].open(outPrefix + "_" + OneComponent::names[i] + ".nii.gz", NiftiImage::NIFTI_WRITE); break;
			case 2:	paramsHdrs[i].open(outPrefix + "_" + TwoComponent::names[i] + ".nii.gz", NiftiImage::NIFTI_WRITE); break;
			case 3:	paramsHdrs[i].open(outPrefix + "_" + ThreeComponent::names[i] + ".nii.gz", NiftiImage::NIFTI_WRITE); break;
		}
		loConstraints[i] = true; hiConstraints[i] = true;
		if (tesla == 3) {
			switch (components) {
				case 1: loBounds[i] = OneComponent::lo3Bounds[i]; hiBounds[i] = OneComponent::hi3Bounds[i]; break;
				case 2: loBounds[i] = TwoComponent::lo3Bounds[i]; hiBounds[i] = TwoComponent::hi3Bounds[i]; break;
				case 3: loBounds[i] = ThreeComponent::lo3Bounds[i]; hiBounds[i] = ThreeComponent::hi3Bounds[i]; break;
			}
		} else if (tesla == 7) {
			switch (components) {
				case 1: loBounds[i] = OneComponent::lo7Bounds[i]; hiBounds[i] = OneComponent::hi7Bounds[i]; break;
				case 2: loBounds[i] = TwoComponent::lo7Bounds[i]; hiBounds[i] = TwoComponent::hi7Bounds[i]; break;
				case 3: loBounds[i] = ThreeComponent::lo7Bounds[i]; hiBounds[i] = ThreeComponent::hi7Bounds[i]; break;
			}
		} else {
			cin >> loBounds[i] >> hiBounds[i];
		}
	}
	// If normalising, don't bother fitting for M0
	if (normalise) {
		loBounds[0] = 1.;
		hiBounds[0] = 1.;
	}
	// If fitting, give a suitable range. Otherwise fix and let the functors
	// pick up the specified value
	if (fitB0) {
		loBounds[nP - 1] = -0.5 / ssfpTR;
		hiBounds[nP - 1] =  0.5 / ssfpTR;
	} else {
		loBounds[nP - 1] = 0.;
		hiBounds[nP - 1] = 0.;
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
	
	for (slice = start_slice; slice < end_slice; slice++)
	{
		if (verbose) cout << "Reading data for slice " << slice << "..." << flush;
		atomic<int> voxCount{0};
		const int sliceOffset = slice * voxelsPerSlice;
		
		// Read data for slices
		for (int i = 0; i < SPGR_files.size(); i++) {
			SPGR_files[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, SPGR_vols[i]);
			if (SPGR_B1_files[i]) SPGR_B1_files[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, SPGR_B1s[i]);
		}
		for (int i = 0; i < SSFP_files.size(); i++) {
			SSFP_files[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, SSFP_vols[i]);
			if (SSFP_B0_files[i]) SSFP_B0_files[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, SSFP_B0s[i]);
			if (SSFP_B1_files[i]) SSFP_B1_files[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, SSFP_B1s[i]);
		}
		if (verbose) cout << "processing..." << flush;
		clock_t loopStart = clock();
		function<void (const int&)> processVox = [&] (const int &vox)
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
				// Add the voxel number to the time to get a decent random seed
				int rSeed = static_cast<int>(time(NULL)) + vox;
				switch (components)
				{
					case 1: {
						OneComponent tc(spgrAngles, SPGR_signals, SPGR_B1, spgrTR,
						                ssfpAngles, ssfpPhases, SSFP_signals, SSFP_B0, SSFP_B1, ssfpTR,
										normalise, fitB0);
						residual = regionContraction<OneComponent>(params, tc, localLo, localHi,
															       loConstraints, hiConstraints,
															       samples, retain, contract, 0.05, expand, rSeed);
						break; }
					case 2: {
						TwoComponent tc(spgrAngles, SPGR_signals, SPGR_B1, spgrTR,
						                ssfpAngles, ssfpPhases, SSFP_signals, SSFP_B0, SSFP_B1, ssfpTR,
										normalise, fitB0);
						residual = regionContraction<TwoComponent>(params, tc, localLo, localHi,
															       loConstraints, hiConstraints,
															       samples, retain, contract, 0.05, expand, rSeed);
						break; }
					case 3: {
						ThreeComponent tc(spgrAngles, SPGR_signals, SPGR_B1, spgrTR,
						                ssfpAngles, ssfpPhases, SSFP_signals, SSFP_B0, SSFP_B1, ssfpTR,
										normalise, fitB0);
						residual = regionContraction<ThreeComponent>(params, tc, localLo, localHi,
																 	loConstraints, hiConstraints,
																 	samples, retain, contract, 0.05, expand, rSeed);
						break; }
				}
			}
			for (int p = 0; p < nP; p++)
				paramsData[p][vox]  = params[p];
			residualData[vox] = residual;
		};
		apply_for(voxelsPerSlice, processVox);
		
		if (verbose) {
			clock_t loopEnd = clock();
			if (voxCount > 0)
				cout << voxCount << " unmasked voxels, CPU time per voxel was "
				          << ((loopEnd - loopStart) / ((float)voxCount * CLOCKS_PER_SEC)) << " s, ";
			cout << "finished." << endl;
		}
		
		for (int p = 0; p < nP; p++)
			paramsHdrs[p].writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, paramsData[p]);
		residualHdr.writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, residualData);
	}
    time_t procEnd = time(NULL);
    struct tm *localEnd = localtime(&procEnd);
	char theTime[MAXSTR];
    strftime(theTime, MAXSTR, "%H:%M:%S", localEnd);
	cout << "Finished processing at " << theTime << "Run-time was " 
	          << difftime(procEnd, procStart) << " s." << endl;
	
	// Clean up memory and close files (automatically done in destructor)
	delete[] paramsHdrs;
	for (int p = 0; p < nP; p++)
		delete[] paramsData[p];
	delete[] paramsData;
	delete[] residualData;
	for (int i = 0; i < SPGR_files.size(); i++)	{
		delete[] SPGR_vols[i];
		if (SPGR_B1_files[i])
			delete[] SPGR_B1s[i];
	}
	for (int i = 0; i < SSFP_files.size(); i++)
	{
		delete[] SSFP_vols[i];
		if (SSFP_B0_files[i])
			delete[] SSFP_B0s[i];
		if (SSFP_B1_files[i])
			delete[] SSFP_B1s[i];
	}
	delete[] M0Data;
	delete[] maskData;
	exit(EXIT_SUCCESS);
}

