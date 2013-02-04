/*
 *  mcdespot_main.cpp
 *
 *  Created by Tobias Wood on 14/02/2012.
 *  Copyright 2012 Tobias Wood. All rights reserved.
 *
 */

#include <time.h>
#include <getopt.h>
#include <signal.h>
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
"mcdespot - written by tobias.wood@kcl.ac.uk. \n\
Acknowledgements greatfully received, grant discussions welcome."
};

const string usage {
"Usage is: mcdespot [options]\n\
\n\
The input must consist of at least one line such as:\n\
SPGR path_to_spgr_file (path_to_B1_file) \n\
SSFP path_to_ssfp_file (path_to_B1) (path_to_B0) \n\
\n\
If -d or --drop is specified further input will be prompted for from stdin.\n\
\n\
Phase-cycling is specified in degrees. If no B0/B1 correction is specified,\n\
default values of B0 = 0 Hz and B1 = 1 will be used. B0 value is only used if\n\
--B0 option specified, otherwise it is a free parameter.\n\
\n\
Options:\n\
	--help, -h		  : Print this message.\n\
	--verbose, -v     : Print extra information.\n\
	--no-prompt, -p   : Don't print prompts for input.\n\
	--mask, -m file   : Mask input with specified file.\n\
	--out, -o path    : Add a prefix to the output filenames.\n\
	--1, --2, --3     : Use 1, 2 or 3 component model (default 2).\n\
	--M0 file         : M0 Map file.\n\
	--B0              : Use the default or specified B0 value (don't fit).\n\
	--start_slice n   : Only start processing at slice n.\n\
	--end_slice n     : Finish at slice n-1.\n\
	--drop, -d        : Drop certain flip-angles (Read from stdin).\n\
	--normalise, -n   : Normalise signals to maximum (Ignore M0).\n\
	--samples, -s n   : Use n samples for region contraction (Default 5000).\n\
	--retain, -r  n   : Retain n samples for new boundary (Default 50).\n\
	--contract, -c n  : Contract a maximum of n times (Default 10).\n\
	--expand, -e n    : Re-expand boundary by percentage n (Default 0).\n\
	--tesla, -t 3     : Boundaries suitable for 3T\n\
	            7     : Boundaries suitable for 7T (default)\n\
	            u     : User specified boundaries from stdin.\n"
};

static int verbose = false, prompt = true,
           normalise = false, drop = false, fitB0 = true,
           start_slice = -1, end_slice = -1, slice = 0,
		   samples = 5000, retain = 50, contract = 10,
		   components = 2, tesla = -1;
static double expand = 0.;
static string outPrefix;
static struct option long_options[] =
{
	{"M0", required_argument, 0, 'M'},
	{"mask", required_argument, 0, 'm'},
	{"out", required_argument, 0, 'o'},
	{"verbose", no_argument, 0, 'v'},
	{"no-prompt", no_argument, 0, 'p'},
	{"start_slice", required_argument, 0, 'S'},
	{"end_slice", required_argument, 0, 'E'},
	{"normalise", no_argument, &normalise, true},
	{"tesla", required_argument, 0, 't'},
	{"samples", required_argument, 0, 's'},
	{"retain", required_argument, 0, 'r'},
	{"contract", required_argument, 0, 'c'},
	{"expand", required_argument, 0, 'e'},
	{"drop", no_argument, &drop, true},
	{"help", no_argument, 0, 'h'},
	{"1", no_argument, &components, 1},
	{"2", no_argument, &components, 2},
	{"3", no_argument, &components, 3},
	{"B0", no_argument, &fitB0, false},
	{0, 0, 0, 0}
};
//******************************************************************************
#pragma mark SIGTERM interrupt handler - for ensuring data is saved on a ctrl-c
//******************************************************************************
NiftiImage *paramsHdrs, residualHdr;
vector<double *> paramsData;
double *residualData;

void int_handler(int sig);
void int_handler(int sig)
{
	fprintf(stdout, "Processing terminated. Writing currently processed data.\n");
	for (int p = 0; p < mcDESPOT::nP(components); p++) {
		paramsHdrs[p].writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, paramsData[p]);
		paramsHdrs[p].close();
	}
	residualHdr.writeSubvolume(0, 0, slice, 0, -1, -1, slice + 1, 1, residualData);
	residualHdr.close();
	exit(EXIT_FAILURE);
}

//******************************************************************************
#pragma mark Read in all required files and data from cin
//******************************************************************************
NiftiImage *parseInput(vector<mcDESPOT::SignalType> &signalTypes, vector<VectorXd> &angles,
                       vector<mcDESPOT::ModelConstants> &consts,
				       vector<NiftiImage *> &signalFiles,
				       vector<NiftiImage *> &B1_files,
				       vector<NiftiImage *> &B0_files);
NiftiImage *parseInput(vector<mcDESPOT::SignalType> &signalTypes, vector<VectorXd> &angles,
                       vector<mcDESPOT::ModelConstants> &consts,
				       vector<NiftiImage *> &signalFiles,
				       vector<NiftiImage *> &B1_files,
				       vector<NiftiImage *> &B0_files) {
	NiftiImage *inHdr, *savedHeader;
	string type, path;
	if (prompt) cout << "Specify next image type (SPGR/SSFP): " << flush;
	while (getline(cin, type) && !type.empty()) {
		
		if (type == "SPGR")
			signalTypes.push_back(mcDESPOT::SignalSPGR);
		else if (type == "SSFP")
			signalTypes.push_back(mcDESPOT::SignalSSFP);
		else {
			cerr << "Unknown signal type: " << type << endl;
			exit(EXIT_FAILURE);
		}
		if (prompt) cout << "Enter image path: " << flush;
		getline(cin, path);
		inHdr = new NiftiImage(path, NiftiImage::NIFTI_READ);
		if (verbose) cout << "Opened " << type << " header: " << path << endl;
		if (signalFiles.size() == 0) // First time through save the header for comparison
			savedHeader = inHdr;
		inHdr->checkVoxelsCompatible(*savedHeader);
		signalFiles.push_back(inHdr);
		double inTR = 0., inPhase = 0.;
		VectorXd inAngles(inHdr->dim(4));
		#ifdef USE_PROCPAR
		ParameterList pars;
		if (ReadProcpar(inHdr->basename() + ".procpar", pars)) {
			inTR = RealValue(pars, "tr");
			for (int i = 0; i < inAngles.size(); i++)
				inAngles[i] = RealValue(pars, "flip1", i);
			if (signalTypes.back() == mcDESPOT::SignalSSFP)
				inPhase = RealValue(pars, "rfphase");
		} else
		#endif
		{
			if (prompt) cout << "Enter TR: " << flush;
			cin >> inTR;
			if (prompt) cout << "Enter Flip-angles: " << flush;
			for (int i = 0; i < inAngles[i]; i++)
				cin >> inAngles[i];
			getline(cin, path); // Just to eat the newline
			if (signalTypes.back() == mcDESPOT::SignalSSFP) {
				if (prompt) cout << "Enter SSFP Phase-Cycling: ";
				cin >> inPhase;
			}
		}
		consts.push_back( { inTR, inPhase * M_PI / 180., 0., 0. } );
		angles.push_back(inAngles * M_PI / 180.);
		
		inHdr = NULL;
		if (prompt) cout << "Enter B1 Map Path (blank for none): " << flush;
		getline(cin, path);
		if (!path.empty()) {
			inHdr = new NiftiImage(path, NiftiImage::NIFTI_READ);
			inHdr->checkVoxelsCompatible(*savedHeader);
			if (verbose) cout << "Opened B1 correction header: " << path << endl;
		}
		B1_files.push_back(inHdr);
		
		inHdr = NULL;
		if (signalTypes.back() == mcDESPOT::SignalSSFP) {
			if (prompt) cout << "Enter path to B0 map: " << flush;
			getline(cin, path);
			if (!path.empty()) {
				inHdr = new NiftiImage(path, NiftiImage::NIFTI_READ);
				inHdr->checkVoxelsCompatible(*savedHeader);
				if (verbose) cout << "Opened B0 correction header: " << path << endl;
			}
		}
		B0_files.push_back(inHdr);
		// Print message ready for next loop
		if (prompt) cout << "Specify next image type (SPGR/SSFP): " << flush;
	}
	if (signalTypes.size() == 0) {
		cerr << "No input images specified." << endl;
		exit(EXIT_FAILURE);
	}
	return savedHeader;
}
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	#pragma mark Argument Processing
	//**************************************************************************
	cout << credit << endl;
	Eigen::initParallel();
	NiftiImage inHeader, *savedHeader;
	double *maskData = NULL, *PDData = NULL;
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvpt:m:o:znds:r:c:e:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'm':
				cout << "Reading mask file " << optarg << endl;
				inHeader.open(optarg, NiftiImage::NIFTI_READ);
				maskData = inHeader.readVolume<double>(0);
				inHeader.close();
				break;
			case 'M':
				cout << "Reading PD file " << optarg << endl;
				inHeader.open(optarg, NiftiImage::NIFTI_READ);
				PDData = inHeader.readVolume<double>(0);
				inHeader.close();
				break;
			case 'o':
				outPrefix = optarg;
				cout << "Output prefix will be: " << outPrefix << endl;
				break;
			case 'v': verbose = true; break;
			case 'p': prompt = false; break;
			case 'S': start_slice = atoi(optarg); break;
			case 'E': end_slice = atoi(optarg); break;
			case 'n': normalise = true; break;
			case 'd': drop = true; break;
			case 's': samples  = atoi(optarg); break;
			case 'r': retain   = atoi(optarg); break;
			case 'c': contract = atoi(optarg); break;
			case 'e': expand   = atof(optarg); break; 
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
			case 0:
				// Just a flag
				break;
			case 'h':
			case '?': // getopt will print an error message
				cout << usage << endl;
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != 0) {
		cerr << usage << endl << "Incorrect number of arguments." << endl;
		return EXIT_FAILURE;
	}

	//**************************************************************************
	#pragma mark  Read input and set up corresponding SPGR & SSFP lists
	//**************************************************************************
	vector<mcDESPOT::SignalType> signalTypes;
	vector<mcDESPOT::ModelConstants> consts;
	vector<VectorXd> angles;
	vector<NiftiImage *> signalFiles, B1_files, B0_files;
	savedHeader = parseInput(signalTypes, angles, consts, signalFiles, B1_files, B0_files);
	//**************************************************************************
	#pragma mark Allocate memory and set up boundaries.
	// Use NULL to indicate that default values should be used -
	// 0 for B0, 1 for B1
	//**************************************************************************
	int voxelsPerSlice = savedHeader->voxelsPerSlice();
	//int voxelsPerVolume = savedHeader.voxelsPerVolume();
	
	vector<double *> signalVolumes(signalFiles.size(), NULL),
	                 B1Volumes(signalFiles.size(), NULL),
				     B0Volumes(signalFiles.size(), NULL);
	for (int i = 0; i < signalFiles.size(); i++) {
		signalVolumes[i] = new double[voxelsPerSlice * angles[i].size()];
		if (B1_files[i]) B1Volumes[i] = new double[voxelsPerSlice];
		if (B0_files[i]) B0Volumes[i] = new double[voxelsPerSlice];
	}
	
	cout << "Using " << components << " component model." << endl;
	const int nP = mcDESPOT::nP(components);
	// The lo/hiBounds methods will make sure these vectors are the right
	// length, even for tesla == -1
	ArrayXd loBounds, hiBounds;
	loBounds = mcDESPOT::defaultLo(components, tesla);
	hiBounds = mcDESPOT::defaultHi(components, tesla);
	vector<string> names = mcDESPOT::names(components);
	residualData = new double[voxelsPerSlice];
	residualHdr = *savedHeader;
	residualHdr.setDim(4, 1); residualHdr.setDatatype(NIFTI_TYPE_FLOAT32);
	residualHdr.open(outPrefix + "MCD_Residual.nii.gz", NiftiImage::NIFTI_WRITE);
	
	paramsData.resize(nP);
	paramsHdrs = new NiftiImage[nP];
	if (prompt && (tesla < 0))
		cout << "Enter " << nP << " parameter pairs (low then high): " << flush;
	for (int i = 0; i < nP; i++) {
		if (tesla < 0)
			cin >> loBounds[i] >> hiBounds[i];
		paramsData[i] = new double[voxelsPerSlice];
		paramsHdrs[i] = *savedHeader;
		paramsHdrs[i].setDim(4, 1); paramsHdrs[i].setDatatype(NIFTI_TYPE_FLOAT32);
		paramsHdrs[i].open(outPrefix + "MCD_" + names[i] + ".nii.gz", NiftiImage::NIFTI_WRITE);
	}
	
	// If fitting, give a suitable range. Otherwise fix and let the functors
	// pick up the specified value (for now assume the last specified file is
	// an SSFP sequence, needs fixing / generalising)
	if (fitB0) {
		loBounds[0] = -0.5 / consts.back().TR;
		hiBounds[0] =  0.5 / consts.back().TR;
	} else {
		loBounds[0] = 0.;
		hiBounds[0] = 0.;
	}
	// If normalising, don't bother fitting for PD
	if (normalise) {
		loBounds[1] = 1.;
		hiBounds[1] = 1.;
	}
	
	if (verbose) {
		cout << "Low bounds: " << loBounds.transpose() << endl;
		cout << "Hi bounds:  " << hiBounds.transpose() << endl;
	}
	//**************************************************************************
	#pragma mark Do the fitting
	//**************************************************************************
    time_t procStart = time(NULL);
	if ((start_slice < 0) || (start_slice >= savedHeader->dim(3)))
		start_slice = 0;
	if ((end_slice < 0) || (end_slice > savedHeader->dim(3)))
		end_slice = savedHeader->dim(3);
	signal(SIGINT, int_handler);	// If we've got here there's actually allocated data to save
	cout << "Starting processing." << endl;
	
	for (slice = start_slice; slice < end_slice; slice++)
	{
		if (verbose) cout << "Reading data for slice " << slice << "..." << flush;
		atomic<int> voxCount{0};
		const int sliceOffset = slice * voxelsPerSlice;
		
		// Read data for slices
		for (int i = 0; i < signalFiles.size(); i++) {
			signalFiles[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, signalVolumes[i]);
			if (B1_files[i]) B1_files[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, B1Volumes[i]);
			if (B0_files[i]) B0_files[i]->readSubvolume<double>(0, 0, slice, 0, -1, -1, slice + 1, -1, B0Volumes[i]);
		}
		if (verbose) cout << "processing..." << flush;
		clock_t loopStart = clock();
		function<void (const int&)> processVox = [&] (const int &vox)
		{
			double residual = 0.;
			ArrayXd params(nP); params.setZero();
			if ((!maskData || (maskData[sliceOffset + vox] > 0.)) &&
			    (!PDData || (PDData[sliceOffset + vox] > 0.))) {
				voxCount++;
				
				vector<VectorXd> signals(signalFiles.size());
				for (int i = 0; i < signalFiles.size(); i++) {
					VectorXd temp(angles[i].size());
					for (int j = 0; j < angles[i].size(); j++) {
						temp[j] = signalVolumes[i][voxelsPerSlice*j + vox];
					}
					if (normalise)
						temp /= temp.mean();
					signals[i] = temp;
					consts[i].B0 = B0Volumes[i] ? B0Volumes[i][vox] : 0.;
					consts[i].B1 = B1Volumes[i] ? B1Volumes[i][vox] : 1.;
				}
				// Add the voxel number to the time to get a decent random seed
				int rSeed = static_cast<int>(time(NULL)) + vox;
				ArrayXd localLo = loBounds, localHi = hiBounds, localP(nP);
				if (PDData) {
					localLo(0) = (double)PDData[sliceOffset + vox];
					localHi(0) = (double)PDData[sliceOffset + vox];
				}
				mcDESPOT mcd(components, signalTypes, angles, signals, consts, normalise, fitB0);
				residual = regionContraction<mcDESPOT>(localP, mcd, localLo, localHi,
															 samples, retain, contract, 0.05, expand, rSeed);
				params = localP;
			}
			for (int p = 0; p < nP; p++) {
				paramsData[p][vox]  = params[p];
			}
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
	char theTime[512];
    strftime(theTime, 512, "%H:%M:%S", localEnd);
	cout << "Finished processing at " << theTime << ". Run-time was " 
	          << difftime(procEnd, procStart) << " s." << endl;
	
	// Clean up memory and close files (automatically done in destructor)
	residualHdr.close();
	delete[] paramsHdrs;
	for (int p = 0; p < nP; p++)
		delete[] paramsData[p];
	delete[] residualData;
	for (int i = 0; i < signalFiles.size(); i++)	{
		delete[] signalVolumes[i];
		if (B1Volumes[i])
			delete[] B1Volumes[i];
		if (B0Volumes[i])
			delete[] B0Volumes[i];
	}
	delete[] PDData;
	delete[] maskData;
	return EXIT_SUCCESS;
}

