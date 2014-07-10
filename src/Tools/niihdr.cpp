//
//  nifti_hdr.cpp
//  NiftiImage
//
//  Created by Tobias Wood on 08/07/2013.
//  Copyright (c) 2013 Tobias Wood. All rights reserved.
//

#include <memory>
#include <vector>
#include <iostream>
#include <getopt.h>

#include "Nifti/Nifti.h"

using namespace std;
using namespace Nifti;

const string usage = "nifti_hdr - A utility for getting information from Nifti headers.\n\
\n\
Usage: nifti_hdr [options] file1 [other files]\n\
By default the abbreviated header will be printed for each file.\n\
Abbreviated, full and compare modes are mutually exclusive.\n\
\n\
Options:\n\
	-d, --dims :  Print the size of each dimension for each file.\n\
	-v, --vox :   Print the voxel sizes for each file (with units). \n\
	-t, --trans : Print the XForm to physical space with precedence.\n\
	-a, --abbrev: Print the abbreviated header.\n\
	-f, --full:   Print the entire header.\n\
	-c, --comp:   Compare first file to all others and print message if the.\n\
	              physical spaces are different.\n\
	-h, --help:   Print this message and quit.\n\
";

enum Modes {
	Nothing = 0,
	Abbreviated,
	Full,
	Compare
};

static int mode = Nothing;
static int printDims = false, printVoxdims = false, printSize = false, printData = false,
           printTransform = false;

static struct option long_options[] =
{
	{"dims",   no_argument, &printDims, true},
	{"vox",    no_argument, &printVoxdims, true},
	{"form",   no_argument, &printTransform, true},
	{"size",   no_argument, &printSize, true},
	{"data",   no_argument, &printData, true},
	{"abbrev", no_argument, &mode, Abbreviated},
	{"full",   no_argument, &mode, Full},
	{"comp",   no_argument, &mode, Compare},
	{"help",   no_argument, 0, 'h'},
	{0, 0, 0, 0}
};

string voxMessage(const Header &hdr) {
	stringstream m;
	m << "Voxel sizes: " << hdr.voxDims().transpose() << " " << hdr.spaceUnits();
	if (hdr.voxDims().rows() > 3) {
		m << "/" << hdr.timeUnits();
	}
	return m.str();
}

string sizeMessage(const Header &hdr) {
	stringstream m;
	m << "Voxels per slice, per volume, total: "
      << hdr.dims().head(2).prod() << ", " << hdr.dims().head(3).prod() << ", " << hdr.dims().prod();
	return m.str();
}

string dataMessage(const Header &hdr) {
	stringstream m;
	m << "Datatype: " << TypeInfo(hdr.datatype()).name << ", size in bytes: " << TypeInfo(hdr.datatype()).size;
	return m.str();
}

int main(int argc, char **argv) {
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "dvtafchse", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'd': printDims = true; break;
			case 'v': printVoxdims = true; break;
			case 't': printTransform = true; break;
			case 's': printSize = true; break;
			case 'a': mode = Abbreviated; break;
			case 'f': mode = Full; break;
			case 'c': mode = Compare; break;
			case '?': // getopt will print an error message
			case 'h':
				cout << usage << endl;
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) <= 0 ) {
		cerr << "No input image file specified." << endl;
		exit(EXIT_FAILURE);
	}
	if (optind == 1) { // No options specified, default is short header
		mode = Abbreviated;
	}
	vector<File> images;
	try {
		images.reserve(argc - optind); // emplace_back can still trigger copies if the vector has to be resized
		for (;optind < argc; optind++) {
			images.emplace_back(argv[optind], Nifti::Mode::Read);
		}
	} catch (exception &e) {
		cerr << e.what() << endl;
	}

	if (mode == Compare) { // Compare first image to all others and check headers are compatible
		for (auto im = images.begin() + 1; im != images.end(); im++) {
			if (!images[0].header().matchesSpace(im->header())) {
				cout << "Header does not match against file: " << im->imagePath() << endl;
			}
		}
	}
	
	for (auto& im: images) {
		Header hdr = im.header();
		if (printDims) cout << im.dims().transpose() << endl;
		if (printVoxdims) cout << voxMessage(hdr) << endl;
		if (printTransform) cout << hdr.transform().matrix() << endl;
		if (printSize) cout << sizeMessage(hdr) << endl;
		if (printData) cout << dataMessage(hdr) << endl;
		
		if (mode == Abbreviated) {
			cout << "Short Nifti Header for file: " << im.imagePath() << endl;
			cout << "Dimensions:  " << hdr.dims().transpose() << endl;
			cout << dataMessage(hdr) << endl;
			cout << voxMessage(hdr) << endl;
			cout << "XForm matrix: " << endl << hdr.transform().matrix() << endl;
			(im.extensions().size() > 0) ? cout << "Has extensions." << endl : cout << "No extensions." << endl;
		} else if (mode == Full) {
			cout << "Full Nifti Header for file: " << im.imagePath() << endl;
			cout << dataMessage(hdr) << endl;
			cout << "Dimensions: " << im.dims().transpose() << endl;
			cout << voxMessage(hdr) << endl;
			
			cout << "Calibration (min, max): " << hdr.calibration_min << ", " << hdr.calibration_max << endl;
			cout << "Scaling (slope, inter): " << hdr.scaling_slope << ", " << hdr.scaling_inter << endl;
			cout << "Dimension labels (Phase, Freq, Slice):   " << hdr.phase_dim << ", " << hdr.freq_dim << ", " << hdr.slice_dim << endl;
			cout << "Slice info (Code, Start, End, Duration): " << ", " << hdr.slice_code << ", " << hdr.slice_start << ", " << hdr.slice_end << ", " << hdr.slice_duration << endl;
			cout << "Slice name: " << hdr.sliceName() << endl;
			cout << "Time offset: " << hdr.toffset << endl;
			
			cout << "Intent name:   " << hdr.intent_name << endl;
			cout << "Intent code:   " << hdr.intentName() << endl;
			cout << "Intent params: " << hdr.intent_p1 << ", " << hdr.intent_p2 << ", " << hdr.intent_p3 << endl;
			cout << "Description: " << hdr.description << endl;
			cout << "Aux File:    " << hdr.aux_file << endl;
			cout << "QForm: " << XFormName(hdr.qcode()) << endl;
			cout << hdr.qform().matrix() << endl;
			cout << "SForm: " << XFormName(hdr.scode()) << endl;
			cout << hdr.sform().matrix() << endl;
			cout << "Number of extensions: " << im.extensions().size() << endl;
		}
	}
	
	return EXIT_SUCCESS;
}

