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

using namespace std;

#include "NiftiImage.h"

const string usage = "nifti_hdr - A utility for getting information from Nifti headers.\n\
\n\
Usage: nifti_hdr [options] file1 [other files]\n\
By default the abbreviated header will be printed for each file.\n\
Abbreviated, full and compare modes are mutually exclusive.\n\
\n\
Options:\n\
	-d, --dims :  Print the size of each dimension for each file.\n\
	-v, --vox :   Print the voxel sizes for each file (with units). \n\
	-t, --trans : Print the transform to physical space with precedence.\n\
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
static int printDims = false, printVoxdims = false, printTransform = false;

static struct option long_options[] =
{
	{"dims",   no_argument, &printDims, true},
	{"vox",    no_argument, &printVoxdims, true},
	{"trans",  no_argument, &printTransform, true},
	{"abbrev", no_argument, &mode, Abbreviated},
	{"full",   no_argument, &mode, Full},
	{"comp",   no_argument, &mode, Compare},
	{"help",   no_argument, 0, 'h'},
	{0, 0, 0, 0}
};

string voxMessage(const NiftiImage& im) {
	stringstream m;
	m << "Voxel sizes: " << im.voxDims().transpose() << " " << im.spaceUnits();
	if (im.voxDims().rows() > 3) {
		m << "/" << im.timeUnits();
	}
	return m.str();
}

int main(int argc, char **argv) {
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "dvtafch", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'd': printDims = true; break;
			case 'v': printVoxdims = true; break;
			case 't': printTransform = true; break;
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
	vector<NiftiImage> images;
	images.reserve(argc - optind); // emplace_back can still trigger copies if the vector has to be resized
	for (;optind < argc; optind++) {
		images.emplace_back(argv[optind], NiftiImage::READ_HEADER);
	}
	
	if (mode == Compare) { // Compare first image to all others and check headers are compatible
		for (auto im = images.begin() + 1; im != images.end(); im++) {
			if (!images[0].matchesSpace(*im)) {
				cout << "Header does not match against file: " << im->imagePath() << endl;
			}
		}
	}
	
	for (auto& im: images) {
		if (printDims) cout << im.dims().transpose() << endl;
		if (printVoxdims) cout << voxMessage(im) << endl;
		if (printTransform) cout << im.ijk_to_xyz() << endl;
		
		if (mode == Abbreviated) {
			cout << "Short Nifti Header for file: " << im.imagePath() << endl;
			cout << "Dimensions:  " << im.dims().transpose() << endl;
			cout << voxMessage(im) << endl;
			cout << "Transform matrix: " << endl << im.ijk_to_xyz() << endl;
		} else if (mode == Full) {
			cout << "Full Nifti Header for file: " << im.imagePath() << endl;
		}
	}
	
	return EXIT_SUCCESS;
}
