//
//  niihdr.cpp
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

const string usage = "niihdr - A utility for getting information from Nifti headers.\n\
\n\
Usage: niihdr [options] file1 [other files]\n\
By default the abbreviated header will be printed for each file.\n\
Alternatively the entire header or specific fields can be printed.\n\
\n\
Options:\n\
	-a, --abbrev : Print the abbreviated header (default).\n\
	-f, --full   : Print the entire header.\n\
	-h, --help   : Print this message and quit.\n\
	--voxdims    : Print the voxel dimensions.\n\
	--tr         : Print the TR.\n\
	--xform      : Print the highest priority transform.\n\
";


static int printAbbrev = false, printFull = false,
           printVoxdims = false, printTR = false, printXform = false;
static struct option long_options[] = {
	{"abbrev",  no_argument, 0, 'a'},
	{"full",    no_argument, 0, 'f'},
	{"help",    no_argument, 0, 'h'},
	{"voxdims", no_argument, &printVoxdims, 1},
	{"tr",      no_argument, &printTR, 1},
	{"xform",   no_argument, &printXform, 1},
	{0, 0, 0, 0}
};

int main(int argc, char **argv) {
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "afch", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'a': printAbbrev = true; printFull = false; break;
			case 'f': printFull = true; printAbbrev = false; break;
			case '?': // getopt will print an error message
			case 'h':
				cout << usage << endl;
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) <= 0 ) {
		cerr << "No input image file specified." << endl;
		cout << usage << endl;
		return EXIT_FAILURE;
	} else if (optind == 1) { // No other arguments specified
		printAbbrev = true;
	}
	vector<Nifti::File> images;
	try {
		images.reserve(argc - optind); // emplace_back can still trigger copies if the vector has to be resized
		for (;optind < argc; optind++) {
			images.emplace_back(argv[optind]);
		}
	} catch (exception &e) {
		cerr << e.what() << endl;
	}

	for (auto& im: images) {
		Nifti::Header hdr = im.header();
		if (printAbbrev) {
			cout << "Image path: " << im.imagePath() << endl;
			cout << "Datatype:   " << hdr.typeInfo().name << endl;
			cout << "Dimensions: " << hdr.dims().transpose() << " (rank " << to_string(hdr.rank()) << ")" << endl;
			cout << "Voxel size: " << hdr.voxDims().transpose() << endl;
			cout << "XForm: " << endl << hdr.transform().matrix() << endl;
			(im.extensions().size() > 0) ? cout << "Has extensions." << endl : cout << "No extensions." << endl;
		}
		if (printFull) {
			cout << "Full Nifti Header for file: " << im.imagePath() << endl;
			cout << hdr << endl;
			cout << "Number of extensions: " << im.extensions().size() << endl;
		}
		if (printVoxdims) {
			if (hdr.rank() > 3) {
				cout << hdr.voxDims().head(3).transpose() << endl;
			} else {
				cout << hdr.voxDims() << endl;
			}
		}
		if (printTR) { cout << hdr.voxDim(4) << endl; }
		if (printXform) { cout << hdr.transform().matrix() << endl; }
	}
	
	return EXIT_SUCCESS;
}

