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

const string usage = "niihdr - A utility for getting information from Nifti headers.\n\
\n\
Usage: nifti_hdr [options] file1 [other files]\n\
By default the abbreviated header will be printed for each file.\n\
Abbreviated, full and compare modes are mutually exclusive.\n\
\n\
Options:\n\
	-a, --abbrev: Print the abbreviated header (default).\n\
	-f, --full:   Print the entire header.\n\
	-c, --comp:   Compare first file to all others and print message if the.\n\
	              physical spaces are different.\n\
	-h, --help:   Print this message and quit.\n\
";

enum class Mode { Abbrev, Full, Cmp };

static Mode mode = Mode::Abbrev;

static struct option long_options[] = {
	{"abbrev", no_argument, 0, 'a'},
	{"full",   no_argument, 0, 'f'},
	{"comp",   no_argument, 0, 'c'},
	{"help",   no_argument, 0, 'h'},
	{0, 0, 0, 0}
};

int main(int argc, char **argv) {
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "afch", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'a': mode = Mode::Abbrev; break;
			case 'f': mode = Mode::Full; break;
			case 'c': mode = Mode::Cmp; break;
			case '?': // getopt will print an error message
			case 'h':
				cout << usage << endl;
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) <= 0 ) {
		cerr << "No input image file specified." << endl;
		cout << usage << endl;
		exit(EXIT_FAILURE);
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

	size_t did_not_match = 0;
	if (mode == Mode::Cmp) { // Compare first image to all others and check headers are compatible
		cout << "Checking headers against: " << images[0].headerPath() << endl;
		for (auto im = images.begin() + 1; im != images.end(); im++) {
			if (!images[0].header().matchesSpace(im->header())) {
				cout << im->headerPath() << " does not match." << endl;
				did_not_match++;
			}
		}
		cout << did_not_match << " headers did not match." << endl;
	}
	
	for (auto& im: images) {
		Nifti::Header hdr = im.header();
		if (mode == Mode::Abbrev) {
			cout << "Image path: " << im.imagePath() << endl;
			cout << "Datatype:   " << hdr.typeInfo().name << endl;
			cout << "Dimensions: " << hdr.dims().transpose() << " (rank " << to_string(hdr.rank()) << ")" << endl;
			cout << "Voxel size: " << hdr.voxDims().transpose() << endl;
			cout << "XForm: " << endl << hdr.transform().matrix() << endl;
			(im.extensions().size() > 0) ? cout << "Has extensions." << endl : cout << "No extensions." << endl;
		} else if (mode == Mode::Full) {
			cout << "Full Nifti Header for file: " << im.imagePath() << endl;
			cout << hdr << endl;
			cout << "Number of extensions: " << im.extensions().size() << endl;
		}
	}
	
	return EXIT_SUCCESS;
}

