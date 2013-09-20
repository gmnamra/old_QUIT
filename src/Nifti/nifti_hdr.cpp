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

#include "Nifti.h"

using namespace std;

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
           printTransform = false, printExtensions = false;

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
	{"ext",    no_argument, 0, 'e'},
	{0, 0, 0, 0}
};

string voxMessage(const Nifti &im) {
	stringstream m;
	m << "Voxel sizes: " << im.voxDims().transpose() << " " << im.spaceUnits();
	if (im.voxDims().rows() > 3) {
		m << "/" << im.timeUnits();
	}
	return m.str();
}

string sizeMessage(const Nifti &im) {
	stringstream m;
	m << "Voxels per slice, per volume, total: "
      << im.voxelsPerSlice() << ", " << im.voxelsPerVolume() << ", " << im.voxelsTotal();
	return m.str();
}

string dataMessage(const Nifti &im) {
	stringstream m;
	m << "Datatype: " << Nifti::TypeInfo(im.datatype()).name << ", size in bytes: " << Nifti::TypeInfo(im.datatype()).size;
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
			case 'e': printExtensions = true; break;
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
	vector<Nifti> images;
	images.reserve(argc - optind); // emplace_back can still trigger copies if the vector has to be resized
	for (;optind < argc; optind++) {
		images.emplace_back(argv[optind], Nifti::Mode::Read);
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
		if (printTransform) cout << im.transform().matrix() << endl;
		if (printSize) cout << sizeMessage(im) << endl;
		if (printData) cout << dataMessage(im) << endl;
		
		if (mode == Abbreviated) {
			cout << "Short Nifti Header for file: " << im.imagePath() << endl;
			cout << "Dimensions:  " << im.dims().transpose() << endl;
			cout << voxMessage(im) << endl;
			cout << "XForm matrix: " << endl << im.transform().matrix() << endl;
			cout << "Number of extensions: " << im.extensions().size() << endl;
		} else if (mode == Full) {
			cout << "Full Nifti Header for file: " << im.imagePath() << endl;
			cout << dataMessage(im) << endl;
			cout << "Dimensions: " << im.dims().transpose() << endl;
			cout << voxMessage(im) << endl;
			
			cout << "Calibration (min, max): " << im.calibration_min << ", " << im.calibration_max << endl;
			cout << "Scaling (slope, inter): " << im.scaling_slope << ", " << im.scaling_inter << endl;
			cout << "Dimension labels (Phase, Freq, Slice):   " << im.phase_dim << ", " << im.freq_dim << ", " << im.slice_dim << endl;
			cout << "Slice info (Code, Start, End, Duration): " << ", " << im.slice_code << ", " << im.slice_start << ", " << im.slice_end << ", " << im.slice_duration << endl;
			cout << "Slice name: " << im.sliceName() << endl;
			cout << "Time offset: " << im.toffset << endl;
			
			cout << "Intent name:   " << im.intent_name << endl;
			cout << "Intent code:   " << im.intentName() << endl;
			cout << "Intent params: " << im.intent_p1 << ", " << im.intent_p2 << ", " << im.intent_p3 << endl;
			cout << "Description: " << im.description << endl;
			cout << "Aux File:    " << im.aux_file << endl;
			cout << "QForm: " << Nifti::XFormName(im.qcode()) << endl;
			cout << im.qform().matrix() << endl;
			cout << "SForm: " << Nifti::XFormName(im.scode()) << endl;
			cout << im.sform().matrix() << endl;
			cout << "Extensions: " << endl;
			for (auto &e : im.extensions()) {
				cout << "Extension Code: " << e.code() << endl;
				string out(e.data().begin(), e.data().end());
				cout << out << endl;
			}
		}
		vector<float> data(8);
		Nifti::ArrayXs start(3); Nifti::ArrayXs size(3);
		start << 4, 2, 10;
		size << 4, 2, 1;
		cout << "Before size: " << data.size() << endl;
		im.readWriteVoxels(start, size, data);
		cout << "After size: " << data.size() << endl;
		for (size_t i = 0; i < 8; i++) cout << data.at(i) << "\t";
	}
	
	return EXIT_SUCCESS;
}

