//
//  nifti_hdr.cpp
//  NiftiImage
//
//  Created by Tobias Wood on 08/07/2013.
//  Copyright (c) 2013 Tobias Wood. All rights reserved.
//

#include <iostream>
#include <getopt.h>

using namespace std;

#include "NiftiImage.h"

const string usage = "nifti_hdr - A utility for getting information from Nifti headers.";

static struct option long_options[] =
{
	{0, 0, 0, 0}
};

enum Modes {
	Abbreviated = 0,
	Full,
	Compare
};

int main(int argc, char **argv) {
	
	int mode = Abbreviated;
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "c", long_options, &indexptr)) != -1) {
		switch (c) {
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
	
	NiftiImage first(argv[optind++], NiftiImage::READ);
	
	if (mode == Abbreviated) {
		cout << "Short Nifti Header for file: " << first.basename() << endl;
		cout << "Dimensions:  " << first.dims().transpose() << endl;
		cout << "Voxel Sizes: " << first.voxDims().transpose() << endl;
		cout << "Transform matrix: " << endl << first.ijk_to_xyz() << endl;
	} else if (mode == Full) {
		cout << "Full Nifti Header for file: " << first.basename() << endl;
	}
	
	return EXIT_SUCCESS;
}
