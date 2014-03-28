/*
 *  crop_main.cpp
 *
 *  Created by Tobias Wood on 02/05/2013.
 *  Copyright (c) 2013 Tobias Wood.
 *
 *  Based in part on work by Sean Deoni
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include "Nifti.h"

const char *usage = "Usage is: crop input_file\n\
\n\
Output is compatible with fslroi";
	
int main(int argc, const char * argv[])
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	if (argc != 2) {
		cerr << usage << endl;
		exit(EXIT_FAILURE);
	}
	
	
	Nifti image;
	image.open(argv[1], Nifti::READ);
	
	float *data = image.readVolume<float>(0);
	float thresh = 0.;
	
	size_t min_i = -1, min_j = -1, min_k = -1;
	size_t max_i = 0, max_j = 0, max_k = 0;
	for (size_t k = 0; k < image.dim(3); k++) {
		for (size_t j = 0; j < image.dim(2); j++) {
			for (size_t i = 0; i < image.dim(1); i++) {
				if (data[((k * image.dim(2)) + j) * image.dim(1) + i] > thresh) {
					if (i < min_i) min_i = i;
					if (j < min_j) min_j = j;
					if (k < min_k) min_k = k;
					if (i > max_i) max_i = i;
					if (j > max_j) max_j = j;
					if (k > max_k) max_k = k;
				}
			}
		}
	}
					
	cout << min_i << " " << (max_i - min_i) << " "
	     << min_j << " " << (max_j - min_j) << " "
		 << min_k << " " << (max_k - min_k) << endl;
	
    return EXIT_SUCCESS;
}



