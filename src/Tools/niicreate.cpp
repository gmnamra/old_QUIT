/*
 *  niicreate.cpp
 *
 *  Created by Tobias Wood on 30/05/2014.
 *  Copyright (c) 2014 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <getopt.h>
#include <iostream>
#include <random>
#include <functional>

#include "Nifti/Nifti.h"
#include "QUIT/QUIT.h"
using namespace std;
using namespace Eigen;
using namespace QUIT;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: niicreate filename dims voxdims [options]\n\
\n\
This is a tool to create Nifti files, either blank headers with orientation\n\
information, e.g. for registration, or files filled with simple patterns of\n\
data e.g. solid values, gradients, or blocks. The default is to create a\n\
3D file filled with zeros. Choose from the options below to create something\n\
else.\n\
\n\
Main Options:\n\
	--help, -h       : Print this message\n\
	--precision, -p F    : Set the datatype to 32 bit float (default)\n\
	                I    : Set the datatype to 16 bit int\n\
	                C    : Set the datatype to 64 bit complex\n\
	                NNN  : Set the datatype to the given valid Nifti datatype\n\
	--xform, -x FILE : Copy header transform from another Nifti\n\
	--rank, -r N     : Set number of dimensions (max 5, default 3)\n\
	--dims, -d \"N N N\"    : Set the dimensions (default 16 16 16)\n\
	--voxdims, -v \"X X X\" : Set the voxel dimensions (default 1mm iso)\n\
	--tr,  -t X             : Set the TR (default 1s)\n\
File Content Options:\n\
	--blank, -b     : Create the header information only\n\
	--fill, -f X    : Fill the entire image with value X (default 0)\n\
	--grad, -g \"D L H\"    : Fill dimension D with a smooth gradient (low, high)\n\
	--step, -s \"D L H S\"  : Fill dimension D with stepped data (low, high, steps)\n"
};

enum class FillType { Fill, Gradient, Steps };
static bool verbose = false, isBlank = false;
static int ndims = 3;
static FillType fillType = FillType::Fill;
static float startVal = 0, deltaVal = 0;
static int stepLength = 1, fillDim = 0;
static Nifti::DataType dType = Nifti::DataType::FLOAT64;
static Eigen::Affine3f xform = Eigen::Affine3f::Identity();
static struct option long_opts[] = {
	{"help",      no_argument,       0, 'h'},
	{"precision", required_argument, 0, 'p'},
	{"xform",     required_argument, 0, 'x'},
	{"rank",      required_argument, 0, 'r'},
	{"dims",      required_argument, 0, 'd'},
	{"voxdims",   required_argument, 0, 'v'},
	{"tr",        required_argument, 0, 't'},
	{"blank",     no_argument,       0, 'b'},
	{"fill",      required_argument, 0, 'f'},
	{"grad",      required_argument, 0, 'g'},
	{"step",      required_argument, 0, 's'},
	{0, 0, 0, 0}
};
static const char* short_opts = "hp:x:r:d:v:t:bf:g:s:";
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	typedef MultiArray<float, 5>::Index ind_t;
	ind_t dims; dims.setOnes();
	ArrayXf voxdims(5); voxdims.setOnes();

	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
			case 'p':
				switch (*optarg) {
					case 'F': dType = Nifti::DataType::FLOAT32; break;
					case 'I': dType = Nifti::DataType::INT16; break;
					case 'C': dType = Nifti::DataType::COMPLEX64; break;
					default: dType = Nifti::DataTypeForCode(atoi(optarg)); break;
				} break;
			case 'x': {
				Nifti::File other(optarg);
				xform = other.header().transform();
			} break;
			case 'r':
				ndims = atoi(optarg);
				if ((ndims < 2) || (ndims > 5)) {
					cerr << "Invalid number of dimensions. Must be 2-5." << endl;
					return EXIT_FAILURE;
				}
				break;
			case 'd': {
				string vals(optarg);
				Array<unsigned long, Dynamic, 1> test(ndims);
				QUIT::ReadEigenFromString(vals, test);
				dims.head(ndims) = test;
			} break;
			case 'v': {
				string vals(optarg);
				QUIT::ReadEigenFromString(vals, voxdims.head(3));
			} break;
			case 't': voxdims[3] = atof(optarg); break;
			case 'b': isBlank = true; break;
			case 'f': fillType = FillType::Fill; startVal = atof(optarg); break;
			case 'g': {
				fillType = FillType::Gradient;
				ArrayXf temp(3);
				QUIT::ReadEigenFromString(string(optarg), temp);
				fillDim = temp[0];
				startVal = temp[1];
				deltaVal = (temp[2] - startVal) / (dims[fillDim] - 1);
				stepLength = 1;
			} break;
			case 's': {
				fillType = FillType::Steps;
				ArrayXf temp(4);
				QUIT::ReadEigenFromString(string(optarg), temp);
				fillDim = temp[0];
				startVal = temp[1];
				float endVal = temp[2];
				int steps = temp[3];
				stepLength = dims[fillDim] / steps;
				deltaVal = (endVal - startVal) / (steps - 1);
			} break;
			case 'h':
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) < 1) {
		cerr << "Missing input filename." << endl;
		cout << usage << endl;
		return EXIT_FAILURE;
	} else if ((argc - optind) > 1) {
		cerr << "Unexpected extra arguments." << endl;
		cout << usage << endl;
		return EXIT_FAILURE;
	}
	string fName(argv[optind++]);

	Nifti::Header hdr(dims, voxdims, dType);
	hdr.setTransform(xform);
	Nifti::File file(hdr, fName);

	if (!isBlank) {
		MultiArray<float, 5> data(dims);
		ind_t start = ind_t::Zero();
		ind_t size = dims;
		size[fillDim] = 0;

		float val = startVal;
		for (int i = 0; i < dims[fillDim]; i++) {
			start[fillDim] = i;
			auto slice = data.slice<4>(start,size);
			for (auto v = slice.begin(); v != slice.end(); ++v) {
				*v = val;
			}
			if ((i % stepLength) == (stepLength - 1)) val += deltaVal;
		}
		file.writeAll(data.begin(), data.end());
	}
	file.close();
	return EXIT_SUCCESS;
}

