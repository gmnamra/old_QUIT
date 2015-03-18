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
	--value, -v D   : Fill dimension D with constant value v (val)\n\
	--slab, -l D    : Fill dimension D with a slab (val, start, end)\n\
	--grad, -g D    : Fill dimension D with a smooth gradient (low, high)\n\
	--step, -s D    : Fill dimension D with stepped data (low, high, steps)\n\
	--uniform, -U D : Fill dimension D with uniform noise (mid, width)\n\
	--gauss, -G D   : Fill dimension D with gaussian noise (mean, std dev)\n"
};

enum class FillType { Zero, Value, Slab, Gradient, Steps, Uniform, Gaussian };
static bool verbose = false, isBlank = false;
static int ndims = 3;
static vector<FillType> fillTypes(4, FillType::Zero);
static Eigen::Array4f startVal = Eigen::Array4f::Zero(), deltaVal = Eigen::Array4f::Zero();
static Eigen::Array4i stepLength = Eigen::Array4i::Ones();
static vector<uniform_real_distribution<float>> uniforms(4);
static vector<normal_distribution<float>> gauss(4);
static Nifti::DataType dType = Nifti::DataType::FLOAT64;
static Eigen::Affine3f xform = Eigen::Affine3f::Identity();
static struct option long_opts[] = {
	{"help",    no_argument,       0, 'h'},
	{"precision", required_argument, 0, 'p'},
	{"xform",   required_argument, 0, 'x'},
	{"rank",    required_argument, 0, 'r'},
	{"dims",    required_argument, 0, 'd'},
	{"voxdims", required_argument, 0, 'v'},
	{"tr",      required_argument, 0, 't'},
	{"blank",   no_argument,       0, 'b'},
	{"zero",    no_argument,       0, 'z'},
	/*{"value",   required_argument, 0, 'v'},*/
	{"slab",    required_argument, 0, 'l'},
	{"grad",    required_argument, 0, 'g'},
	{"step",    required_argument, 0, 's'},
	{"uniform", required_argument, 0,'U'},
	{"gauss",   required_argument, 0, 'G'},
	{0, 0, 0, 0}
};
static const char* short_opts = "p:r:d:t:x:bzv:l:g:s:U:G:h";
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	size_t expected_extra_args = 1;
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
				stringstream valstream(vals);
				Array<unsigned long, Dynamic, 1> test(ndims);
				QUIT::Read<Array<unsigned long, Dynamic, 1>>::FromLine(valstream, test);
				dims.head(ndims) = test;
			} break;
			case 'v': {
				string vals(optarg);
				stringstream valstream(vals);
				QUIT::Read<ArrayXf>::FromLine(valstream, voxdims.head(3));
			} break;
			case 'b': isBlank = true; break;
			case 't': voxdims[3] = atof(optarg); break;
			case 'l': fillTypes.at(atoi(optarg)) = FillType::Slab; expected_extra_args += 3; break;
			case 'g': fillTypes.at(atoi(optarg)) = FillType::Gradient; expected_extra_args += 2; break;
			case 's': fillTypes.at(atoi(optarg)) = FillType::Steps; expected_extra_args += 3; break;
			case 'U': fillTypes.at(atoi(optarg)) = FillType::Uniform; expected_extra_args += 2; break;
			case 'G': fillTypes.at(atoi(optarg)) = FillType::Gaussian; expected_extra_args += 2; break;
			case 'h':
			case '?': // getopt will print an error message
				return EXIT_FAILURE;
		}
	}
	if ((argc - optind) != expected_extra_args) {
		cerr << "Wrong number of arguments." << endl;
		cout << usage << endl;
		return EXIT_FAILURE;
	}
	string fName(argv[optind++]);

	Nifti::Header hdr(dims, voxdims, dType);
	hdr.setTransform(xform);
	Nifti::File file(hdr, fName);

	if (!isBlank) {
		MultiArray<float, 5> data(dims);
		ind_t starts; starts.setZero();
		ind_t ends = dims;
		for (size_t d = 0; d < ndims; d++) {
			switch (fillTypes.at(d)) {
				case FillType::Zero:
					stepLength[d] = 1;
					break;
				case FillType::Value:
					startVal[d] = atof(argv[optind++]);
					deltaVal[d] = 0;
					stepLength[d] = 1;
					break;
				case FillType::Slab:
					startVal[d] = atof(argv[optind++]);
					deltaVal[d] = 0;
					stepLength[d] = 1;
					starts[d] = atoi(argv[optind++]);
					ends[d] = atoi(argv[optind++]);
					if (starts[d] >= data.dims()[d]) {
						cerr << "Invalid slab start slice." << endl;
						return EXIT_FAILURE;
					}
					if (ends[d] > data.dims()[d]) {
						cerr << "Invalid slab end slice." << endl;
						return EXIT_FAILURE;
					}
					break;
				case FillType::Gradient:
					startVal[d] = atof(argv[optind++]);
					deltaVal[d] = (atof(argv[optind++]) - startVal[d]) / (dims[d] - 1);
					stepLength[d] = 1;
					break;
				case FillType::Steps: {
					startVal[d] = atof(argv[optind++]);
					float endVal = atof(argv[optind++]);
					int steps = atoi(argv[optind++]);
					stepLength[d] = dims[d] / steps;
					deltaVal[d] = (endVal - startVal[d]) / (steps - 1);
				} break;
				case FillType::Uniform: {
					float lo = atof(argv[optind++]);
					float hi = atof(argv[optind++]);
					uniforms[d] = uniform_real_distribution<float>(lo, hi);
				} break;
				case FillType::Gaussian: {
					float mean = atof(argv[optind++]);
					float std  = atof(argv[optind++]);
					gauss[d] = normal_distribution<float>(mean, std);
				} break;
			}
		}
		random_device seed;
		mt19937_64 twist(seed());

		ind_t index; index.setZero();
		function<void (const int&, const float&)> processDim = [&] (const int &d, const float& outVal) {
			float inVal = outVal + startVal[d];
			for (index[d] = starts[d];index[d] < ends[d]; index[d]++) {
				if (fillTypes[d] == FillType::Uniform)
					inVal = outVal + uniforms[d].operator()(twist);
				else if (fillTypes[d] == FillType::Gaussian)
					inVal = outVal + gauss[d].operator()(twist);
				if (d > 0) {
					processDim(d - 1, inVal);
				} else {
					data[index] = inVal;
				}
				if ((index[d] % stepLength[d]) == (stepLength[d] - 1)) inVal += deltaVal[d];
			}
		};
		processDim(ndims - 1, 0.);
		file.writeAll(data.begin(), data.end());
	}
	file.close();
	return EXIT_SUCCESS;
}

