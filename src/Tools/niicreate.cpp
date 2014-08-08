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
using namespace QUIT;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: niicreate [options] filename dims voxdims [extra arguments]\n\
\n\
This is a tool to create Nifti files, either blank headers with orientation\n\
information, e.g. for registration, or files filled with simple patterns of\n\
data e.g. solid values, gradients, or blocks. The default is to create a\n\
3D file filled with zeros. Choose from the options below to create something\n\
else.\n\
\n\
Main Options:\n\
	--help, -h       : Print this message\n\
	--dtype, -t F    : Set the datatype to 64 bit float (default)\n\
	            I    : Set the datatype to 16 bit int\n\
	            C    : Set the datatype to 128 bit complex\n\
	            NNN  : Set the datatype to the given valid Nifti datatype\n\
	--xform, -x FILE : Copy header transform from another Nifti\n\
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
static vector<FillType> fillTypes(4, FillType::Zero);
static Eigen::Array4f startVal = Eigen::Array4f::Zero(), deltaVal = Eigen::Array4f::Zero();
static Eigen::Array4i stepLength = Eigen::Array4i::Ones();
static vector<uniform_real_distribution<float>> uniforms(4);
static vector<normal_distribution<float>> gauss(4);
static Nifti::DataType dType = Nifti::DataType::FLOAT64;
static Eigen::Affine3f xform = Eigen::Affine3f::Identity();
static struct option long_options[] = {
	{"help", no_argument, 0, 'h'},
	{"dtype", required_argument, 0, 't'},
	{"xform", required_argument, 0, 'x'},
	{"blank", no_argument, 0, 'b'},
	{"zero", no_argument, 0, 'z'},
	{"value", required_argument, 0, 'v'},
	{"slab", required_argument, 0, 'l'},
	{"grad", required_argument, 0, 'g'},
	{"step", required_argument, 0, 's'},
	{"uniform", required_argument, 0,'U'},
	{"gauss", required_argument, 0, 'G'},
	{0, 0, 0, 0}
};
//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	size_t expected_extra_args = 9;
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "d:t:x:bzv:l:g:s:U:G:h", long_options, &indexptr)) != -1) {
		switch (c) {
			case 't':
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
			case 'b': isBlank = true; break;
			case 'v': fillTypes.at(atoi(optarg)) = FillType::Value; expected_extra_args += 1; break;
			case 'l': fillTypes.at(atoi(optarg)) = FillType::Slab; expected_extra_args += 3; break;
			case 'g': fillTypes.at(atoi(optarg)) = FillType::Gradient; expected_extra_args += 2; break;
			case 's': fillTypes.at(atoi(optarg)) = FillType::Steps; expected_extra_args += 3; break;
			case 'U': fillTypes.at(atoi(optarg)) = FillType::Uniform; expected_extra_args += 2; break;
			case 'G': fillTypes.at(atoi(optarg)) = FillType::Gaussian; expected_extra_args += 2; break;
			case 'h':
			case '?': // getopt will print an error message
				exit(EXIT_FAILURE);
		}
	}
	if ((argc - optind) != expected_extra_args) {
		cerr << "Wrong number of arguments." << endl;
		cout << usage << endl;
		exit(EXIT_FAILURE);
	}

	string fName(argv[optind++]);
	fName += OutExt();
	MultiArray<float, 4>::Index dims;
	Eigen::ArrayXf vdims(4);
	for (size_t i = 0; i < 4; i++) { dims[i] = atoi(argv[optind++]); }
	for (size_t i = 0; i < 4; i++) { vdims[i] = atof(argv[optind++]); }
	Nifti::Header hdr(dims, vdims, dType);
	hdr.setTransform(xform);
	Nifti::File file(hdr, fName);

	if (!isBlank) {
		MultiArray<float, 4> data(dims);
		MultiArray<float, 4>::Index starts; starts.setZero();
		MultiArray<float, 4>::Index ends = data.dims();
		for (size_t d = 0; d < 4; d++) {
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
						exit(EXIT_FAILURE);
					}
					if (ends[d] > data.dims()[d]) {
						cerr << "Invalid slab end slice." << endl;
						exit(EXIT_FAILURE);
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

		MultiArray<float, 4>::Index index;
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
		processDim(3, 0.);
		file.writeAll(data.begin(), data.end());
	}
	file.close();
	return EXIT_SUCCESS;
}

