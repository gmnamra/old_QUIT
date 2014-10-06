/*
 *  niicomplex.cpp
 *
 *  Created by Tobias Wood on 22/04/2014.
 *  Copyright (c) 2014 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <getopt.h>
#include <iostream>

#include "Nifti/Nifti.h"
#include "QUIT/MultiArray.h"
using namespace std;
using namespace Nifti;

//******************************************************************************
// Arguments / Usage
//******************************************************************************
const string usage {
"Usage is: niicomplex [options] inputs outputs\n\
\n\
Default mode is to convert a magnitude/phase image pair into a real/imaginary \n\
image pair. If you have/want different inputs/outputs, then specify the -i/-o \n\
options. The correct number of input/output names must be given as arguments.\n\
\n\
Options:\n\
	--help, -h        : Print this message\n\
	--verbose, -v     : Print more information\n\
	--input, -i m     : Input is magnitude/phase (default)\n\
	            r     : Input is real/imaginary\n\
	            c     : Input is complex\n\
	--output, -o m    : Output will be magnitude/phase images\n\
	             r    : Output will be real/imaginary images (default)\n\
	             c    : Output will be a single complex image\n\
	--dtype, -d f     : Force output datatype to float\n\
	            d     : Force output datatype to double\n\
	            l     : Force output datatype to long double\n"
};

enum class Type { MagPhase, RealImag, Complex };
static bool verbose = false, forceDType = false;
static Type inputType = Type::MagPhase, outputType = Type::RealImag;
static DataType precision = DataType::FLOAT32;
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"input", required_argument, 0, 'i'},
	{"output", required_argument, 0, 'o'},
	{"dtype", required_argument, 0, 'd'},
	{0, 0, 0, 0}
};

template<typename T>
void mag_to_cmp(Nifti::File &in1, Nifti::File &in2, size_t vol,
                vector<T> &v1, vector<T> &v2,
                vector<complex<T>> &c)
{
	if (verbose) cout << "Reading magnitude volume " << vol << endl;
	in1.readVolumes(v1.begin(), v1.end(), vol, 1);
	if (verbose) cout << "Reading phase volume " << vol << endl;
	in2.readVolumes(v2.begin(), v2.end(), vol, 1);
	for (size_t i = 0; i < v1.size(); i++) {
		c[i] = polar(v1[i], v2[i]);
	}
}

template<typename T>
void re_im_to_cmp(Nifti::File &in1, Nifti::File &in2, size_t vol,
                  vector<T> &v1, vector<T> &v2,
                  vector<complex<T>> &c)
{
	if (verbose) cout << "Reading real volume " << vol << endl;
	in1.readVolumes(v1.begin(), v1.end(), vol, 1);
	if (verbose) cout << "Reading imaginary volume " << vol << endl;
	in2.readVolumes(v2.begin(), v2.end(), vol, 1);
	for (size_t i = 0; i < v1.size(); i++) {
		c[i] = complex<T>(v1[i], v2[i]);
	}
}

template<typename T>
void cmp_to_mag(Nifti::File &out1, Nifti::File &out2, size_t vol,
                vector<T> &v1, vector<T> &v2,
                vector<complex<T>> &c)
{
	for (size_t i = 0; i < v1.size(); i++) {
		v1[i] = abs(c[i]); v2[i] = arg(c[i]);
	}
	if (verbose) cout << "Writing magnitude volume " << vol << endl;
	out1.writeVolumes(v1.begin(), v1.end(), vol, 1);
	if (verbose) cout << "Writing phase volume " << vol << endl;
	out2.writeVolumes(v2.begin(), v2.end(), vol, 1);
}

template<typename T>
void cmp_to_re_im(Nifti::File &out1, Nifti::File &out2, size_t vol,
                  vector<T> &v1, vector<T> &v2,
                  vector<complex<T>> &c)
{
	for (size_t i = 0; i < v1.size(); i++) {
		v1[i] = c[i].real(); v2[i] = c[i].imag();
	}
	if (verbose) cout << "Writing real volume " << vol << endl;
	out1.writeVolumes(v1.begin(), v1.end(), vol, 1);
	if (verbose) cout << "Writing imaginary volume " << vol << endl;
	out2.writeVolumes(v2.begin(), v2.end(), vol, 1);
}

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	//**************************************************************************
	// Argument Processing
	//**************************************************************************
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hvi:o:d:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'v': verbose = true; break;
			case 'i':
				switch (*optarg) {
					case 'm': inputType = Type::MagPhase;  cout << "Input is magnitude and phase." << endl; break;
					case 'r': inputType = Type::RealImag; cout << "Input is real and imaginary." << endl; break;
					case 'c': inputType = Type::Complex; cout << "Input is complex." << endl; break;
					default:
						cerr << "Unknown input type " << optarg << endl;
						return EXIT_FAILURE;
				} break;
			case 'o':
				switch (*optarg) {
					case 'm': outputType = Type::MagPhase;  cout << "Output will be magnitude and phase." << endl; break;
					case 'r': outputType = Type::RealImag; cout << "Output will be real and imaginary." << endl; break;
					case 'c': outputType = Type::Complex; cout << "Output will be complex." << endl; break;
					default:
						cerr << "Unknown output type " << optarg << endl;
						return EXIT_FAILURE;
				} break;
			case 'd':
				switch (*optarg) {
					case 'f': precision = DataType::FLOAT32; break;
					case 'd': precision = DataType::FLOAT64; break;
					case 'l': precision = DataType::FLOAT128; break;
					default:
						cerr << "Unknown precision type " << optarg << endl;
						return EXIT_FAILURE;
				} break;
			case 'h':
			case '?': // getopt will print an error message
				cout << usage << endl;
				return EXIT_SUCCESS;
		}
	}

	size_t expected_number_of_arguments = 0;
	switch (inputType) {
		case Type::MagPhase: expected_number_of_arguments += 2; break;
		case Type::RealImag: expected_number_of_arguments += 2; break;
		case Type::Complex:  expected_number_of_arguments += 1; break;
	}
	switch (outputType) {
		case Type::MagPhase: expected_number_of_arguments += 2; break;
		case Type::RealImag: expected_number_of_arguments += 2; break;
		case Type::Complex: expected_number_of_arguments += 1; break;
	}
	if (expected_number_of_arguments != (argc - optind)) {
		cout << "Expected " << expected_number_of_arguments<<  " filenames, but " << (argc - optind) << " were given." << endl << usage << endl;
		return EXIT_FAILURE;
	}

	File in1, in2;
	if (verbose) cout << "Opening input file: " << argv[optind] << endl;
	in1.open(argv[optind++], Mode::Read);
	if (inputType != Type::Complex) {
		if (verbose) cout << "Opening input file: " << argv[optind] << endl;
		in2.open(argv[optind++], Mode::Read);
		if (!in2.header().matchesSpace(in1.header())) {
			cerr << "Input files are incompatible." << endl;
			return EXIT_FAILURE;
		}
	}

	Header outHdr = in1.header();
	if (outputType == Type::Complex) {
		switch (precision) {
			case DataType::FLOAT32: outHdr.setDatatype(DataType::COMPLEX64); break;
			case DataType::FLOAT64: outHdr.setDatatype(DataType::COMPLEX128); break;
			case DataType::FLOAT128: outHdr.setDatatype(DataType::COMPLEX256); break;
			default: throw(std::logic_error("Invalid precision type."));
		}
	} else {
		outHdr.setDatatype(precision);
	}
	File out1, out2;
	if (verbose) cout << "Opening output file: " << argv[optind] << endl;
	out1.setHeader(outHdr);
	out1.open(argv[optind++], Nifti::Mode::Write);
	if (outputType != Type::Complex) {
		if (verbose) cout << "Opening output file: " << argv[optind] << endl;
		out2.setHeader(outHdr);
		out2.open(argv[optind++], Nifti::Mode::Write);
	}

	// Allocate different types to save memory
	vector<complex<float>> cmp_flt;
	vector<complex<double>> cmp_dbl;
	vector<complex<long double>> cmp_ldbl;
	vector<float> flt1, flt2;
	vector<double> dbl1, dbl2;
	vector<long double> ldbl1, ldbl2;

	size_t nEl = in1.matrix().prod();
	switch (precision) {
		case DataType::FLOAT32:  cmp_flt.resize(nEl);  flt1.resize(nEl);  flt2.resize(nEl); break;
		case DataType::FLOAT64:  cmp_dbl.resize(nEl);  dbl1.resize(nEl);  dbl2.resize(nEl); break;
		case DataType::FLOAT128: cmp_ldbl.resize(nEl); ldbl1.resize(nEl); ldbl2.resize(nEl); break;
		default:
			break; // We have checked that this can't happen earlier (famous last words)
	}

	for (size_t vol = 0; vol < in1.dim(4); vol++) {
		if (verbose) cout << "Converting volume " << vol << "..." << endl;
		switch (inputType) {
			case Type::MagPhase: {
				switch (precision) {
					case DataType::FLOAT32:  mag_to_cmp<float>(in1, in2, vol, flt1, flt2, cmp_flt); break;
					case DataType::FLOAT64:  mag_to_cmp<double>(in1, in2, vol, dbl1, dbl2, cmp_dbl); break;
					case DataType::FLOAT128: mag_to_cmp<long double>(in1, in2, vol, ldbl1, ldbl2, cmp_ldbl); break;
					default: break;
				}
			} break;
			case Type::RealImag: {
				switch (precision) {
					case DataType::FLOAT32:  re_im_to_cmp<float>(in1, in2, vol, flt1, flt2, cmp_flt); break;
					case DataType::FLOAT64:  re_im_to_cmp<double>(in1, in2, vol, dbl1, dbl2, cmp_dbl); break;
					case DataType::FLOAT128: re_im_to_cmp<long double>(in1, in2, vol, ldbl1, ldbl2, cmp_ldbl); break;
					default: break;
				}
			} break;
			case Type::Complex : {
				if (verbose) cout << "Reading complex volume " << vol << endl;
				switch (precision) {
					case DataType::FLOAT32: in1.readVolumes(cmp_flt.begin(), cmp_flt.end(), vol, 1); break;
					case DataType::FLOAT64: in1.readVolumes(cmp_dbl.begin(), cmp_dbl.end(), vol, 1); break;
					case DataType::FLOAT128: in1.readVolumes(cmp_ldbl.begin(), cmp_ldbl.end(), vol, 1); break;
					default: break;
				}
			} break;
		}
		switch (outputType) {
			case Type::MagPhase: {
				switch (precision) {
					case DataType::FLOAT32:  cmp_to_mag<float>(out1, out2, vol, flt1, flt2, cmp_flt); break;
					case DataType::FLOAT64:  cmp_to_mag<double>(out1, out2, vol, dbl1, dbl2, cmp_dbl); break;
					case DataType::FLOAT128: cmp_to_mag<long double>(out1, out2, vol, ldbl1, ldbl2, cmp_ldbl); break;
					default: break;
				}
			} break;
			case Type::RealImag: {
				switch (precision) {
					case DataType::FLOAT32:  cmp_to_re_im<float>(out1, out2, vol, flt1, flt2, cmp_flt); break;
					case DataType::FLOAT64:  cmp_to_re_im<double>(out1, out2, vol, dbl1, dbl2, cmp_dbl); break;
					case DataType::FLOAT128: cmp_to_re_im<long double>(out1, out2, vol, ldbl1, ldbl2, cmp_ldbl); break;
					default: break;
				}
			} break;
			case Type::Complex : {
				if (verbose) cout << "Writing complex volume " << vol << endl;
				switch (precision) {
					case DataType::FLOAT32: out1.writeVolumes(cmp_flt.begin(), cmp_flt.end(), vol, 1); break;
					case DataType::FLOAT64: out1.writeVolumes(cmp_dbl.begin(), cmp_dbl.end(), vol, 1); break;
					case DataType::FLOAT128: out1.writeVolumes(cmp_ldbl.begin(), cmp_ldbl.end(), vol, 1); break;
					default: break;
				}
			} break;
		}
	}
	out1.close();
	if (outputType != Type::Complex) {
		out2.close();
	}

	return EXIT_SUCCESS;
}
