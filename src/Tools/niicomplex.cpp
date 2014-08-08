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
static DataType outDType = DataType::COMPLEX64;
static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"verbose", no_argument, 0, 'v'},
	{"input", required_argument, 0, 'i'},
	{"output", required_argument, 0, 'o'},
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
						exit(EXIT_FAILURE);
						break;
				} break;
			case 'o':
				switch (*optarg) {
					case 'm': outputType = Type::MagPhase;  cout << "Output will be magnitude and phase." << endl; break;
					case 'r': outputType = Type::RealImag; cout << "Output will be real and imaginary." << endl; break;
					case 'c': outputType = Type::Complex; cout << "Output will be complex." << endl; break;
					default:
						cerr << "Unknown output type " << optarg << endl;
						exit(EXIT_FAILURE);
						break;
				} break;
			case 'h':
			case '?': // getopt will print an error message
				exit(EXIT_FAILURE);
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
		case Type::Complex:  expected_number_of_arguments += 1; break;
	}
	if (expected_number_of_arguments != (argc - optind)) {
		cout << "Expected " << expected_number_of_arguments<<  " filenames, but " << (argc - optind) << " were given." << endl << usage << endl;
		exit(EXIT_FAILURE);
	}

	File file1, file2;
	if (verbose) cout << "Opening input file: " << argv[optind] << endl;
	file1.open(argv[optind++], Mode::Read);
	size_t nEl = file1.dims().prod();
	vector<complex<long double>> complexData(nEl);

	switch (inputType) {
		case Type::MagPhase: {
			if (verbose) cout << "Opening input file: " << argv[optind] << endl;
			file2.open(argv[optind++], Mode::Read);
			if (!file2.header().matchesSpace(file1.header())) {
				cerr << "Magnitude and phase files are incompatible." << endl;
				exit(EXIT_FAILURE);
			}
			vector<long double> magData(nEl), phaseData(nEl);
			file1.readAll(magData.begin(), magData.end());
			file2.readAll(phaseData.begin(), phaseData.end());
			for (size_t i = 0; i < nEl; i++) {
				complexData[i] = polar(magData[i], phaseData[i]);
			}
			file2.close();
		} break;
		case Type::RealImag: {
			if (verbose) cout << "Opening input file: " << argv[optind] << endl;
			file2.open(argv[optind++], Nifti::Mode::Read);
			if (!file2.header().matchesSpace(file1.header())) {
				cerr << "Real and imaginary files are incompatible." << endl;
				exit(EXIT_FAILURE);
			}
			vector<long double> realData(nEl), imagData(nEl);
			file1.readAll(realData.begin(), realData.end());
			file2.readAll(imagData.begin(), imagData.end());
			for (size_t i = 0; i < nEl; i++) {
				complexData[i] = complex<long double>(realData[i], imagData[i]);
			}
			file2.close();
		} break;
		case Type::Complex : {
			file1.readAll(complexData.begin(), complexData.end());
		}
	}

	file1.close();
	Header outHdr = file1.header();
	if (forceDType) {
		outHdr.setDatatype(outDType);
	} else {

	}

	switch (outputType) {
		case Type::MagPhase: {
			switch (outHdr.datatype()) {
				case (Nifti::DataType::FLOAT32) : case (Nifti::DataType::COMPLEX64) :
					outHdr.setDatatype(Nifti::DataType::FLOAT32);  break;
				case (Nifti::DataType::FLOAT64) : case (Nifti::DataType::COMPLEX128) :
					outHdr.setDatatype(Nifti::DataType::FLOAT64); break;
				case (Nifti::DataType::FLOAT128) : case (Nifti::DataType::COMPLEX256) :
					outHdr.setDatatype(Nifti::DataType::FLOAT128); break;
				default: outHdr.setDatatype(Nifti::DataType::FLOAT64); break;
			}
			vector<long double> absData(nEl), argData(nEl);
			for (size_t i = 0; i < nEl; i++) {
				absData[i] = abs(complexData[i]);
				argData[i] = arg(complexData[i]);
			}
			if (verbose) cout << "Writing magnitude file: " << argv[optind] << endl;
			File out1(outHdr, argv[optind++]);
			out1.writeAll(absData.begin(), absData.end());
			out1.close();
			if (verbose) cout << "Writing phase file: " << argv[optind] << endl;
			File out2(outHdr, argv[optind++]);
			out2.writeAll(argData.begin(), argData.end());
			out2.close();
		} break;
		case Type::RealImag: {
			switch (outHdr.datatype()) {
				case (Nifti::DataType::FLOAT32) : case (Nifti::DataType::COMPLEX64) :
					outHdr.setDatatype(Nifti::DataType::FLOAT32);  break;
				case (Nifti::DataType::FLOAT64) : case (Nifti::DataType::COMPLEX128) :
					outHdr.setDatatype(Nifti::DataType::FLOAT64); break;
				case (Nifti::DataType::FLOAT128) : case (Nifti::DataType::COMPLEX256) :
					outHdr.setDatatype(Nifti::DataType::FLOAT128); break;
				default: outHdr.setDatatype(Nifti::DataType::FLOAT128); break;
			}
			vector<long double> realData(nEl), imagData(nEl);
			for (size_t i = 0; i < nEl; i++) {
				realData[i] = real(complexData[i]);
				imagData[i] = imag(complexData[i]);
			}
			if (verbose) cout << "Writing real file: " << argv[optind] << endl;
			File out1(outHdr, argv[optind++]);
			out1.writeAll(realData.begin(), realData.end());
			out1.close();
			if (verbose) cout << "Writing imaginary file: " << argv[optind] << endl;
			File out2(outHdr, argv[optind++]);
			out2.writeAll(imagData.begin(), imagData.end());
			out2.close();

		} break;
		case Type::Complex : {
			switch (outHdr.datatype()) {
				case (Nifti::DataType::FLOAT32) :    outHdr.setDatatype(Nifti::DataType::COMPLEX64);  break;
				case (Nifti::DataType::FLOAT64) :    outHdr.setDatatype(Nifti::DataType::COMPLEX128); break;
				case (Nifti::DataType::FLOAT128) :   outHdr.setDatatype(Nifti::DataType::COMPLEX256); break;
				default: outHdr.setDatatype(Nifti::DataType::COMPLEX64); break;
			}
			if (verbose) cout << "Writing complex file: " << argv[optind] << endl;
			File out(outHdr, argv[optind++]);
			out.writeAll(complexData.begin(), complexData.end());
			out.close();
		} break;
	}
	return EXIT_SUCCESS;
}
