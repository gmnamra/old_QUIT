/*
 *  niinudge.cpp
 *
 *  Copyright (c) 2014 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <getopt.h>

#include <Eigen/Core>
#include "Nifti/Nifti.h"
#include "QUIT/MultiArray.h"
#include "QUIT/Util.h"

using namespace std;
using namespace Eigen;
using namespace QUIT;

const string usage = "niinudge - A utility for moving Nifti images in physical space.\n\
\n\
Usage: niinudge [options] infile outfile\n\
By default nothing happens. Specify one of the options to move your image.\n\
Many of the options require a 3 dimensional vector argument. Valid formats for\n\
this are:\n\
X Y Z - Make sure you encase this format in quotes (\" \")!\n\
\n\
Options:\n\
	--nudge, -n \"X Y Z\"    : Nudge the image (X Y Z added to current offset)\n\
	--offset, -f \"X Y Z\"   : Set the offset to (X,Y,Z)\n\
	--cog, -c                : Move the Center of Gravity to the origin of the\n\
	                           first image, and make subsequent images match\n\
	--rotate, -r \"X/Y/Z a\" : Rotate about specifed axis by specified angle\n\
	--verbose, -v            : Print out what the program is doing\n\
	-h, --help               :   Print this message and quit.\n\
";

static const struct option long_opts[] = {
	{"nudge",  required_argument, 0, 'n'},
	{"offset", required_argument, 0, 'f'},
	{"cog",    no_argument, 0, 'c'},
	{"rotate", required_argument, 0, 'r'},
	{"verbose", no_argument, 0, 'v'},
	{"help",   no_argument, 0, 'h'},
	{0, 0, 0, 0}
};
static const char *short_opts = "n:o:f:cr:vh";
static string prefix;
static int verbose = false, output_transform = false;

Vector3f calc_cog(Nifti::File &f);
Vector3f calc_cog(Nifti::File &f) {
	const auto dims = f.matrix();
	MultiArray<float, 3> a(dims);
	f.open(f.imagePath(), Nifti::Mode::Read);
	f.readVolumes(a.begin(), a.end(), 0, 1);

	Vector3f cog = Vector3f::Zero();
	float mass = 0.;
	for (size_t k = 0; k < dims[2]; k++) {
		for (size_t j = 0; j < dims[1]; j++) {
			for (size_t i = 0; i < dims[0]; i++) {
				float val = a[{i,j,k}];
				if (isfinite(val)) {
					cog += val * Vector3f(i,j,k);
					mass += val;
				}
			}
		}
	}
	cog /= mass;
	if (verbose) cout << "CoG in voxels: " << cog.transpose() << ", mass: " << mass << endl;
	cog = f.header().transform() * cog;
	if (verbose) cout << "CoG in space:  " << cog.transpose() << endl;
	f.close();
	return cog;
}

void write_transform(const Affine3f &tfm, const string path);
void write_transform(const Affine3f &tfm, const string path){
	ofstream file(path);

	IOFormat fmt(StreamPrecision, DontAlignCols);
	file << "#Insight Transform File V1.0" << endl;
	file << "# Transform 0" << endl;
	file << "Transform: MatrixOffsetTransformBase_double_3_3" << endl;
	file << "Parameters: " << tfm.matrix().block<1,3>(0,0).format(fmt)
	     << " " << tfm.matrix().block<1,3>(1,0).format(fmt)
	     << " " << tfm.matrix().block<1,3>(2,0).format(fmt)
	     << " " << tfm.translation().matrix().transpose().format(fmt) << endl;
	file << "FixedParameters: 0 0 0" << endl;
	file.close();
}

int main(int argc, char **argv) {
	// Make a first pass to permute the options and get filenames at the end
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
		case 'v': verbose = true; break;
		case '?': // getopt will print an error message
		case 'h':
			cout << usage << endl;
			return EXIT_FAILURE;
		default:
			break;
		}
	}
	if ((argc - optind) <= 0 ) {
		cerr << "No input image file specified." << endl;
		cout << usage << endl;
		return EXIT_FAILURE;
	}
	string inFilename(argv[optind++]);
	string outFilename(argv[optind++]);
	if (verbose) cout << "Opening input file: " << inFilename << endl;
	Nifti::File inFile(inFilename);
	// Now reset optind and process options
	optind = 1;
	Nifti::Header hdr = inFile.header();
	Affine3f xfm = hdr.transform();
	if (verbose) cout << "Initial transform: " << endl << xfm.matrix() << endl;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
		case 'n': {
			Vector3f nudge;
			ReadEigen(optarg, nudge);
			if (verbose) cout << "Nudging by: " << nudge.transpose() << endl;
			xfm = Translation3f(nudge) * xfm;
		} break;
		case 'f': {
			Vector3f nudge;
			ReadEigen(optarg, nudge);
			if (verbose) cout << "Setting offset to: " << nudge.transpose() << endl;
			xfm.translation() = nudge;
		} break;
		case 'c': {
			if (verbose) cout << "Aligning origin to CoG." << endl;
			Vector3f CoG = calc_cog(inFile);
			xfm = Translation3f(-CoG) * xfm;
		} break;
		case 'r': {
			stringstream args{string{optarg}};
			string axis; args >> axis;
			float angle; args >> angle;
			if (verbose) cout << "Rotating image by " << angle << " around " << axis << " axis." << endl;
			if (axis == "X") {
				xfm = AngleAxisf(angle * M_PI / 180., Vector3f::UnitX()) * xfm;
			} else if (axis == "Y") {
				xfm = AngleAxisf(angle * M_PI / 180., Vector3f::UnitY()) * xfm;
			} else if (axis == "Z") {
				xfm = AngleAxisf(angle * M_PI / 180., Vector3f::UnitZ()) * xfm;
			} else {
				throw(runtime_error("Invalid axis specification: " + axis));
			}
		} break;
		}
	}
	if (verbose) cout <<  "Final transform: " << endl << xfm.matrix() << endl;
	vector<char> data(inFile.dataSize());
	inFile.readBytes(data);
	if (verbose) cout << "Writing file: " << outFilename << endl;
	hdr.setTransform(xfm);
	Nifti::File outFile(hdr, outFilename, inFile.extensions());
	outFile.writeBytes(data);
	inFile.close();
	outFile.close();
	return EXIT_SUCCESS;
}


