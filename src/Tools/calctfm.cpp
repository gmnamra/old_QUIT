/*
 *  calctfm.cpp
 *
 *  Created by Tobias Wood on 19/11/2014.
 *  Copyright (c) 2014 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <iostream>
#include <fstream>
#include <Eigen/Dense>

#include <getopt.h>

#include "QUIT/Util.h"

using namespace Eigen;
using namespace std;

const string usage = "calctfm - A utility for calculating transforms that I needed.\n\
\n\
Usage: calctfm [options] [--] X Y Z filename\n\
Calculates a transform to align the CoG to center and rotate it the right way.\n\
WARNING - If your CoG has negative values make sure you put -- before it.\n\
\n\
Options:\n\
	--tfm, -t   : Output an Insight Transform file for ANTs\n\
	--mat, -m   : Output a .mat file FSL (default)\n\
	--comp, -p  : Write the \"compass point\" to stdout\n\
	--corax, -c : Data has been rotated to match human coronal definition\n\
	--orig, -o \"X Y Z\" : Specify a different origin for angle calculation\n\
	--in, -i    : Bottom of skull faces in towards origin (default opposite)\n";

static struct option long_opts[] =
{
	{"tfm", no_argument, 0, 't'},
	{"mat", no_argument, 0, 'm'},
	{"comp", no_argument, 0, 'p'},
	{"corax", no_argument, 0, 'c'},
	{"orig", required_argument, 0, 'o'},
	{"in", no_argument, 0, 'i'},
	{0, 0, 0, 0}
};
static const char *short_opts = "tmphvco:i";
enum class Format { FSL, ANTs };
Format output = Format::FSL;
Vector3f origin = Vector3f::Zero();
bool verbose = false, compass = false, corax = false, inwards = false;

int main(int argc, char **argv) {
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
		case 't': output = Format::ANTs; break;
		case 'm': output = Format::FSL; break;
		case 'p': compass = true; break;
		case 'v': verbose = true; break;
		case 'c': corax = true; break;
		case 'o': origin = QUIT::parse_vector(optarg); break;
		case 'i': inwards = true; break;
		case '?': // getopt will print an error message
		case 'h':
			cout << usage << endl;
			return EXIT_FAILURE;
		default:
			break;
		}
	}

	if ((argc - optind) != 4) {
		cerr << "Must have 4 arguments (X, Y, Z, filename)" << endl << usage << endl;
		return EXIT_FAILURE;
	}

	float x = atof(argv[optind++]);
	float y = atof(argv[optind++]);
	float z = atof(argv[optind++]);
	string filename(argv[optind++]);

	Vector3f CoG{x, y, z};
	float angle;
	if (corax) {
		angle = atan2(z - origin(2), x - origin(0));
	} else {
		angle = atan2(y - origin(1), x - origin(0));
	}

	if (compass) {
		int octant = fmod(4 * angle / M_PI + 4.5, 8);
		static const array<string, 8> compass_points{"W","SW","S","SE","E","NE","N","NW"};
		cout << compass_points[octant] << endl;
	}
	// Now we want to rotate the head so that the line to the CoG becomes parallel to (0, 1)
	// or (0, -1) depending on whether the heads were scanned with bottom of skull towards center
	// or not
	if (inwards)
		angle = (M_PI / 2.) - angle;
	else
		angle =  (3. * M_PI / 2.) - angle;
	Affine3f transform;
	transform = Translation3f(-CoG) * AngleAxisf(-angle, Vector3f::UnitZ());
	IOFormat fmt(StreamPrecision, DontAlignCols);
	ofstream file(filename);
	switch (output) {
	case (Format::FSL):
		file << transform.matrix() << endl;
		break;
	case (Format::ANTs):
		file << "#Insight Transform File V1.0" << endl;
		file << "# Transform 0" << endl;
		file << "Transform: MatrixOffsetTransformBase_double_3_3" << endl;
		file << "Parameters: " << transform.matrix().block(0, 0, 1, 3).format(fmt) << " "
		                       << transform.matrix().block(1, 0, 1, 3).format(fmt) << " "
		                       << transform.matrix().block(2, 0, 1, 3).format(fmt) << " "
		                       << transform.matrix().block(0, 3, 3, 1).transpose().format(fmt) << endl;
		file << "FixedParameters: 0 0 0" << endl;
		break;
	}
	return EXIT_SUCCESS;
}
