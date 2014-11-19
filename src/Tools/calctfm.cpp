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
#include <Eigen/Dense>

#include <getopt.h>

using namespace Eigen;
using namespace std;

const string usage = "calctfm - A utility for calculating transforms that I needed.\n\
\n\
Usage: calctfm [options] X Y Z\n\
Calculates a transform to align the CoG to center and rotate it the right way.\n\
\n\
Options:\n\
	--tfm, -t : Output an Insight Transform file for ANTs\n\
	--mat, -m : Output a .mat file FSL (default)\n\
";

static struct option long_opts[] =
{
	{"tfm", no_argument, 0, 't'},
	{"mat", no_argument, 0, 'm'},
	{0, 0, 0, 0}
};
static const char *short_opts = "tmh";
enum class Format { FSL, ANTs };
Format output = Format::FSL;

int main(int argc, char **argv) {
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, short_opts, long_opts, &indexptr)) != -1) {
		switch (c) {
		case 't': output = Format::ANTs; break;
		case 'm': output = Format::FSL; break;
		case '?': // getopt will print an error message
		case 'h':
			cout << usage << endl;
			return EXIT_FAILURE;
		default:
			break;
		}
	}

	if ((argc - optind) != 3) {
		cerr << "Must have 3 arguments (X, Y, Z)" << endl;
		return EXIT_FAILURE;
	}

	float x = atof(argv[optind++]);
	float y = atof(argv[optind++]);
	float z = atof(argv[optind++]);

	Vector3f CoG{x, y, z};
	Vector3f tgt{0, 1, 0};

	float angle = acos(CoG.dot(tgt) / (CoG.norm() * tgt.norm()));
	Affine3f transform;
	transform = Translation3f(-CoG) * AngleAxisf(angle, Vector3f::UnitZ());

	IOFormat fmt(StreamPrecision, DontAlignCols);
	switch (output) {
	case (Format::FSL):
		cout << transform.matrix() << endl;
		break;
	case (Format::ANTs):
		cout << "#Insight Transform File V1.0" << endl;
		cout << "# Transform 0" << endl;
		cout << "Transform: MatrixOffsetTransformBase_double_3_3" << endl;
		cout << "Parameters: " << transform.matrix().block(0, 0, 1, 3).format(fmt) << " "
		                       << transform.matrix().block(1, 0, 1, 3).format(fmt) << " "
		                       << transform.matrix().block(2, 0, 1, 3).format(fmt) << " "
		                       << transform.matrix().block(0, 3, 3, 1).transpose().format(fmt) << endl;
		cout << "FixedParameters: 0 0 0" << endl;
		break;
	}
	return EXIT_SUCCESS;
}
