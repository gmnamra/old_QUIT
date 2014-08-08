/*
 *  niiext.cpp
 *
 *  Created by Tobias Wood on 19/05/2014.
 *  Copyright (c) 2014 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include <memory>
#include <algorithm>
#include <iostream>
#include <getopt.h>

#include "Nifti/Nifti.h"

using namespace std;
using namespace Nifti;

const string usage = "niiext - A utility for manipulating Nifti extensions.\n\
\n\
Usage: niiext [options] nifti_file\n\
By default any extensions will be summarised and printed to stdout.\n\
\n\
Options:\n\
	--del, -d [N]    : Delete extension number N from file, or all if no\n\
	                   argument.\n\
	--sum, -s [N]    : Summarise (type, length) extension N on stdout, or\n\
	                   all if no argument.\n\
	--ext, -e [N]    : Extract raw bytes of extension N to stdout, or all\n\
	                   if no argument.\n\
	-h, --help:   Print this message and quit.\n";

static struct option long_options[] =
{
	{"add", optional_argument, 0, 'a'},
	{"del", optional_argument, 0, 'd'},
	{"ext", optional_argument, 0, 'e'},
	{"sum", optional_argument, 0, 's'},
	{0, 0, 0, 0}
};

int main(int argc, char **argv) {
	int indexptr = 0, c;

	// Scan once to get the non-optional nifti file name
	while (getopt_long(argc, argv, "d?e?s?h", long_options, &indexptr) != -1) {}
	if ((argc - optind) != 1) {
		cerr << "Must have exactly one Nifti filename to process." << endl << usage << endl;
		return EXIT_FAILURE;
	}
	string fileName(argv[optind]);
	File input(fileName);
	auto exts = input.extensions();
	optind = 1;
	while ((c = getopt_long(argc, argv, "d?e?s?h", long_options, &indexptr)) != -1) {
		switch (c) {
			//case 'a':
			case 'd': {
				if (optarg) {
					int N = atoi(optarg);
					if (N >= 0) {
						auto el = exts.begin();
						advance(el, N);
						exts.erase(el);
					}
				} else {
					exts.clear();
				}
			}	break;
			case 's': {
				int startN = 0, endN = exts.size();
				if (optarg) {
					startN = atoi(optarg);
					endN = 1;
				}
				auto start = exts.begin(); advance(start, startN);
				auto end   = start; advance(end, endN);
				for_each(start, end, [](Extension &e){ cout << "Extension Code: " << e.code() << ", size: " << e.size() << " bytes." << endl; });
			}	break;
			case '?': // getopt will print an error message
			case 'h':
				cout << usage << endl;
				return EXIT_FAILURE;
		}
	}
}
