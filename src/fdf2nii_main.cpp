//
//  fdf2nii_main.cpp
//
//  Created by Tobias Wood on 29/08/2013.
//
//

#include <string>
#include <iostream>
#include <getopt.h>

#include "fdf.h"
#include "Nifti.h"

using namespace std;

static bool zip = false;
static double scale = 1.;
static string outPrefix;
static struct option long_options[] =
{
	{"scale", required_argument, 0, 's'},
	{"out", required_argument, 0, 'o'},
	{"zip", no_argument, 0, 'z'},
	{0, 0, 0, 0}
};

const string usage {
"fdf2nii - A utility to convert Agilent fdf files to nifti.\n\
\n\
Usage: fdf2nii [opts] image1 image2 ... imageN\n\
image1 to imageN are paths to the Agilent .img folders, not individual .fdf\n\
files\n\
Options:\n\
 -s, --scale:  Scale factor for image dimensions (set to 10 for use with SPM).\n\
 -o, --out:    Specify an output prefix.\n\
 -z, --zip:    Create .nii.gz files\n"
};

int main(int argc, char **argv) {
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "s:o:z", long_options, &indexptr)) != -1) {
		switch (c) {
			case 0: break; // It was an option that just sets a flag.
			case 's': scale = atof(optarg); break;
			case 'o': outPrefix = string(optarg); break;
			case 'z': zip = true; break;
			default: cout << "Unknown option " << optarg << endl;
		}
	}
	
	if ((argc - optind) <= 0) {
		cout << "No input images specified." << endl << usage << endl;
		exit(EXIT_FAILURE);
	}
	while (optind < argc) {
		string path(argv[optind]);
		string outPath = path.substr(0, path.find_last_of(".")) + ".nii";
		if (zip)
			outPath += ".gz";
		cout << "Converting " << path << " to " << outPath << endl;
		Agilent::fdfImage input(path);
		Nifti::File output(input.dim(0), input.dim(1), input.dim(2), input.dim(3),
						   input.voxdim(0), input.voxdim(1), input.voxdim(2), 1.,
			               DT_FLOAT32, input.ijk_to_xyz().cast<float>());
		output.open(outPath, Nifti::Modes::Write);
		for (size_t v = 0; v < input.dim(3); v++) {
			output.writeVolume<float>(v, input.readVolume<float>(v));
		}
		output.close();
		cout << input.ijk_to_xyz() << endl << output.ijk_to_xyz() << endl;
		optind++;
	}
    return EXIT_SUCCESS;
}
