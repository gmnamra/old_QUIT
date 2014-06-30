//
//  fdfbval.cpp
//
//  Created by Tobias Wood on 19/05/2014.
//
//  Copyright Tobias Wood 2014

#include <string>
#include <iostream>
#include <getopt.h>

#include "QUIT/Util.h"
#include "Nifti/Nifti.h"
#include "Nifti/ExtensionCodes.h"

using namespace std;
using namespace QUIT;
using namespace Nifti;

int main(int argc, char **argv) {

	if (argc != 2) {
		cout << "Exactly one input file (procpar or nifti) required." << endl;
		exit(EXIT_FAILURE);
	}

	string iname(argv[1]), outprefix;
	cout << "Reading diffusion directions from: " << iname << endl;
	Agilent::ProcPar pp;
	try {
		Nifti1 f(iname, Mode::Read);
		ReadPP(f, pp);
		outprefix = f.basePath();
	} catch (runtime_error &e) {
		ifstream f(iname);
		f >> pp;
		outprefix = iname;
	}

	auto image = pp.realValues("image");
	auto bval  = pp.realValues("bvalue");
	auto dro   = pp.realValues("dro");
	auto dpe   = pp.realValues("dpe");
	auto dsl   = pp.realValues("dsl");
	size_t nvol = (image == 1).cast<int>().sum();
	Eigen::ArrayXf outval(nvol), outro(nvol), outpe(nvol), outsl(nvol);
	size_t j = 0;
	for (size_t i = 0; i < image.rows(); i++) {
		if (image[i] == 1) {
			outval[j] = bval[i];
			outro[j]  = dro[i];
			outpe[j]  = dpe[i];
			outsl[j]  = dsl[i];
			j++;
		}
	}

	cout << "Writing diffusion values to:     " << outprefix+".bval" << endl;
	cout << "Writing diffusion directions to: " << outprefix+".bvec" << endl;
	ofstream ovals(outprefix + ".bval");
	ofstream ovecs(outprefix + ".bvec");

	ovals << outval.transpose() << endl;
	ovecs << outro.transpose() << endl;
	ovecs << outpe.transpose() << endl;
	ovecs << outsl.transpose() << endl;

	ovals.close();
	ovecs.close();
	cout << "Finished." << endl;
}
