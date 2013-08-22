//
//  main.c
//  procparse
//
//  Created by Tobias Wood on 10/07/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#include <iostream>
#include <getopt.h>
#include "procpar.h"

using namespace std;
using namespace Recon;

static int full = false, partial = false, verbose = false;

const string usage = "procparse - A utility to find interesting information in Agilent procpar files.\n\
\n\
Usage: procparse [opts] file1 par1 par2 ... parN\n\
par1 to parN are parameter names to search for in procpar. If none are specified \
then the whole file will be listed.\n\
Options:\n\
 -f, --full:       Print the full parameter information, not a shortened version.\n\
 -p, --partial:    Print parameters that are partial matches.\n\
 -c, --cmp file2:  Compare to file 2. List parameters that differ or are missing.\n\
 -v, --verbose:    Print more information.\n";

int main(int argc, char **argv)
{
	static struct option long_options[] =
	{
		{"full", no_argument, &full, true},
		{"partial", no_argument, &partial, true},
		{"cmp",  required_argument, 0, 'c'},
		{"verbose", no_argument, &verbose, true},
		{0, 0, 0, 0}
	};
	
	int indexptr = 0, c;
	string procparFile = "", cmpFile = "";
	while ((c = getopt_long(argc, argv, "fphvc:", long_options, &indexptr)) != -1) {
		switch (c) {
			case 0: break; // It was an option that just sets a flag.
			case 'f': full = true; break;
			case 'p': partial = true; break;
			case 'h': cout << usage << endl; break;
			case 'v': verbose = true; break;
			case 'c': cmpFile = string(optarg); break;
			default: cout << "Unknown option " << optarg << endl;
		}
	}
	
	ParameterList pars, cmpPars;
	ParameterList::iterator par, cmpPar;
	
	if ((argc - optind) <= 0) {
		cout << "No procpar file specified." << endl << usage << endl;
		exit(EXIT_FAILURE);
	}
	
	procparFile = string(argv[optind++]);
	if (verbose)
		cout << "Reading procpar file: " << procparFile << endl;
	
	if (!ReadProcpar(procparFile, pars)) {
		cout << "Could not find procpar file.";
		exit(EXIT_FAILURE);
	}
	if (cmpFile != "") {
		if (verbose)
			cout << "Reading comparison procpar file: " << cmpFile;
		if (!ReadProcpar(cmpFile, cmpPars)) {
			cout << "Could not find procpar file.";
			exit(EXIT_FAILURE);
		}
	}
	
	if (optind == argc && cmpFile != "") {
		// Look for parameters present in one file but not the other
		for (par = pars.begin(); par != pars.end(); par++) {
			cmpPar = cmpPars.find(par->first);
			if (cmpPar == cmpPars.end()) {
				cout << "Parameter " << par->first << " missing in " << cmpFile << endl;
			} else {
				if (par->second != cmpPar->second) {
					cout << "Parameter " << par->first << ": " << par->second.print_values() << " vs " << cmpPar->second.print_values() << endl;
				}
			}
		}
		for (cmpPar = cmpPars.begin(); cmpPar != cmpPars.end(); cmpPar++) {
			par = pars.find(cmpPar->first);
			if (par == pars.end())
				cout << "Parameter " << cmpPar->first << " missing in " << procparFile << endl;
		}
	} else {
		while (optind < argc) {
			if (partial) {
				if (verbose)
					cout << "Listing partial matches for: " << argv[optind] << endl;
				int matches = 0;
				for (par = pars.begin(); par != pars.end(); par++) {
					if (par->first.find(string(argv[optind])) != string::npos) {
						if (full)
							cout << par->second << endl;
						else
							cout << par->first << ": " << par->second.print_values() << endl;
						matches++;
						
						if (cmpFile != "") {
							cmpPar = cmpPars.find(par->first);
							if (cmpPar == cmpPars.end())
								cout << "Not present in " << cmpFile << endl;
							else if (!(cmpPar->second == par->second)) {
								cout << "Differs in " << cmpFile << endl;
								if (full)
									cout << cmpPar->second << endl;
								else
									cout << cmpPar->first << ": " << cmpPar->second.print_values() << endl;
							}
						}
					}
				}
				if (verbose)
					cout << "Found " << matches << " matches." << endl;
			} else {
				string search(argv[optind]);
				par = pars.find(search);
				if (par != pars.end()) {
					if (full)
						cout << par->second << endl;
					else
						cout << par->first << ": " << par->second.print_values() << endl;
				}
				if (cmpFile != "") {
					cmpPar = cmpPars.find(search);
					if (cmpPar != cmpPars.end()) {
						cout << "In " << cmpFile << endl;
						if (full)
							cout << cmpPar->second << endl;
						else
							cout << cmpPar->first << ": " << cmpPar->second.print_values() << endl;
					}
				}
			}
			optind++;
		}
	}
    return EXIT_SUCCESS;
}

