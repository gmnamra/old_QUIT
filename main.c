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
using namespace ProcPar;

static int fullPrint = false, all = false, list = false, verbose = false;

void listMatchingPars(const ParameterList &pars, const vector<string> &patterns);
void listMatchingPars(const ParameterList &pars, const vector<string> &patterns)
{
	for (int i = 0; i < patterns.size(); i++)
	{
		if (verbose)
			cout << "Listing parameter names containing: " << patterns[i] << endl;
		
		ParameterList::const_iterator it;
		int matches = 0;
		for (it = pars.begin(); it != pars.end(); it++) {
			if (it->first.find(patterns[i]) != string::npos) {
				cout << it->first;
				matches++;
			}
		}
		if (verbose)
			cout << "Found " << matches << " matches." << endl;
	}
}

const string usage = "procparse - A utility to find interesting information in Agilent procpar files.\n\
\n\
Usage: procparse [opts] file1 par1 par2 ... parN\n\
par1 to parN are parameter names to search for in procpar. If none are specified \
they can be entered at via stdin.\n\
Options:\n\
 --full:     Print the full parameter information, not a shortened version.\n\
 --all:      Print all parameters in file.\n\
 --list:     Print parameter names that contain any of par1 ... parN.\n\
 --cmp file2:	Open file2 and list only parameters that differ.\n";

int main(int argc, char **argv)
{
	static struct option long_options[] =
	{
		{"full", no_argument, &fullPrint, true},
		{"all",  no_argument, &all, true},
		{"list", no_argument, &list, true},
		{"cmp",  required_argument, 0, 'c'},
		{0, 0, 0, 0}
	};
	
	int indexptr = 0, c;
	string procparFile = "", cmpFile = "";
	while ((c = getopt_long(argc, argv, "hv", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 0:
				// It was an option that just sets a flag.
				break;
			case 'c':
				cmpFile = string(optarg);
				break;
			case 'v':
				verbose = true;
				break;
			case 'h':
				cout << usage << endl;
				break;
			default:
				cout << "Unknown option " << optarg << endl;
		}
	}
	
	ParameterList pars, cmpPars;
	
	procparFile = string(argv[optind++]);
	if (verbose)
		cout << "Reading procpar file: " << procparFile << endl;
	pars = ReadProcpar(procparFile);
	if (cmpFile != "") {
		if (verbose)
			cout << "Reading comparison procpar file: " << cmpFile;
		cmpPars = ReadProcpar(cmpFile);
	}
	
	Parameter par, cmpPar, tmpPar;
	
	if (optind < argc) {
		if (list) {
			vector<string> patterns(argc - optind);
			for (int i = 0; i < argc - optind; i++)
				patterns[i] = string(argv[optind + i]);
			listMatchingPars(pars, patterns);
		} else {
			if (verbose)
				cout << "Searching for " << (argc - optind) << " parameters." << endl;
			for (int i = optind; i < argc; i++) {
				string search(argv[i]);
				ParameterList::iterator it = pars.find(search);
				if ((it != pars.end()) && (cmpFile != "")) {
					ParameterList::iterator cmpit = cmpPars.find(search);
					if (cmpit != cmpPars.end()) {
						if (it->second == cmpit->second)
							cout << "Parameter " << search << " is the same in both files." << endl;
						else {
							cout << "Parameter " << search << "differs: " << endl
							     << procparFile << ": " << it->second << endl
								 << cmpFile << ": " << it->second << endl;
						}
					} else {
						cout << "Parameter " << search << " present in "
						     << procparFile << " but not " << cmpFile << endl;
					}
				}
				else
					cout << it->second << endl;
			}
		}
	} else if (all)	{
		if (list) {
			if (verbose)
				cout << "Listing all parameter names." << endl;
			ParameterList::iterator it;
			for (it = pars.begin(); it != pars.end(); it++)
				cout << it->second.name() << endl;
		} else {
			if (verbose)
				cout << "Searching all parameters in " << procparFile << endl;
			ParameterList::iterator it, cmpit;
			for (it = pars.begin(); it != pars.end(); it++) {
				if (cmpFile != "") {
					cmpit = cmpPars.find(it->first);
					if (cmpit != cmpPars.end()) {
						if (it->second == cmpit->second)
							cout << "Parameter " << it->first << " is the same in both files." << endl;
						else {
							cout << "Parameter " << it->first << "differs: " << endl
							     << procparFile << ": " << it->second << endl
								 << cmpFile << ": " << it->second << endl;
						}
					} else {
						cout << "Parameter " << it->first << " present in "
						     << procparFile << " but not " << cmpFile << endl;
					}
				}
				else
					cout << it->second << endl;
			}
		}
	} else {
		string search;
		ParameterList::iterator it;
		while (true) {
			if (verbose)
				cout << "Enter parameter name (. to exit): ";
			if (!(cin >> search))
				break;
			if (search == ".")
				break;
			if (list) {
				vector<string> patterns;
				patterns.push_back(search);
				listMatchingPars(pars, patterns);
			} else {
				it = pars.find(search);
				if (it == pars.end())
					cout << search << " not found." << endl;
				else if (fullPrint)
					cout << it->second;
				else {
					cout << it->second.name() << " type: " << it->second.subtype_name() << " with " << it->second.nvals() << " values: " << endl
					     << it->second.print_values() << endl;
				}
			}
		}
	}
	
    return EXIT_SUCCESS;
}

