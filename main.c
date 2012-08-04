//
//  main.c
//  procparse
//
//  Created by Tobias Wood on 10/07/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#include <stdio.h>
#include <stdbool.h>
#include <getopt.h>
#include "procpar.h"

static int fullPrint = false, all = false, list = false, verbose = false;

bool comparePars(par_t *par, par_t *cmpPar)
{
	if (par->nvals != cmpPar->nvals)
	{
		fprintf(stdout, "Parameter %s differs in number of values: %d vs. %d.\n",
				par->name, par->nvals, cmpPar->nvals);
		return false;
	}
	else
	{
		switch (par->subtype)
		{
			case 1: case 3:	case 5:	case 6:
			{
				double *pvals = (double *)par->vals;
				double *cmpvals = (double *)cmpPar->vals;
				for (int v = 0; v < par->nvals; v++)
				{
					if (pvals[v] != cmpvals[v])
					{
						fprintf(stdout, "Parameter %s[%d] differs: %f vs. %f.\n",
								par->name, v, pvals[v], cmpvals[v]);
						return false;
					}
				}
			}
			break;
			case 2:	case 4:
			{
				char **pvals = (char **)par->vals;
				char **cmpvals = (char **)cmpPar->vals;
				for (int v = 0; v < par->nvals; v++)
				{
					if (strcmp(pvals[v], cmpvals[v]) != 0)
					{
						fprintf(stdout, "Parameter %s[%d] differs: \"%s\" vs. \"%s\"\n",
								par->name, v, pvals[v], cmpvals[v]);
						return false;
					}
				}
			}
			break;
			case 7:
			{
				int *pvals = (int *)par->vals;
				int *cmpvals = (int *)cmpPar->vals;
				for (int v = 0; v < par->nvals; v++)
				{
					if (pvals[v] != cmpvals[v])
					{
						fprintf(stdout, "Parameter %s[%d] differs: %d vs. %d.\n",
								par->name, v, pvals[v], cmpvals[v]);
						return false;
					}
				}
			}
			break;
		}
	}
	return true;
}

void listMatchingPars(par_t *pars, char **patterns, int n);
void listMatchingPars(par_t *pars, char **patterns, int n)
{
	for (int i = 0; i < n; i++)
	{
		if (verbose)
			fprintf(stdout, "Listing parameter names containing \"%s\":\n", patterns[i]);
		par_t *current, *tmp;
		int matched = 0;
		HASH_ITER(hh, pars, current, tmp)
		{
			if (strstr(current->name, patterns[i]))
			{
				fprintf(stdout, "%s\n", current->name);
				matched++;
			}
		}
		if (verbose)
			fprintf(stdout, "Found %d matches.\n\n", matched);
	}
}

const char* usage = "procparse - A utility to find interesting information in Agilent procpar files.\n\
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
	
	char *procparFile = NULL, *cmpFile = NULL;
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "hv", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 0:
				// It was an option that just sets a flag.
				break;
			case 'c':
				cmpFile = optarg;
				break;
			case 'v':
				verbose = true;
				break;
			default:
				fprintf(stdout, "Unknown option %s.\n", long_options[indexptr].name);
			case 'h':
				fprintf(stdout, "%s", usage);
				break;
		}
	}
	
	par_t *pars, *cmpPars;
	
	procparFile = argv[optind++];
	if (verbose)
		fprintf(stdout, "Reading procpar file: %s\n", procparFile);
	pars = readProcpar(procparFile);
	if (!pars)
		exit(EXIT_FAILURE);
	if (cmpFile)
	{
		if (verbose)
			fprintf(stdout, "Reading comparison procpar file: %s\n", cmpFile);
		cmpPars = readProcpar(cmpFile);
		if (!cmpPars)
			exit(EXIT_FAILURE);
	}
	
	par_t *par, *cmpPar, *tmp;
	
	if (optind < argc)
	{
		if (list)
			listMatchingPars(pars, argv + optind, (argc - optind));
		else
		{
			if (verbose)
				fprintf(stdout, "Searching for %d parameters.\n", argc - optind);
			for (int i = optind; i < argc; i++)
			{
				HASH_FIND_STR(pars, argv[i], par);
				if (cmpFile)
				{
					HASH_FIND_STR(cmpPars, par->name, cmpPar);
					if (cmpPar)
						comparePars(par, cmpPar);
					else
						fprintf(stdout, "Parameter %s present in %s but not %s.\n",
								par->name, procparFile, cmpFile);
				}
				else
					fprintPar(stdout, par);
			}
		}
	}
	else if (all)
	{
		if (list)
		{
			if (verbose)
				fprintf(stdout, "Listing all parameter names.\n");
			HASH_ITER(hh, pars, par, tmp)
				fprintf(stdout, "%s\n", par->name);
		}
		else
		{
			if (verbose)
				fprintf(stdout, "Searching all parameters in %s.\n", procparFile);
			HASH_ITER(hh, pars, par, tmp)
			{
				if (cmpFile)
				{
					HASH_FIND_STR(cmpPars, par->name, cmpPar);
					if (cmpPar)
						comparePars(par, cmpPar);
					else
						fprintf(stdout, "Parameter %s present in %s but not %s.\n",
								par->name, procparFile, cmpFile);
				}
				else
					fprintPar(stdout, par);
			}
		}
	}
	else
	{
		char *search = malloc(MAXSTR * sizeof(char));
		par_t *par;
		while (true)
		{
			if (verbose)
				fprintf(stdout, "Enter parameter name (. to exit): ");
			if (fscanf(stdin, "%s", search) == EOF)
				break;
			if (strcmp(search, ".") == 0)
				break;
			if (list)
			{
				char **patterns = &search;
				listMatchingPars(pars, patterns, 1);
			}
			else
			{
				HASH_FIND_STR(pars, search, par);
				if (!par)
					fprintf(stdout, "%s not found.\n", search);
				else if (fullPrint)
					fprintPar(stdout, par);
				else
				{
					fprintf(stdout, "%s type: %s with %d values: ", par->name, parTypeStr(par->type), par->nvals);
					fprintVals(stdout, par);
				}
			}
		}
		free(search);
	}
	
    return EXIT_SUCCESS;
}

