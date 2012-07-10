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

int main(int argc, char **argv)
{
	static int fullPrint, all;
	static struct option long_options[] =
	{
		{"full", no_argument, &fullPrint, true},
		{"all",  no_argument, &all, true},
		{"cmp",  required_argument, 0, 'c'},
		{0, 0, 0, 0}
	};
	
	char *procparFile = NULL, *cmpFile = NULL;
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "m:z", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 'c':
				cmpFile = optarg;
			break;
		}
	}
	
	par_t *pars, *cmpPars;
	
	procparFile = argv[optind++];
	fprintf(stdout, "Reading procpar file: %s\n", procparFile);
	pars = readProcpar(procparFile);
	if (cmpFile)
	{
		fprintf(stdout, "Reading comparison procpar file: %s\n", cmpFile);
		cmpPars = readProcpar(cmpFile);
	}
	
	par_t *par, *cmpPar, *tmp;
	
	if (optind < argc)
	{
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
	else if (all)
	{
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
	else
		fprintf(stdout, "No parameters specified.\n");
	
    return EXIT_SUCCESS;
}

