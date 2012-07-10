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

int main(int argc, char **argv)
{
	static int fullPrint, allPrint;
	static struct option long_options[] =
	{
		{"full", no_argument, &fullPrint, true},
		{"all",  no_argument, &allPrint, true},
		{0, 0, 0, 0}
	};
	
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "m:z", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
		}
	}
	
	fprintf(stdout, "Reading procpar file: %s\n", argv[optind]);
	par_t *pars = readProcpar(argv[optind]);
	int npar = HASH_COUNT(pars);
	fprintf(stdout, "Read %d parameters.\n", npar);
	if (allPrint)
	{
		par_t *par, *tmp;
		HASH_ITER(hh, pars, par, tmp)
		{
			fprintPar(stdout, par);
		}
	}
	else
	{
		char search[MAXSTR];
		par_t *par;
		while (true)
		{
			fprintf(stdout, "Enter parameter name (. to exit): ");
			fscanf(stdin, "%s", search);
			if (strcmp(search, ".") == 0)
				break;
			HASH_FIND_STR(pars, search, par);
			if (!par)
				fprintf(stdout, "%s not found.\n", search);
			else if (fullPrint)
				fprintPar(stdout, par);
			else
			{
				fprintf(stdout, "Type: %s with %d values\n", parTypeStr(par), par->nvals);
				fprintVals(stdout, par);
			}

		}
	}
	
    return EXIT_SUCCESS;
}

