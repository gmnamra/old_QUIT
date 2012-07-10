//
//  procpar.h
//  procparse
//
//  Created by Tobias Wood on 10/07/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#ifndef procparse_procpar_h
#define procparse_procpar_h

#include <stdio.h>
#include "uthash.h"

#define MAXSTR 1024

typedef struct par_s
{
	char *name;
	int  subtype;
	int  type;
	double  max;
	double  min;
	double  step;
	int  ggroup;
	int  dgroup;
	int  protection;
	int  active;
	int  intptr;
	void *vals;
	void *allowed;
	int nvals;
	int nallowed;
	UT_hash_handle hh;
} par_t;

par_t *readProcpar(const char *filename);
par_t *readPar(FILE *in);
void fprintPar(FILE *f, par_t *p);
void fprintVals(FILE *f, par_t *p);
void fprintAllowedVals(FILE *f, par_t *p);
const char *parTypeStr(par_t *par);

double realVal(par_t *pars, char *name, int i);
int    intVal(par_t *pars, char *name, int i);
char   *stringVal(par_t *pars, char*name, int i);

double *realVals(par_t *pars, char *name, int *nvals);
int    *intVals(par_t *pars, char *name, int *nvals);
char   **stringVals(par_t *pars, char *name, int *nvals);

#endif
