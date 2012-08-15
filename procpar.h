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

#ifdef __cplusplus
extern "C"{
#endif
#define MAXSTR 1024

enum PAR_TYPES
{
	TYPE_REAL = 1,
	TYPE_STRING
};

enum PAR_SUBTYPES
{
	PAR_REAL = 1,
	PAR_STRING,
	PAR_DELAY,
	PAR_FLAG,
	PAR_FREQ,
	PAR_PULSE,
	PAR_INT
};

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
void  freeProcpar(par_t *p);
par_t *readPar(FILE *in);
void fprintPar(FILE *f, par_t *p);
void fprintVals(FILE *f, par_t *p);
void fprintAllowedVals(FILE *f, par_t *p);
const char *parTypeStr(int type);
const char *parSubTypeStr(int subtype);
double realVal(par_t *pars, char *name, int i);
char   *stringVal(par_t *pars, char*name, int i);

double *realVals(par_t *pars, char *name, int *nvals);
char   **stringVals(par_t *pars, char *name, int *nvals);

#ifdef __cplusplus
}
#endif
#endif
