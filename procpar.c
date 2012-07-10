//
//  procpar.c
//  procparse
//
//  Created by Tobias Wood on 10/07/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#include "procpar.h"

par_t *readProcpar(const char *filename)
{
	FILE *procpar = fopen(filename, "r");
	if (!procpar)
	{
		fprintf(stderr, "Cannot open procpar file: %s\n", filename);
		return NULL;
	}
	par_t *par, *pars = NULL;
	while ((par = readPar(procpar)))
		HASH_ADD_KEYPTR(hh, pars, par->name, strlen(par->name), par);
	return pars;
}



par_t *readPar(FILE *in)
{
	par_t par;	// Local, make sure we have a complete par before returning
	char temp[MAXSTR];
	
	if (fscanf(in, "%s", temp) == EOF)
		return NULL;
	par.name = malloc((strlen(temp) + 1) * sizeof(char));
	strcpy(par.name, temp);
	
	if (fscanf(in, "%d %d %lf %lf %lf %d %d %d %d %d",
	            &(par.subtype), &(par.type),
		        &(par.max), &(par.min), &(par.step),
		        &(par.ggroup), &(par.dgroup),
		        &(par.protection), &(par.active), &(par.intptr)) == EOF)
		return NULL;
	
	// Make sure to eat the space after nvals, or string reading breaks
	if (fscanf(in, "%d ", &(par.nvals)) == EOF)
		return NULL;
	
	switch (par.subtype)
	{
		case PAR_REAL: case PAR_DELAY: case PAR_FREQ: case PAR_PULSE:
		{
			par.vals = malloc(par.nvals * sizeof(double));
			double inVal;
			for (size_t i = 0; i < par.nvals; i++)
			{
				fscanf(in, "%lf", &inVal);
				((double *)par.vals)[i] = inVal;
			}
			if (feof(in))
				return NULL;
		}
		break;
		case PAR_STRING: case PAR_FLAG:
			par.vals = malloc(par.nvals * sizeof(char *));
			char nextString[MAXSTR];
			for (size_t i = 0; i < par.nvals; i++)
			{
				fgets(nextString, MAXSTR, in);
				char *nextQuote = strtok(nextString, "\"");
				// nextQuote = strtok(NULL, "\""); // Eat first "
				// strtok does not neatly detect null tokens (ie "")
				// It sets the first delimiter to \0 and moves the pointer
				// on to the first non-delimiter character. For vals, which are
				// written 1 per line, can just check for the trailing newline
				if (*nextQuote == '\n') // We have an empty string
					((char **)par.vals)[i] = "";
				else
				{
					((char **)par.vals)[i] = malloc((strlen(nextQuote) + 1) * sizeof(char));
					strcpy(((char **)par.vals)[i], nextQuote);
				}
			}
		break;
		case PAR_INT:
		{
			par.vals = malloc(par.nvals * sizeof(int));
			int inVal;
			for (size_t i = 0; i < par.nvals; i++)
			{
				fscanf(in, "%d", &inVal);
				((int *)par.vals)[i] = inVal;	
			}
			if (feof(in))
				return NULL;
		}
		break;
	}
	
	if (fscanf(in, "%d ", &(par.nallowed)) == EOF)
		return NULL;
	if (par.nallowed > 0)
	{
		switch (par.subtype)
		{
			case PAR_REAL: case PAR_DELAY: case PAR_FREQ: case PAR_PULSE:
			{
				par.allowed = malloc(par.nallowed * sizeof(double));
				double inVal;
				for (size_t i = 0; i < par.nallowed; i++)
				{
					fscanf(in, "%lf", &inVal);
					((double *)par.allowed)[i] = inVal;
				}
				if (feof(in))
					return NULL;
			}
			break;
			case PAR_STRING: case PAR_FLAG:
			{
				par.allowed = malloc(par.nallowed * sizeof(char *));
				char restOfLine[MAXSTR];
				fgets(restOfLine, MAXSTR, in);
				char *nextQuote = strtok(restOfLine, "\"");
				char *lastQuote = nextQuote;
				for (size_t i = 0; i < par.nallowed; i++)
				{
					// See above about strtok
					// If we've got more than 3 chars difference, there was an
					// empty string	
					if (lastQuote - nextQuote > 3)
						((char **)par.allowed)[i] = "";
					else
					{
						((char **)par.allowed)[i] = malloc((strlen(nextQuote) + 1) * sizeof(char));
						strcpy(((char **)par.allowed)[i], nextQuote);
					}
					// Eat spaces in between "
					lastQuote = nextQuote;
					nextQuote = strtok(NULL, "\" ");
				}
			}
			break;
			case PAR_INT: // Integer
			{
				par.allowed = malloc(par.nallowed * sizeof(int));
				int inVal;
				for (size_t i = 0; i < par.nallowed; i++)
				{
					fscanf(in, "%d", &inVal);
					((int *)par.allowed)[i] = inVal;
				}
				if (feof(in))
					return NULL;
			}
			break;
		}
	}
	
	// Got to here without hitting EOF, have a complete parameter
	// so can safely allocate some memory without leaking
	par_t *p = malloc(sizeof(par_t));
	memcpy(p, &par, sizeof(par_t));
	return p;
}

void fprintPar(FILE* f, par_t *p)
{
	fprintf(f, "%s %d %d %g %g %g %d %d %d %d %d\n%d ",
	        p->name, p->subtype, p->type,
		    p->max, p->min, p->step,
		    p->ggroup, p->dgroup,
		    p->protection, p->active, p->intptr,
			p->nvals);
	fprintVals(f, p);
	fprintf(f, "%d ", p->nallowed);
	fprintAllowedVals(f, p);
}

void fprintVals(FILE* f, par_t *p)
{
	for (size_t i = 0; i < p->nvals; i++)
	{
		switch (p->subtype)
		{
			case PAR_REAL: case PAR_DELAY: case PAR_FREQ: case PAR_PULSE:
				fprintf(f, "%g ", ((double*)p->vals)[i]);
				break;
			case PAR_STRING: case PAR_FLAG:
				fprintf(f, "\"%s\"\n", ((char **)p->vals)[i]);
				break;
			case PAR_INT:
				fprintf(f, "%d ", ((int *)p->vals)[i]);
				break;	
		}
	}
	if (p->subtype != 2 && p->subtype != 4)
		fprintf(f, "\n");
}

void fprintAllowedVals(FILE *f, par_t *p)
{
	for (size_t i = 0; i < p->nallowed; i++)
	{
		switch (p->subtype)
		{
			case PAR_REAL: case PAR_DELAY: case PAR_FREQ: case PAR_PULSE:
				fprintf(f, "%g ", ((double*)p->allowed)[i]);
				break;
			case PAR_STRING: case PAR_FLAG:
				fprintf(f, "\"%s\" ", ((char **)p->allowed)[i]);
				break;
			case PAR_INT:
				fprintf(f, "%d ", ((int *)p->allowed)[i]);
				break;	
		}
	}
	fprintf(f, "\n");
}

const char *parTypeStr(int subtype)
{
	switch (subtype)
	{
		case PAR_REAL:
    		return "real";
		case PAR_STRING:
			return "string";
		case PAR_DELAY:
			return "delay";
		case PAR_FLAG:
			return "flag";
		case PAR_FREQ:
			return "frequency";
		case PAR_PULSE:
			return "pulse";
		case PAR_INT:
			return "integer";
		default:
			return "INVALID";
	}
}

double realVal(par_t *pars, char *name, int i)
{
	par_t *p;
	HASH_FIND_STR(pars, name, p);
	if (p)
	{
		if (i < p->nvals)
			return ((double *)p->vals)[i];
		else
			fprintf(stderr, "Tried to access element %d of parameter %s, only %d values.\n",
		    	            i, name, p->nvals);
	}
	else
		fprintf(stderr, "Parameter %s not found.\n", name);
	return 0;
}

int intVal(par_t *pars, char *name, int i)
{
	par_t *p;
	HASH_FIND_STR(pars, name, p);
	if (p)
	{
		if (i < p->nvals)
			return ((int *)p->vals)[i];
		else
			fprintf(stderr, "Tried to access element %d of parameter %s, only %d values.\n",
		    	            i, name, p->nvals);
	}
	else
		fprintf(stderr, "Parameter %s not found.\n", name);
	return 0;
}

char *stringVal(par_t *pars, char*name, int i)
{
	par_t *p;
	HASH_FIND_STR(pars, name, p);
	if (p)
	{
		if (i < p->nvals)
			return ((char **)p->vals)[i];
		else
			fprintf(stderr, "Tried to access element %d of parameter %s, only %d values.\n",
		    	            i, name, p->nvals);
	}
	else
		fprintf(stderr, "Parameter %s not found.\n", name);
	return 0;
}

double *realVals(par_t *pars, char *name, int *nvals)
{
	par_t *p;
	HASH_FIND_STR(pars, name, p);
	if (p)
	{
		if ((p->subtype == PAR_REAL) ||
		    (p->subtype == PAR_DELAY) ||
			(p->subtype == PAR_FREQ) ||
			(p->subtype == PAR_PULSE))
		{
			if (nvals)
				*nvals = p->nvals;
			return (double *)p->vals;
		}
		fprintf(stderr, "Parameter %s is of type %s, not %s, %s, %s or %s.\n",
		        name, parTypeStr(p->subtype), parTypeStr(PAR_REAL),
				parTypeStr(PAR_DELAY), parTypeStr(PAR_FREQ), parTypeStr(PAR_PULSE));
	}
	else
		fprintf(stderr, "Parameter %s not found.\n", name);
	return NULL;
}

int *intVals(par_t *pars, char *name, int *nvals)
{
	par_t *p;
	HASH_FIND_STR(pars, name, p);
	if (p)
	{
		if (p->subtype == PAR_INT)
		{
			if (nvals)
				*nvals = p->nvals;
			return (int *)p->vals;
		}
		fprintf(stderr, "Parameter %s is of type %s, not %s.\n",
		        p->name, parTypeStr(p->subtype), parTypeStr(PAR_INT));
	}
	else
		fprintf(stderr, "Parameter %s not found.\n", name);
	return NULL;
}

char **stringVals(par_t *pars, char *name, int *nvals)
{
	par_t *p;
	HASH_FIND_STR(pars, name, p);
	if (p)
	{
		if ((p->subtype == PAR_STRING) || (p->subtype == PAR_FLAG))
		{
			if (nvals)
				*nvals = p->nvals;
			return (char **)p->vals;
		}
		fprintf(stderr, "Parameter %s is of type %s, not %s or %s.\n",
		        name, parTypeStr(p->subtype),
				parTypeStr(PAR_STRING), parTypeStr(PAR_FLAG));
	}
	else
		fprintf(stderr, "Parameter %s not found.\n", name);
	return NULL;
}