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

void freeProcpar(par_t *p)
{
	if (p)
	{
		par_t *current, *tmp;
		HASH_ITER(hh, p, current, tmp)
		{
			free(current->vals);
			if (current->nallowed > 0)
				free(current->allowed);
			HASH_DEL(p, current); 
		}
		p = NULL;
	}
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
	
	switch (par.type)
	{
		// Even ints have to be read as reals (doubles) because VnmrJ doesn't
		// output them as ints.
		case TYPE_REAL:
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
		case TYPE_STRING:
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
	}
	
	if (fscanf(in, "%d ", &(par.nallowed)) == EOF)
		return NULL;
	if (par.nallowed > 0)
	{
		switch (par.type)
		{
			case TYPE_REAL:
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
			case TYPE_STRING:
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
		}
	}
	
	// Got to here without hitting EOF, have a complete parameter
	// so can safely allocate some memory without leaking
	par_t *p = malloc(sizeof(par_t));
	memcpy(p, &par, sizeof(par_t));
	return p;
}

par_t *createPar(const char *name, int subtype, int nvals, void *vals)
{
	par_t *par = calloc(1, sizeof(par_t));	// Local, make sure we have a complete par before returning
	par->name = malloc((strlen(name) + 1) * sizeof(char));
	strcpy(par->name, name);
	par->subtype = subtype;
	switch (subtype)
	{
		case PAR_REAL:
		case PAR_DELAY:
		case PAR_FREQ:
		case PAR_PULSE:
		case PAR_INT:
			par->type = TYPE_REAL;
			break;
		case PAR_STRING:
		case PAR_FLAG:
			par->type = TYPE_STRING;
			break;
	}
	setPar(par, nvals, vals);
	return par;
}

void setPar(par_t *par, int nvals, void *vals)
{
	par->nvals = nvals;
	switch (par->type)
	{
		case TYPE_REAL:
		{
			if (par->vals)
				free(par->vals);
			par->vals = malloc(par->nvals * sizeof(double));
			memcpy(par->vals, vals, par->nvals * sizeof(double));
		}
		break;
		case TYPE_STRING:
			if (par->vals)
			{
				for (size_t i = 0; i < par->nvals; i++)
					free(((char **)par->vals)[i]);
				free(par->vals);
			}
			par->vals = malloc(par->nvals * sizeof(char *));
			for (size_t i = 0; i < par->nvals; i++)
			{
				((char **)par->vals)[i] = malloc((strlen(((char **)vals)[i]) + 1) * sizeof(char));
				strcpy(((char **)par->vals)[i], ((char **)vals)[i]);
			}
		break;
	}
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
		switch (p->type)
		{
			case TYPE_REAL:
				fprintf(f, "%g ", ((double*)p->vals)[i]);
				break;
			case TYPE_STRING:
				fprintf(f, "\"%s\"\n", ((char **)p->vals)[i]);
				break;	
		}
	}
	if (p->type == TYPE_REAL)
		fprintf(f, "\n");
}

void fprintAllowedVals(FILE *f, par_t *p)
{
	for (size_t i = 0; i < p->nallowed; i++)
	{
		switch (p->type)
		{
			case TYPE_REAL:
				fprintf(f, "%g ", ((double*)p->allowed)[i]);
				break;
			case TYPE_STRING:
				fprintf(f, "\"%s\" ", ((char **)p->allowed)[i]);
				break;	
		}
	}
	fprintf(f, "\n");
}

const char *parTypeStr(int type)
{
	switch (type)
	{
		case TYPE_REAL:
		    return "real";
		case TYPE_STRING:
			return "string";
		default:
			return "INVALID";
	    break;
	}
}

const char *parSubTypeStr(int subtype)
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

double realVal(par_t *pars, const char *name, int i)
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

char *stringVal(par_t *pars, const char*name, int i)
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

double *realVals(par_t *pars, const char *name, int *nvals)
{
	par_t *p;
	HASH_FIND_STR(pars, name, p);
	if (p)
	{
		if (p->type == TYPE_REAL)
		{
			if (nvals)
				*nvals = p->nvals;
			return (double *)p->vals;
		}
		fprintf(stderr, "Parameter %s is not a real.\n", name);
	}
	else
		fprintf(stderr, "Parameter %s not found.\n", name);
	return NULL;
}

char **stringVals(par_t *pars, const char *name, int *nvals)
{
	par_t *p;
	HASH_FIND_STR(pars, name, p);
	if (p)
	{
		if (p->type == TYPE_STRING)
		{
			if (nvals)
				*nvals = p->nvals;
			return (char **)p->vals;
		}
		fprintf(stderr, "Parameter %s is not a string.\n", name);
	}
	else
		fprintf(stderr, "Parameter %s not found.\n", name);
	return NULL;
}