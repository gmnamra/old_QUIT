/*
 *  mathsOptimisers.h
 *
 *  Created by Tobias Wood on 01/02/2012.
 *  Copyright 2012 Tobias Wood. All rights reserved.
 *
 */
#ifndef MATHS_OPTIM
#define MATHS_OPTIM

#include <stdlib.h>
#include <time.h>

#include "mathsUtil.h"
#include "mathsMatrix.h"

//******************************************************************************
#pragma mark Fitting Library
//******************************************************************************

typedef double (eval_type) (double x, double *parameters, double *constants);
typedef void (eval_array_type) (double *x, double *p, double *c, double *y, size_t n);
typedef void (jacob_type) (double *data, int nD, double *parameters, double *constants, double *result);
double calcResiduals(double *parameters, double *constants,
                     double *dataX, double *dataY, size_t nD,
	 	             eval_type *function, double *residuals, bool norm);
double calcAResiduals(double *parameters, double *constants,
                      double *dataX, double *dataY, size_t nD,
	 	              eval_array_type *function, double *residuals, bool norm);
double calcMResiduals(double *params, double **consts, size_t nM,
					  double **dataX, double **dataY, size_t *nD,
					  eval_type **funcs, double **residuals, bool norm);
double calcMAResiduals(double *params, double **consts, size_t nM,
					   double **dataX, double **dataY, size_t *nD,
					   eval_array_type **funcs, double **residuals,
					   bool norm);
void linearLeastSquares(double *X, double *Y, int nD,
						double *slope, double *inter, double *res);
int goldenSection(double *parameters, int nP, int P,
                  double loP, double hiP, double *constants,
                  double *dataX, double *dataY, int nD,
				  eval_type *function, double *finalResidue);
int levMar(double *parameters, size_t nP, double **constants,
           double **dataX, double **dataY, eval_array_type **funcs,
		   jacob_type **jacFuncs, size_t *nD, size_t nF,
		   double *loBounds, double *hiBounds,
		   bool normalise, double *finalResidue);
int simplex(double *params, size_t nP, double **consts, size_t nM,
		    double **dataX, double **dataY, size_t *nD,
		    eval_type **funcs, double ** initial, double *finalResidue);
void regionContraction(double *params, size_t nP, double **consts, size_t nM,
		               double **dataX, double **dataY, size_t *nD, bool norm,
					   eval_array_type **funcs, double **initBounds,
					   bool **constrained, size_t nS, size_t nR,
					   size_t maxContractions, double thresh, double expand,
					   double *finalResidue);
#endif