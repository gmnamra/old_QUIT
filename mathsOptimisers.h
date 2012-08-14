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

#include "mathsMatrix.h"
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort_vector.h>

//******************************************************************************
#pragma mark Fitting Library
//******************************************************************************
typedef void (eval_type)  (gsl_vector *x, gsl_vector *p, void *c, gsl_vector *y);
typedef void (jacob_type) (gsl_vector *x, gsl_vector *p, void *c, gsl_matrix *dydx);
double calcResiduals(gsl_vector *parameters, void *constants,
					 gsl_vector *dataX, gsl_vector *dataY,
					 eval_type  *function, gsl_vector *residuals);
double calcMResiduals(gsl_vector *params, size_t nF, void **consts,
                      gsl_vector **dataX, gsl_vector **dataY,
					  eval_type  **funcs, gsl_vector **residuals);
void linearLeastSquares(double *X, double *Y, int nD,
						double *slope, double *inter, double *res);
int goldenSection(gsl_vector *parameters, int P,
                  double loP, double hiP, void *constants,
                  gsl_vector *dataX, gsl_vector *dataY,
				  eval_type *function, double *finalResidue);
extern int LEV_DEBUG;
int levenbergMarquardt(gsl_vector *parameters, size_t nF, void **constants,
                       gsl_vector **dataX, gsl_vector **dataY,
					   eval_type  **funcs, jacob_type **jacFuncs,
					   double *finalResidue);
int simplex(gsl_vector *params, size_t nF, void **consts,
		    gsl_vector **dataX, gsl_vector **dataY,
		    eval_type  **funcs, gsl_vector **initial, double *finalResidue);
void regionContraction(gsl_vector *params, size_t nF, void **consts,
					   gsl_vector **dataX, gsl_vector **dataY,
					   eval_type  **funcs, gsl_vector *initBounds[2],
					   bool **constrained, size_t nS, size_t nR,
					   size_t maxContractions, double thresh, double expand,
					   double *finalResidue);
#endif