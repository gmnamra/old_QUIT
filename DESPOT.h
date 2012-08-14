/*
 *  DESPOT1.h
 *  MacRI
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright 2011 Tobias Wood. All rights reserved.
 *
 */
#ifndef __DESPOT__
#define __DESPOT__

#include "mathsMatrix.h"
#include "mathsOptimisers.h"

#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
using namespace Eigen;
int tests();
double SPGR(double flipAngle, double *p, double *c);
void SPGR_Jacobian(double *angles, int nD, double *p, double *c, double *result);
double IRSPGR(double TI, double *p, double *c);
void IRSPGR_Jacobian(double *data, int nD, double *par, double *c, double *result);

typedef struct
{
	double TR, M0, B1;
} SPGR_constants;

typedef struct
{
	double TR, M0, T1, B0, B1, rfPhase;
} SSFP_constants;

eval_type a1cSSFP;
eval_type a2cSPGR;
eval_type a2cSSFP;

double calcDESPOT1(double *flipAngles, double *spgrVals, int n,
				   double TR, double B1, double *M0, double *T1);
double classicDESPOT2(gsl_vector *flipAngles, gsl_vector *ssfpVals,
                      double TR, double T1, double B1, double *M0, double *T2);
void simplexDESPOT2(size_t nPhases, size_t *nD, double *phases, double **angles, double **ssfp,
					double TR, double T1, double B1, double *p);
#endif