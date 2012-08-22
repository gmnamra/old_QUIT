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

#include <iostream>

#include <Eigen/Dense>
using namespace Eigen;

//******************************************************************************
#pragma mark Convenience stuff
//******************************************************************************
#ifndef TRUE
	#define TRUE 1
#endif
#ifndef FALSE
	#define FALSE 0
#endif
#ifndef M_PI
	#define M_PI 3.14159265358979323846264338327950288
#endif
#ifndef M_SQRT2
	#define M_SQRT2     1.41421356237309504880168872420969808
#endif
extern int MATHS_DEBUG;

double radians(double degrees);
double degrees(double radians);

double clamp(double value, double low, double high);

void linearLeastSquares(double *X, double *Y, int nD,
						double *slope, double *inter, double *res);
double classicDESPOT1(const ArrayXd &flipAngles, const ArrayXd &spgrVals,
				      double TR, double B1, double *M0, double *T1);
double classicDESPOT2(const ArrayXd &flipAngles, const ArrayXd &ssfpVals,
                      double TR, double T1, double B1, double *M0, double *T2);

#endif
