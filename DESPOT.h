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
double clamp(double value, double low, double high);

void linearLeastSquares(double *X, double *Y, long nD,
						double *slope, double *inter, double *res);
double classicDESPOT1(const ArrayXd &flipAngles, const ArrayXd &spgrVals,
				      double TR, double B1, double *M0, double *T1);
double classicDESPOT2(const ArrayXd &flipAngles, const ArrayXd &ssfpVals,
                      double TR, double T1, double B1, double *M0, double *T2);

#endif
