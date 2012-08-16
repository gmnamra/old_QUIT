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
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/NonLinearOptimization>

#include "DESPOT_Functors.h"

#include <stdio.h>

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
//******************************************************************************
#pragma mark Parameter Files
//******************************************************************************
int fgetArray(FILE *in, char type, size_t n, void *array);

void linearLeastSquares(double *X, double *Y, int nD,
						double *slope, double *inter, double *res);
double classicDESPOT1(const ArrayXd &flipAngles, const ArrayXd &spgrVals,
				      double TR, double B1, double *M0, double *T1);
double classicDESPOT2(const ArrayXd &flipAngles, const ArrayXd &ssfpVals,
                      double TR, double T1, double B1, double *M0, double *T2);

typedef Array<bool, Dynamic, Dynamic> ArrayXb;
typedef std::pair<int, double> argsort_pair;

bool argsort_comp(const argsort_pair& left, const argsort_pair& right);

template<typename Derived>
ArrayXi arg_partial_sort(const ArrayBase<Derived> &x, int middle)
{
    ArrayXi indices(middle);
    std::vector<argsort_pair> data(x.size());
    for(int i=0;i<x.size();i++)
	{
        data[i].first = i;
        data[i].second = x(i);
    }
    std::partial_sort(data.begin(), data.begin() + middle, data.end(), argsort_comp);
    for(int i=0;i<middle;i++) {
        indices(i) = data[i].first;
    }    
    return indices;
}

double regionContraction(ArrayXd &params, SPGR_2c SPGR, SSFP_2c SSFP,
                         const ArrayXd &loStart, const ArrayXd &hiStart,
					     const ArrayXi &loConstrained, const ArrayXi &hiConstrained,
					     const int nS, const int nR, const int maxContractions,
						 const double thresh, const double expand);

#endif
