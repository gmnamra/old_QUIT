/*
 *  recovery.c
 *  MacRI
 *
 *  Created by Tobias Wood on 07/11/2011.
 *  Copyright 2011 Tobias Wood. All rights reserved.
 *
 */

#include "recovery.h"
#include "mathUtil.h"

inline double saturation_Mz(const double TR, const double* par, const double *c)
{
	double M0 = par[0];
	double T1 = par[1];
	double val = M0 * (1. - exp(-TR / T1));
	return val;
}

inline double inversion_Mz(const double TR, const double *par, const double *c)
{
	double M0 = par[0];
	double T1 = par[1];
	double val = M0 * (1. - 2. * exp(-TR / T1));
	return val;
}

inline double recovery_Mz(const double TR, const double *par, const double *c)
{
	double M0 = par[0];
	double T1 = par[1];
	double alpha = par[2];
	double val = M0 * (1. - alpha * (exp(-TR / T1)));
	return val;
}

inline double recovery_dMzdM0(const double TR, const double *par, const double *c)
{
	//double M0 = par[0];
	double T1 = par[1];
	double alpha = par[2];
	double val = (1. - alpha * (exp(-TR / T1)));
	return val;
}

inline double recovery_dMzdT1(const double TR, const double *par, const double *c)
{
	double M0 = par[0];
	double T1 = par[1];
	double alpha = par[2];
	double val = - M0 * alpha * TR * (exp(-TR / T1)) / (T1 * T1);
	return val;
}

inline double recovery_dMzdalpha(const double TR, const double *par, const double *c)
{
	double M0 = par[0];
	double T1 = par[1];
	//double alpha = par[2];
	double val = - M0 * (exp(-TR / T1));
	return val;
}

extern int MATH_LEVENBERG_DEBUG;
double calcRecovery(double *vals, double* TR, int n, double *M0out, double *T1out, double *alpha)
{
	// Initial guesses of M0 and T1
	// First value should be close to last Mz, T1 is just a guess
	int n_par = 3;
	double par[3] = {1.2 * vals[n - 1], 1200., 1.}; 
	eval_type *gradients[3] = { recovery_dMzdM0,
	                            recovery_dMzdT1,
								recovery_dMzdalpha};
	double finalRes;
	MATH_LEVENBERG_DEBUG = 0;
	levMar(par, n_par, NULL, TR, vals, n, &recovery_Mz, gradients, &finalRes);
	*M0out = par[0];
	*T1out = par[1];
	*alpha = par[2];
	return finalRes;
}