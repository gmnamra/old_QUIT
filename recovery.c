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

inline float saturation_Mz(const float TR, const float* par, const float *c)
{
	float M0 = par[0];
	float T1 = par[1];
	float val = M0 * (1. - exp(-TR / T1));
	return val;
}

inline float inversion_Mz(const float TR, const float *par, const float *c)
{
	float M0 = par[0];
	float T1 = par[1];
	float val = M0 * (1. - 2. * exp(-TR / T1));
	return val;
}

inline float recovery_Mz(const float TR, const float *par, const float *c)
{
	float M0 = par[0];
	float T1 = par[1];
	float alpha = par[2];
	float val = M0 * (1. - alpha * (exp(-TR / T1)));
	return val;
}

inline float recovery_dMzdM0(const float TR, const float *par, const float *c)
{
	//float M0 = par[0];
	float T1 = par[1];
	float alpha = par[2];
	float val = (1. - alpha * (exp(-TR / T1)));
	return val;
}

inline float recovery_dMzdT1(const float TR, const float *par, const float *c)
{
	float M0 = par[0];
	float T1 = par[1];
	float alpha = par[2];
	float val = - M0 * alpha * TR * (exp(-TR / T1)) / (T1 * T1);
	return val;
}

inline float recovery_dMzdalpha(const float TR, const float *par, const float *c)
{
	float M0 = par[0];
	float T1 = par[1];
	//float alpha = par[2];
	float val = - M0 * (exp(-TR / T1));
	return val;
}

extern int MATH_LEVENBERG_DEBUG;
float calcRecovery(float *vals, float* TR, int n, float *M0out, float *T1out, float *alpha)
{
	// Initial guesses of M0 and T1
	// First value should be close to last Mz, T1 is just a guess
	int n_par = 3;
	float par[3] = {1.2 * vals[n - 1], 1200., 1.}; 
	eval_type *gradients[3] = { recovery_dMzdM0,
	                            recovery_dMzdT1,
								recovery_dMzdalpha};
	float finalRes;
	MATH_LEVENBERG_DEBUG = 0;
	//levMar(par, n_par, NULL, TR, vals, n, &recovery_Mz, gradients, &finalRes);
	*M0out = par[0];
	*T1out = par[1];
	*alpha = par[2];
	return finalRes;
}