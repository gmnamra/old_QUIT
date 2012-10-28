/*
 *  DESPOT1.c
 *  MacRI
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright 2011 Tobias Wood. All rights reserved.
 *
 */

#include "DESPOT.h"

double clamp(double value, double low, double high)
{
	if (value < low)
		return low;
	if (value > high)
		return high;
	return value;
}

void apply_for(const int max, const function<void(int)> f, const int num_threads) {
	vector<thread> pool;
	
	function<void(int)> worker = [&max, &f, &num_threads](int local) {
		while (local < max) {
			f(local);
			local += num_threads;
		}
	};
	
	for (int t = 0; t < num_threads; t++)
		pool.push_back(thread(worker, t));
	for (int t = 0; t < pool.size(); t++)
		pool[t].join();

}

//******************************************************************************
// Basic least squares fitting
//******************************************************************************
void linearLeastSquares(double *X, double *Y, long nD,
						double *slope, double *inter, double *res)
{
	double sumX, sumY, sumXX, sumXY;
	sumX = sumY = sumXX = sumXY = 0.;
	for (int i = 0; i < nD; i++)
	{
		double x = X[i];
		double y = Y[i];
		
		sumX  += x;
		sumY  += y;
		sumXX += (x*x);
		sumXY += (x*y);
	}
	
	*slope = (sumXY - (sumX * sumY)/nD) / (sumXX - (sumX * sumX)/nD);
	*inter = (sumY - (*slope) * sumX) / nD;
	
	if (res)
	{
		*res = 0.;
		double m = *slope; double c = *inter;
		for (int i = 0; i < nD; i++)
			*res += pow(Y[i] - (m*X[i] + c), 2.);
	}
}

double classicDESPOT1(const ArrayXd &flipAngles, const ArrayXd &spgrVals,
				      double TR, double B1, double *M0, double *T1)
{
	// Linearise the data, then least-squares
	long n = flipAngles.size();
	double X[n], Y[n], slope, inter, res;
	for (long i = 0; i < flipAngles.size(); i++)
	{
		X[i] = spgrVals[i] / tan(flipAngles[i] * B1);
		Y[i] = spgrVals[i] / sin(flipAngles[i] * B1);
	}
	linearLeastSquares(X, Y, n, &slope, &inter, &res);	
	*T1 = -TR / log(slope);
	*M0 = inter / (1. - slope);
	return res;
}

double classicDESPOT2(const ArrayXd &flipAngles, const ArrayXd &ssfpVals,
                      double TR, double T1, double B1, double *M0, double *T2)
{
	// As above, linearise, then least-squares
	long n = flipAngles.size();
	double X[n], Y[n], slope, inter, residual;
	for (long i = 0; i < flipAngles.size(); i++)
	{
		X[i] = ssfpVals[i] / tan(flipAngles[i] * B1);
		Y[i] = ssfpVals[i] / sin(flipAngles[i] * B1);
	}
	linearLeastSquares(X, Y, n, &slope, &inter, &residual);
	double eT1 = exp(-TR / T1);
	*T2 = -TR / log((eT1 - slope) / (1. - slope * eT1));
	double eT2 = exp(-TR / (*T2));
	*M0 = inter * (1. - eT1 * eT2) / (eT2 * (1. - eT1));
	return residual;
}

//******************************************************************************
// Stochastic Region Contraction
// nS = number of points to sample in parameter space at each step
// nR = number of best points to retain and use for contracting the bounds
//******************************************************************************
