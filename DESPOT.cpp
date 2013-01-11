/*
 *  DESPOT.cpp
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
				      double TR, double B1, double &M0, double &T1)
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
	T1 = -TR / log(slope);
	M0 = inter / (1. - slope);
	return res;
}

double classicDESPOT2(const ArrayXd &flipAngles, const ArrayXd &ssfpVals,
                      double TR, double T1, double B1, double &M0, double &T2) {
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
	T2 = -TR / log((eT1 - slope) / (1. - slope * eT1));
	double eT2 = exp(-TR / T2);
	M0 = inter * (1. - eT1 * eT2) / (eT2 * (1. - eT1));
	return residual;
}

ArrayXd SPGR(const ArrayXd &flip, const double &B1, const double &TR, const double &M0, const double &T1)
{
	double e1 = exp(-TR / T1);
	ArrayXd spgr = M0 * (1. - e1) * sin(B1 * flip) / (1. - (e1 * cos(B1 * flip)));
	return spgr;
}

ArrayXd IRSPGR(const ArrayXd &TI, const double &B1, const double &TR,
               const double &flip, const double &eff,
			   const double &M0, const double &T1)
{
	double irEfficiency;
	// Adiabatic pulses aren't susceptible to B1 problems.
	if (eff > 0)
		irEfficiency = cos(eff * M_PI) - 1;
	else
		irEfficiency = cos(B1 * M_PI) - 1;
	
	ArrayXd eTI = exp(-TI / T1);
	ArrayXd eFull = exp(-(TI + TR) / T1);

	VectorXd irspgr = ((M0 * sin(B1 * flip) * (1. + irEfficiency * eTI + eFull))).abs();
	return irspgr;
}

double calcHIFI(const ArrayXd &flipAngles, const ArrayXd &spgrVals, const double spgrTR,
				const ArrayXd &TI, const ArrayXd &irVals, const double irFlipAngle,
				const double irTR, const double nReadout, const double eff,
                double &M0, double &T1, double &B1) {
	// Golden Section Search to find B1
	// From www.mae.wvu.edu/~smirnov/nr/c10-1.pdf
	double R = 0.61803399; // Golden ratio - 1
	double C = 1 - R;
	double precision = 0.001;
    
	// Set up initial bracket using some guesses
	double B1_0 = 0.3; double B1_3 = 1.8; double B1_1, B1_2;
	
	B1 = B1_0;
	classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
	double res1 = (spgrVals - SPGR(flipAngles, spgrTR, B1, M0, T1)).square().sum() +
	              (irVals - IRSPGR(TI, irTR, B1, irFlipAngle, eff, M0, T1)).square().sum();
	B1 = B1_3;
	classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
	double res2 = (spgrVals - SPGR(flipAngles, spgrTR, B1, M0, T1)).square().sum() +
	              (irVals - IRSPGR(TI, irTR, B1, irFlipAngle, eff, M0, T1)).square().sum();
    
	if (res1 < res2) {
		B1_1 = B1_0 + 0.2;
		B1_2 = B1_1 + C * (B1_3 - B1_1);
	} else {
		B1_2 = B1_3 - 0.2;
		B1_1 = B1_2 - C * (B1_2 - B1_0);
	}
    
	B1 = B1_1;
	classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
	res1 = (spgrVals - SPGR(flipAngles, spgrTR, B1, M0, T1)).square().sum() +
	       (irVals - IRSPGR(TI, irTR, B1, irFlipAngle, eff, M0, T1)).square().sum();
	B1 = B1_2;
	classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
	res2 = (spgrVals - SPGR(flipAngles, spgrTR, B1, M0, T1)).square().sum() +
	       (irVals - IRSPGR(TI, irTR, B1, irFlipAngle, eff, M0, T1)).square().sum();
	while ( fabs(B1_3 - B1_0) > precision * (fabs(B1_1) + fabs(B1_2))) {
		if (res2 < res1) {
			B1_0 = B1_1; B1_1 = B1_2;
			B1_2 = R * B1_1 + C * B1_3;
			res1 = res2;
			B1 = B1_2;
			classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
			res2 = (spgrVals - SPGR(flipAngles, spgrTR, B1, M0, T1)).square().sum() +
			       (irVals - IRSPGR(TI, irTR, B1, irFlipAngle, eff, M0, T1)).square().sum();
		} else {
			B1_3 = B1_2; B1_2 = B1_1;
			B1_1 = R * B1_2 + C * B1_0;
			res2 = res1;
			B1 = B1_1;
			classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
			res1 = (spgrVals - SPGR(flipAngles, spgrTR, B1, M0, T1)).square().sum() +
			       (irVals - IRSPGR(TI, irTR, B1, irFlipAngle, eff, M0, T1)).square().sum();
		}
	}
	// Best value for B1
	if (res1 < res2) {
		B1 = B1_1;
		return res1;
	} else {
		B1 = B1_2;
		return res2;
	}
}