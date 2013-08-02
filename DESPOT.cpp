/*
 *  DESPOT.cpp
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright (c) 2011-2013 Tobias Wood.
 *
 *  Based in part on work by Sean Deoni
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
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

//******************************************************************************
// Basic least squares fitting
//******************************************************************************
void linearLeastSquares(const ArrayXd &X, const ArrayXd &Y,
                        double &slope, double &inter, double &res)
{
	eigen_assert(X.size() == Y.size());
	double sumX, sumY, sumXX, sumXY;
	sumX = sumY = sumXX = sumXY = 0.;
	for (int i = 0; i < X.size(); i++) {
		double x = X[i];
		double y = Y[i];
		
		sumX  += x;
		sumY  += y;
		sumXX += (x*x);
		sumXY += (x*y);
	}
	
	slope = (sumXY - (sumX * sumY) / X.size()) /
	        (sumXX - (sumX * sumX) / X.size());
	inter = (sumY - (slope) * sumX) / X.size();
	
	res = (Y - (X*slope + inter)).square().sum();
}

double classicDESPOT1(const ArrayXd &flipAngles, const ArrayXd &spgrVals,
				      double TR, double B1, double &M0, double &T1)
{
	// Linearise the data, then least-squares
	long n = flipAngles.size();
	ArrayXd X(n), Y(n);
	double slope, inter, res;
	for (long i = 0; i < flipAngles.size(); i++) {
		X[i] = spgrVals[i] / tan(flipAngles[i] * B1);
		Y[i] = spgrVals[i] / sin(flipAngles[i] * B1);
	}
	linearLeastSquares(X, Y, slope, inter, res);	
	T1 = -TR / log(slope);
	M0 = inter / (1. - slope);
	return res;
}

double classicDESPOT2(const ArrayXd &flipAngles, const ArrayXd &ssfpVals,
                      double TR, double T1, double B1, double &M0, double &T2) {
	// As above, linearise, then least-squares
	long n = flipAngles.size();
	ArrayXd X(n), Y(n);
	double slope, inter, residual;
	for (long i = 0; i < flipAngles.size(); i++) {
		X[i] = ssfpVals[i] / tan(flipAngles[i] * B1);
		Y[i] = ssfpVals[i] / sin(flipAngles[i] * B1);
	}
	linearLeastSquares(X, Y, slope, inter, residual);
	double eT1 = exp(-TR / T1);
	T2 = -TR / log((eT1 - slope) / (1. - slope * eT1));
	double eT2 = exp(-TR / T2);
	M0 = inter * (1. - eT1 * eT2) / (eT2 * (1. - eT1));
	return residual;
}

ArrayXd SPGR(const ArrayXd &flip, const double &TR, const double &B1, const double &M0, const double &T1)
{
	double e1 = exp(-TR / T1);
	ArrayXd spgr = M0 * (1. - e1) * sin(B1 * flip) / (1. - (e1 * cos(B1 * flip)));
	return spgr;
}

ArrayXd IRSPGR(const ArrayXd &TI, const double &TR, const double &B1,
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

double HIFIResidual(const ArrayXd &flipAngles, const ArrayXd &spgrVals, const double spgrTR,
				const ArrayXd &TI, const ArrayXd &irVals, const double irFlipAngle,
				const double irTR, const double nReadout, const double eff,
                double &M0, double &T1, double &B1) {
	ArrayXd st = SPGR(flipAngles, spgrTR, B1, M0, T1);
	ArrayXd it = IRSPGR(TI, irTR, B1, irFlipAngle, eff, M0, T1);
	double res = (spgrVals - st).square().sum() +
	              (irVals - it).square().sum();
	return res;
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
	double res1 = HIFIResidual(flipAngles, spgrVals, spgrTR, TI, irVals, irFlipAngle, irTR, nReadout, eff, M0, T1, B1);
	B1 = B1_3;
	classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
	double res2 = HIFIResidual(flipAngles, spgrVals, spgrTR, TI, irVals, irFlipAngle, irTR, nReadout, eff, M0, T1, B1);
	if (res1 < res2) {
		B1_1 = B1_0 + 0.2;
		B1_2 = B1_1 + C * (B1_3 - B1_1);
	} else {
		B1_2 = B1_3 - 0.2;
		B1_1 = B1_2 - C * (B1_2 - B1_0);
	}
	
	B1 = B1_1;
	classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
	res1 = HIFIResidual(flipAngles, spgrVals, spgrTR, TI, irVals, irFlipAngle, irTR, nReadout, eff, M0, T1, B1);
	B1 = B1_2;
	classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
	res2 = HIFIResidual(flipAngles, spgrVals, spgrTR, TI, irVals, irFlipAngle, irTR, nReadout, eff, M0, T1, B1);
	while ( fabs(B1_3 - B1_0) > precision * (fabs(B1_1) + fabs(B1_2))) {
		if (res2 < res1) {
			B1_0 = B1_1; B1_1 = B1_2;
			B1_2 = R * B1_1 + C * B1_3;
			res1 = res2;
			B1 = B1_2;
			classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
			res2 = HIFIResidual(flipAngles, spgrVals, spgrTR, TI, irVals, irFlipAngle, irTR, nReadout, eff, M0, T1, B1);
		} else {
			B1_3 = B1_2; B1_2 = B1_1;
			B1_1 = R * B1_2 + C * B1_0;
			res2 = res1;
			B1 = B1_1;
			classicDESPOT1(flipAngles, spgrVals, spgrTR, B1, M0, T1);
			res1 = HIFIResidual(flipAngles, spgrVals, spgrTR, TI, irVals, irFlipAngle, irTR, nReadout, eff, M0, T1, B1);
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
	std::cout << "Finished Golden Ratio Search" << std::endl<< std::endl<< std::endl<< std::endl;
}
