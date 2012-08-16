/*
 *  DESPOT1.c
 *  MacRI
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright 2011 Tobias Wood. All rights reserved.
 *
 */

#include "DESPOT.h"

// Normalised, 3 component versions
/*	Full parameter vector is
	0 - T1_a
	1 - T1_b
	2 - T1_c
	3 - T2_a
	4 - T2_b
	5 - T2_c
	6 - f_a
	7 - f_b
	8 - tau_a
	9 - tau_b
	10 - dw
	Relationships:
	f_b + f_a = 1. (Only two components)
	k_ = 1. / tau_ (Exchange is inverse of lifetime)
	f_b / tau_b = f_a / tau_a (Exchange equilibrium)
	Constants vector:
	0 - TR
	1 - B1
	2 - Phase cycle/offset
void a3cSPGR(double *alpha, double *p, double *c, double *signal, size_t nA)
{
	double T1_a = p[0], T1_b = p[1], T1_c = p[2],
		   f_a = p[6], f_b = p[7],
		   tau_a = p[8], tau_b = p[9],
		   TR = c[0], B1 = c[1];
	double f_c = 1. - f_a - f_b;
	double tau_c = f_c * tau_b / f_a;
	
	double M0[3] = {f_a, f_b, f_c}, Mobs[3];
	double eT1_a = -TR * (1. / T1_a + 1 / tau_a);
	double eT1_b = -TR * (1. / T1_b + 1 / tau_b);
	double eT1_c = -TR * (1. / T1_c + 1 / tau_c);
	double k_a = TR / tau_a;
	double k_b = TR / tau_b;
	double k_c = TR / tau_c;
	double A[9]  = { eT1_a,   k_b,   k_c,
	                   k_a, eT1_b,   k_c,
					   k_a,   k_b, eT1_c }; 
	double eye[9] = { 1., 0., 0.,
					  0., 1., 0.,
					  0., 0., 1. };
	
	matrix_exp(A, 3);
	double top[9], bottom[9], all[9];
	arraySub(top, eye, A, 9);
	for (size_t n = 0; n < nA; n++)
	{
		arrayScale(bottom, A, cos(B1 * alpha[n]), 9);
		arraySub(bottom, eye, bottom, 9);
	
		matrixSolve(all, bottom, top, 3, 3);
		arrayScale(all, all, sin(B1 * alpha[n]), 9);
		matrix_mult(Mobs, all, M0, 3, 3, 1);
		signal[n] = Mobs[0] + Mobs[1] + Mobs[2];
	}
}

void a3cSSFP(double *alpha, double *p, double *c, double *signal, size_t nA)
{
	double T1_a = p[0], T1_b = p[1], T1_c = p[2],
	       T2_a = p[3], T2_b = p[4], T2_c = p[5],
		   f_a = p[6], f_b = p[7],
		   tau_a = p[8], tau_b = p[9], dO = p[10],
		   TR = c[0], B1 = c[1], rfPhase = c[2];
		   
	double f_c = 1. - f_a - f_b;
	double tau_c = f_b * tau_a / f_a;
	double eT2_a = -TR * (1./T2_a + 1./tau_a);
	double eT2_b = -TR * (1./T2_b + 1./tau_b);
	double eT2_c = -TR * (1./T2_c + 1./tau_c);	
	double eT1_a = -TR * (1./T1_a + 1./tau_a);
	double eT1_b = -TR * (1./T1_b + 1./tau_b);
	double eT1_c = -TR * (1./T1_c + 1./tau_c);
	double k_a   = TR / tau_a;
	double k_b   = TR / tau_b;
	double k_c   = TR / tau_c;
	double phase = rfPhase + (dO * TR * 2. * M_PI);

	double eye[81] = { 1., 0., 0., 0., 0., 0., 0., 0., 0.,
	                   0., 1., 0., 0., 0., 0., 0., 0., 0.,
					   0., 0., 1., 0., 0., 0., 0., 0., 0.,
					   0., 0., 0., 1., 0., 0., 0., 0., 0.,
					   0., 0., 0., 0., 1., 0., 0., 0., 0.,
					   0., 0., 0., 0., 0., 1., 0., 0., 0.,
					   0., 0., 0., 0., 0., 0., 1., 0., 0.,
					   0., 0., 0., 0., 0., 0., 0., 1., 0.,
					   0., 0., 0., 0., 0., 0., 0., 0., 1. };
	double A[81] = {  eT2_a,    k_b,    k_c, phase,    0.,    0.,    0.,    0.,    0.,
					    k_a,  eT2_b,    k_c,    0., phase,    0.,    0.,    0.,    0.,
					    k_a,    k_b,  eT2_c,    0.,    0., phase,    0.,    0.,    0.,
					 -phase,     0.,     0., eT2_a,   k_b,   k_c,    0.,    0.,    0.,
					     0., -phase,     0.,   k_a, eT2_b,   k_c,    0.,    0.,    0.,
					     0.,     0., -phase,   k_a,   k_b, eT2_c,    0.,    0.,    0.,
					     0.,     0.,     0.,    0.,    0.,    0., eT1_a,   k_b,   k_c,
					     0.,     0.,     0.,    0.,    0.,    0.,   k_a, eT1_b,   k_c,
					     0.,     0.,     0.,    0.,    0.,    0.,   k_a,   k_b, eT1_c };
	double R[81] = { 1., 0., 0., 0., 0., 0., 0., 0., 0.,
				     0., 1., 0., 0., 0., 0., 0., 0., 0.,
					 0., 0., 1., 0., 0., 0., 0., 0., 0.,
					 0., 0., 0., 0., 0., 0., 0., 0., 0.,
					 0., 0., 0., 0., 0., 0., 0., 0., 0.,
					 0., 0., 0., 0., 0., 0., 0., 0., 0.,
					 0., 0., 0., 0., 0., 0., 0., 0., 0.,
					 0., 0., 0., 0., 0., 0., 0., 0., 0.,
					 0., 0., 0., 0., 0., 0., 0., 0., 0. };

	double top[81], bottom[81], all[81];
	double M0[9] = { 0., 0., 0., 0., 0., 0., f_a , f_b, f_c }, Mobs[9];
	
	matrix_exp(A, 9);
	// Top of matrix divide
	arraySub(top, eye, A, 81);

	for (size_t n = 0; n < nA; n++)
	{
		double ca = cos(B1 * alpha[n]), sa = sin(B1 * alpha[n]);
		R[30] = ca; R[40] = ca; R[50] = ca; R[60] = ca; R[70] = ca; R[80] = ca;
		R[33] = sa; R[43] = sa; R[53] = sa; R[57] = -sa; R[67] = -sa; R[77] = -sa;
		// Inverse bracket term (i.e. bottom of matrix divide)
		matrix_mult(bottom, A, R, 9, 9, 9);
		arraySub(bottom, eye, bottom, 81);
		
		// Matrix 'divide'
		matrixSolve(all, bottom, top, 9, 9);
		matrix_mult(Mobs, all, M0, 9, 9, 1);
		signal[n] =  sqrt(pow(Mobs[0] + Mobs[1] + Mobs[2], 2.) +
					      pow(Mobs[3] + Mobs[4] + Mobs[5], 2.));
	}
}*/

double radians(double degrees)
{	return degrees * M_PI / 180.;	}
double degrees(double radians)
{	return radians * 180. / M_PI;	}

double clamp(double value, double low, double high)
{
	if (value < low)
		return low;
	if (value > high)
		return high;
	return value;
}

//******************************************************************************
#pragma mark Console Input
//******************************************************************************
int fgetArray(FILE *in, char type, size_t n, void *array)
{
		for (size_t i = 0; i < n; i++)
		{
			switch (type)
			{
				case 'i':
				{	int inVal; fscanf(in, "%d", &inVal);
					((int *)array)[i] = inVal;
				} break;
				case 'f':
				{	float inVal; fscanf(in, "%f", &inVal);
					((float *)array)[i] = inVal;
				} break;
				case 'd':
				{	double inVal; fscanf(in, "%lf", &inVal);
					((double *)array)[i] = inVal;
				} break;
				case 's':
				{	char inVal[1024]; fscanf(stdin, "%s", &inVal);
					((char **)array)[i] = (char *)malloc(strlen(inVal) * sizeof(char));
					strcpy(((char **)array)[i], inVal); } break;
			}
		}
	return n;
}

//******************************************************************************
// Basic least squares fitting
//******************************************************************************
void linearLeastSquares(double *X, double *Y, int nD,
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
	
	*slope = (nD * sumXY - (sumX * sumY)) / (nD * sumXX - (sumX * sumX));
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
	int n = flipAngles.size();
	double X[n], Y[n], slope, inter, res;
	for (int i = 0; i < flipAngles.size(); i++)
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
	// p[0] = M0, p[1] = T2
	int n = flipAngles.size();
	double X[n], Y[n], slope, inter, residual;
	for (int i = 0; i < flipAngles.size(); i++)
	{
		X[i] = ssfpVals[i] / tan(flipAngles[i] * B1);
		Y[i] = ssfpVals[i] / sin(flipAngles[i] * B1);
	}
	linearLeastSquares(X, Y, n, &slope, &inter, &residual);
	double eT1 = exp(-TR / T1);
	*T2 = -TR / log((eT1 - slope) / (1. - slope * eT1));
	double eT2 = exp(-TR / *T2);
	*M0 = inter * (1. - eT1 * eT2) / (1. - eT1);
	return residual;
}
