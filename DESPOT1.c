/*
 *  DESPOT1.c
 *  MacRI
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright 2011 Tobias Wood. All rights reserved.
 *
 */

#include "math3d.h"
#include "DESPOT1.h"
#include "stdio.h"

int tests()
{
	float testSPGR[2], testFlips[2];
	float spgrTR = 5, irTR = 5, irFlip = radians(5.);
	int nReadout = 128;
	float T1 = 1500, M0 = 10000., B1 = .8;

	fprintf(stdout, "Running DESPOT1 Tests. Parameters are\n");
	fprintf(stdout, "General - M0: %f B1: %f T1: %f\n", M0, B1, T1);
	fprintf(stdout, "SPGR - TR: %f\n",  spgrTR);
	fprintf(stdout, "SPGR-IR - TR: %f nReadout: %d\n", irTR, nReadout);
	testFlips[0] = radians(5.); testFlips[1] = radians(10.);
	testSPGR[0] = SPGR(M0, B1, testFlips[0], T1, spgrTR);
	testSPGR[1] = SPGR(M0, B1, testFlips[1], T1, spgrTR);
	float testIR[2], testTI[2];
	testTI[0] = 350; testTI[1] = 450;
	// The 0.9 is a scale factor in Sean's code
	testIR[0] = IRSPGR(M0, B1, irFlip, T1, testTI[0] * 0.9, irTR, nReadout);
	testIR[1] = IRSPGR(M0, B1, irFlip, T1, testTI[1] * 0.9, irTR, nReadout);
	
	float spgrRes = calcSPGRResiduals(testSPGR, testFlips, 2, spgrTR, B1, T1, M0);
	float irRes   = calcIRSPGRResiduals(testSPGR, testFlips, 2, spgrTR, testIR, testTI, 2, irFlip, irTR, nReadout, B1);
	fprintf(stdout, "Residuals - SPGR = %f, SPGR-IR = %f\n", spgrRes, irRes);
	
	float checkT1, checkM0, checkB1;
	calcDESPOT1(testSPGR, testFlips, 2, spgrTR, B1, &checkT1, &checkM0);
	fprintf(stdout, "DESPOT1 Result - T1 = %f, M0 = %f\n", checkT1, checkM0);
	checkB1 = calcHIFI(testSPGR, testFlips, 2, spgrTR, 
	                   testIR, testTI, 2, irFlip, irTR, nReadout);
	calcDESPOT1(testSPGR, testFlips, 2, spgrTR, checkB1, &checkT1, &checkM0);					   
	fprintf(stdout, "HIFI Result - T1 = %f, M0 = %f\n, B1 = %f\n", checkT1, checkM0, checkB1);
	
	fprintf(stdout, "\nTesting Noisy Data\n");
	testSPGR[0] = testSPGR[0] * 0.92;
	testSPGR[1] = testSPGR[1] * 1.07;
	testIR[0] = testIR[0] * 1.02;
	testIR[1] = testIR[1] * 0.93;
	spgrRes = calcSPGRResiduals(testSPGR, testFlips, 2, spgrTR, B1, T1, M0);
	irRes   = calcIRSPGRResiduals(testSPGR, testFlips, 2, spgrTR, testIR, testTI, 2, irFlip, irTR, nReadout, B1);
	fprintf(stdout, "Residuals - SPGR = %f, SPGR-IR = %f\n", spgrRes, irRes);
	calcDESPOT1(testSPGR, testFlips, 2, spgrTR, B1, &checkT1, &checkM0);
	fprintf(stdout, "DESPOT1 Result - T1 = %f, M0 = %f\n", checkT1, checkM0);
	checkB1 = calcHIFI(testSPGR, testFlips, 2, spgrTR, 
	                   testIR, testTI, 2, irFlip, irTR, nReadout);
	calcDESPOT1(testSPGR, testFlips, 2, spgrTR, checkB1, &checkT1, &checkM0);	
	fprintf(stdout, "HIFI Result - T1 = %f, M0 = %f, B1 = %f\n", checkT1, checkM0, checkB1);
	fprintf(stdout, "\nTests finished.\n");	
		
	if ((spgrRes == 0.) && (irRes == 0.))
		return 1;
	else
		return 0;
}

float SPGR(float M0, float B1, float flipAngle, float T1, float TR)
{
	float e1 = exp(-TR / T1);
	float spgr = M0 * (1. - e1) * sin(flipAngle * B1) /
				      (1. - e1 * cos(flipAngle * B1));
	return spgr;
}

float IRSPGR(float M0, float B1, float flipAngle, float T1, float TI, float TR, float nReadout)
{
	float M0scale = 0.975;
	float irEfficiency = cos(B1 * M_PI) - 1;

	float fullRepTime = TI + (nReadout * TR);
	float eTI = exp(-TI / T1);
	float eFull = exp(-fullRepTime / T1);

	float irspgr = fabs(M0scale * M0 * sin(B1 * flipAngle) *
					       (1. + irEfficiency * eTI + eFull));
	return irspgr;
}

float calcSPGRResiduals(float *spgrVals, float *flipAngles,
					   int n, float TR, float B1,
					   float T1, float M0)
{	
	float cumulativeR = 0.;
	for (int a = 0.; a < n; a++)
	{
		float theory = SPGR(M0, B1, flipAngles[a], T1, TR);
		float diff = theory - spgrVals[a];
		cumulativeR += sqrt(pow(diff, 2.));
	}
	
	return cumulativeR;
}

float calcIRSPGRResiduals(float *spgrVals, float *flipAngles, int nSPGR, float spgrTR,
						 float *irVals, float *TI, int nIR, float irFlipAngle, float irTR, float nReadout,
						 float B1)
{
	
	float T1, M0;
	calcDESPOT1(spgrVals, flipAngles, nSPGR, spgrTR, B1, &T1, &M0);
	
	float guessSPGR[nSPGR];
	float guessIR[nIR];
	
	float spgrR = 0.; float irR = 0.;
	
	if (!isnan(T1))
	{
		for (int i = 0; i < nSPGR; i++)
		{
			if (T1 > 0.)
				guessSPGR[i] = SPGR(M0, B1, flipAngles[i], T1, spgrTR);
			else
				guessSPGR[i] = 0.;
			float diff = guessSPGR[i] - spgrVals[i];
			spgrR += pow(diff, 2.);
		}
		for (int i = 0; i < nIR; i++)
		{
			if (T1 > 0.)
				guessIR[i] = IRSPGR(M0, B1, irFlipAngle, T1, TI[i], irTR, nReadout);
			else
				guessIR[i] = 0.;
			float diff = guessIR[i] - irVals[i];
			irR += pow(diff, 2.);
		}
		float sumR = irR + spgrR;
		return sumR;
	}
	else
		return NAN;
}

void calcDESPOT1(float *spgrVals, float *flipAngles,
				 int n, float TR, float B1,
				 float *T1, float *M0)
{
	float sumX, sumY, sumXX, sumXY;
	sumX = sumY = sumXX = sumXY = 0.;
	for (int i = 0; i < n; i++)
	{
		float x = spgrVals[i] / tan(flipAngles[i] * B1);
		float y = spgrVals[i] / sin(flipAngles[i] * B1);
		
		sumX  += x;
		sumY  += y;
		sumXX += (x*x);
		sumXY += (x*y);
	}
	
	float slope = (n * sumXY - (sumX * sumY)) / (n * sumXX - (sumX * sumX));
	float inter = (sumY - slope * sumX) / n;
	
	*T1 = -TR / log(slope);
	*M0 = inter / (1 - slope);
	
	if ((*T1 < 0.) || (*M0 < 0.))
	{
		*T1 = 0.; *M0 = 0.;
	}
}

void calcDESPOT1d(double *flipAngles, double *spgrVals, int n,
				  double TR, double B1, double *M0, double *T1)
{
	float sumX, sumY, sumXX, sumXY;
	sumX = sumY = sumXX = sumXY = 0.;
	for (int i = 0; i < n; i++)
	{
		float x = spgrVals[i] / tan(flipAngles[i] * B1);
		float y = spgrVals[i] / sin(flipAngles[i] * B1);
		
		sumX  += x;
		sumY  += y;
		sumXX += (x*x);
		sumXY += (x*y);
	}
	
	float slope = (n * sumXY - (sumX * sumY)) / (n * sumXX - (sumX * sumX));
	float inter = (sumY - slope * sumX) / n;
	
	*T1 = -TR / log(slope);
	*M0 = inter / (1 - slope);
}

double SPGRd(double flipAngle, double *p, double *c)
{
	double M0 = p[0], T1 = p[1], B1 = p[2], TR = c[0];
	double e1 = exp(-TR / T1);
	double spgr = M0 * (1. - e1) * sin(flipAngle * B1) /
				      (1. - e1 * cos(flipAngle * B1));
	return spgr;
}

double IRSPGRd(double TI, double *p, double *c)
{
	double M0 = p[0], T1 = p[1], B1 = p[2];
	double flipAngle = c[0], TR = c[1], nReadout = c[2];

	double M0scale = 0.975;
	double irEfficiency = cos(B1 * M_PI) - 1;

	double fullRepTime = TI + (nReadout * TR);
	double eTI = exp(-TI / T1);
	double eFull = exp(-fullRepTime / T1);

	double irspgr = fabs(M0scale * M0 * sin(B1 * flipAngle) *
					       (1. + irEfficiency * eTI + eFull));
	return irspgr;
}


double calcCombined(const double X, double *p, double *c)
{
	int nSPGR = (int)c[1];
	int nIR   = (int)c[5];
	int current = (int)c[6];
	
	if (current < nSPGR)
	{
		c[6] = (double)(current + 1);
		return SPGRd(X, p, c);
	}
	else if (current < nIR + nSPGR)
	{
		c[6] = (double)(current + 1);
		return IRSPGRd(X, p, &(c[2])); // Cut off the SPGR constants
	}
	else
	{
		fprintf(stdout, "Called calcCombined too many times.\n");
		exit(EXIT_FAILURE);
	}
}

void calcHIFId(double *flipAngles, double *spgrVals, int nSPGR, double spgrTR,
			   double *TI, double *irVals, int nIR, double irFlipAngle, double irTR, double nReadout,
			   double *M0, double *T1, double *B1)
{
	// Get an initial estimate from DESPOT1
	*B1 = 1.;
	calcDESPOT1d(flipAngles, spgrVals, nSPGR, spgrTR, *B1, M0, T1);
	// Set up combined parameter, constants etc.
	double pars[3] = {*M0, *T1, *B1};
	// Oh dear lord this just got hacky
	// Last constant is the current iteration
	double constants[7] = {spgrTR, (double)nSPGR,
						   irFlipAngle, irTR, nReadout, (double)nIR, (double)0};
	double dataX[nSPGR + nIR], dataY[nSPGR + nIR];
	for (int i = 0; i < nSPGR; i++)
	{
		dataX[i] = flipAngles[i];
		dataY[i] = spgrVals[i];
	}
	for (int i = 0; i < nIR; i++)
	{
		dataX[i + nSPGR] = TI[i];
		dataY[i + nSPGR] = irVals[i];
	}
	double res;
	levMar(pars, 3, constants, dataX, dataY, nIR + nSPGR, calcCombined, <#eval_type **derivatives#>, <#double *finalResidue#>)
	*M0 = pars[1];
	*T1 = pars[2];
	*B1 = pars[3];
}
	


float calcHIFI(float *spgrVals, float *flipAngles, int nSPGR, float spgrTR,
              float *irVals, float *TI, int nIR, float irFlipAngle, float irTR, float nReadout)
{	
	// Golden Section Search to find B1	
	// From www.mae.wvu.edu/~smirnov/nr/c10-1.pdf
	float R = 0.61803399; // Golden ratio - 1
	float C = 1 - R;
	float precision = 0.001;	
	
	// Set up initial bracket using some guesses
	float B1_0 = 0.3; float B1_3 = 1.8; float B1_1, B1_2;
	float res1 = calcIRSPGRResiduals(spgrVals, flipAngles, nSPGR, spgrTR,
	                                   irVals, TI, nIR, irFlipAngle, irTR, nReadout,
									   B1_0);
	float res2 = calcIRSPGRResiduals(spgrVals, flipAngles, nSPGR, spgrTR,
	                                   irVals, TI, nIR, irFlipAngle, irTR, nReadout,
									   B1_3);
	
	
	if (res1 < res2)
	{
		B1_1 = B1_0 + 0.2;
		B1_2 = B1_1 + C * (B1_3 - B1_1);
	}
	else
	{
		B1_2 = B1_3 - 0.2;
		B1_1 = B1_2 - C * (B1_2 - B1_0);
	}
	
	res1 = calcIRSPGRResiduals(spgrVals, flipAngles, nSPGR, spgrTR,
							   irVals, TI, nIR, irFlipAngle, irTR, nReadout,
							   B1_1);
	res2 = calcIRSPGRResiduals(spgrVals, flipAngles, nSPGR, spgrTR,
							   irVals, TI, nIR, irFlipAngle, irTR, nReadout,
							   B1_2);
	
	while ( fabs(B1_3 - B1_0) > precision * (fabs(B1_1) + fabs(B1_2)))
	{
		if (res2 < res1)
		{
			B1_0 = B1_1; B1_1 = B1_2;
			B1_2 = R * B1_1 + C * B1_3;
			res1 = res2;
			res2 = calcIRSPGRResiduals(spgrVals, flipAngles, nSPGR, spgrTR,
	                                   irVals, TI, nIR, irFlipAngle, irTR, nReadout,
									   B1_2);
		}
		else
		{
			B1_3 = B1_2; B1_2 = B1_1;
			B1_1 = R * B1_2 + C * B1_0;
			res2 = res1;
			res1 = calcIRSPGRResiduals(spgrVals, flipAngles, nSPGR, spgrTR,
	                                   irVals, TI, nIR, irFlipAngle, irTR, nReadout,
									   B1_1);
		}
	}
	
	// Best value for B1
	if (res1 < res2)
		return B1_1;
	else
		return B1_2;
}