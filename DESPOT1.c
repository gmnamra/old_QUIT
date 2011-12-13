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

/*int tests()
{
	double testSPGR[2], testFlips[2];
	double spgrTR = 5, irTR = 5, irFlip = radians(5.);
	int nReadout = 128;
	double T1 = 1500, M0 = 10000., B1 = .8;

	fprintf(stdout, "Running DESPOT1 Tests. Parameters are\n");
	fprintf(stdout, "General - M0: %f B1: %f T1: %f\n", M0, B1, T1);
	fprintf(stdout, "SPGR - TR: %f\n",  spgrTR);
	fprintf(stdout, "SPGR-IR - TR: %f nReadout: %d\n", irTR, nReadout);
	testFlips[0] = radians(5.); testFlips[1] = radians(10.);
	testSPGR[0] = SPGR(M0, B1, testFlips[0], T1, spgrTR);
	testSPGR[1] = SPGR(M0, B1, testFlips[1], T1, spgrTR);
	double testIR[2], testTI[2];
	testTI[0] = 350; testTI[1] = 450;
	// The 0.9 is a scale factor in Sean's code
	testIR[0] = IRSPGR(M0, B1, irFlip, T1, testTI[0] * 0.9, irTR, nReadout);
	testIR[1] = IRSPGR(M0, B1, irFlip, T1, testTI[1] * 0.9, irTR, nReadout);
	
	double spgrRes = calcSPGRResiduals(testSPGR, testFlips, 2, spgrTR, B1, T1, M0);
	double irRes   = calcIRSPGRResiduals(testSPGR, testFlips, 2, spgrTR, testIR, testTI, 2, irFlip, irTR, nReadout, B1);
	fprintf(stdout, "Residuals - SPGR = %f, SPGR-IR = %f\n", spgrRes, irRes);
	
	double checkT1, checkM0, checkB1;
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
}*/

double SPGR(double flipAngle, double *p, double *c)
{
	double M0 = p[0], T1 = p[1], B1 = p[2], TR = c[0];
	double e1 = exp(-TR / T1);
	double spgr = M0 * (1. - e1) * sin(flipAngle * B1) /
				      (1. - e1 * cos(flipAngle * B1));
	return spgr;
}

double IRSPGR(double TI, double *p, double *c)
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

void calcDESPOT1(double *flipAngles, double *spgrVals, int n,
				 double TR, double B1, double *M0, double *T1)
{
	double sumX, sumY, sumXX, sumXY;
	sumX = sumY = sumXX = sumXY = 0.;
	for (int i = 0; i < n; i++)
	{
		double x = spgrVals[i] / tan(flipAngles[i] * B1);
		double y = spgrVals[i] / sin(flipAngles[i] * B1);
		
		sumX  += x;
		sumY  += y;
		sumXX += (x*x);
		sumXY += (x*y);
	}
	
	double slope = (n * sumXY - (sumX * sumY)) / (n * sumXX - (sumX * sumX));
	double inter = (sumY - slope * sumX) / n;
	
	*T1 = -TR / log(slope);
	*M0 = inter / (1 - slope);
}

double calcHIFI(double *flipAngles, double *spgrVals, int nSPGR, double spgrTR,
				double *TI, double *irVals, int nIR, double irFlipAngle, double irTR, double nReadout,
				double *M0, double *T1, double *B1)
{
	// Golden Section Search to find B1	
	// From www.mae.wvu.edu/~smirnov/nr/c10-1.pdf
	double R = 0.61803399; // Golden ratio - 1
	double C = 1 - R;
	double precision = 0.001;	
	
	// Set up initial bracket using some guesses
	double B1_0 = 0.3; double B1_3 = 1.8; double B1_1, B1_2;
	
	// Assemble parameters
	double par[3] = { *M0, *T1, B1_1 };
	double spgrConstants[1] = { spgrTR };
	double irConstants[3] = { irFlipAngle, irTR, nReadout };
	double spgrRes[nSPGR], irRes[nIR];
	
	par[2] = B1_0;
	calcDESPOT1(flipAngles, spgrVals, nSPGR, spgrTR, par[2], &(par[0]), &(par[1]));
	double res1 = calcResiduals(par, spgrConstants, flipAngles, spgrVals, nSPGR, &SPGR, spgrRes) +
	              calcResiduals(par, irConstants, TI, irVals, nIR, &IRSPGR, irRes);
	par[2] = B1_3;
	calcDESPOT1(flipAngles, spgrVals, nSPGR, spgrTR, par[2], &(par[0]), &(par[1]));
	double res2 = calcResiduals(par, spgrConstants, flipAngles, spgrVals, nSPGR, SPGR, spgrRes) +
	              calcResiduals(par, irConstants, TI, irVals, nIR, IRSPGR, irRes);
	
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
	
	par[2] = B1_1;
	calcDESPOT1(flipAngles, spgrVals, nSPGR, spgrTR, par[2], &(par[0]), &(par[1]));
	res1 = calcResiduals(par, spgrConstants, flipAngles, spgrVals, nSPGR, SPGR, spgrRes) +
	       calcResiduals(par, irConstants, TI, irVals, nIR, IRSPGR, irRes);
	par[2] = B1_2;
	calcDESPOT1(flipAngles, spgrVals, nSPGR, spgrTR, par[2], &(par[0]), &(par[1]));
	res2 = calcResiduals(par, spgrConstants, flipAngles, spgrVals, nSPGR, SPGR, spgrRes) +
	       calcResiduals(par, irConstants, TI, irVals, nIR, IRSPGR, irRes);
	
	while ( fabs(B1_3 - B1_0) > precision * (fabs(B1_1) + fabs(B1_2)))
	{
		if (res2 < res1)
		{
			B1_0 = B1_1; B1_1 = B1_2;
			B1_2 = R * B1_1 + C * B1_3;
			res1 = res2;
			par[2] = B1_2;
			calcDESPOT1(flipAngles, spgrVals, nSPGR, spgrTR, par[2], &(par[0]), &(par[1]));
			res2 = calcResiduals(par, spgrConstants, flipAngles, spgrVals, nSPGR, SPGR, spgrRes) +
	               calcResiduals(par, irConstants, TI, irVals, nIR, IRSPGR, irRes);
		}
		else
		{
			B1_3 = B1_2; B1_2 = B1_1;
			B1_1 = R * B1_2 + C * B1_0;
			res2 = res1;
			par[2] = B1_1;
			calcDESPOT1(flipAngles, spgrVals, nSPGR, spgrTR, par[2], &(par[0]), &(par[1]));
			res1 = calcResiduals(par, spgrConstants, flipAngles, spgrVals, nSPGR, SPGR, spgrRes) +
	               calcResiduals(par, irConstants, TI, irVals, nIR, IRSPGR, irRes);
		}
	}
	
	// Best value for B1
	if (res1 < res2)
	{
		*B1 = B1_1;
		return res1;
	}
	else
	{
		*B1 = B1_2;
		return res2;
	}
}