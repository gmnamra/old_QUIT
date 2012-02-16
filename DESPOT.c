/*
 *  DESPOT1.c
 *  MacRI
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright 2011 Tobias Wood. All rights reserved.
 *
 */

#include "DESPOT.h"
#include "math3d.h"
#include "stdio.h"
#include "cblas.h"
#include "clapack.h"

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

void SPGR_Jacobian(double *angles, int nD, double *p, double *c, double *result)
{
	double M0 = p[0], T1 = p[1], B1 = p[2], TR = c[0];
	double eTR = exp(-TR / T1);
	for (int d = 0; d < nD; d++)
	{
		double alpha = angles[d];
		
		double denom = (1. - eTR * cos(B1 * alpha));
		
		double dMzM0 = (1 - eTR * sin(B1 * alpha)) / denom;
		double dMzT1 = (M0 * TR * sin(B1 * alpha) * eTR * (cos(B1 * alpha) - 1.)) /
		               (T1 * T1 * denom * denom);
		double dMzB1 = (M0 * B1 * eTR * (1. - eTR + cos(B1 * alpha))) / (denom * denom);
		result[0 * nD + d] = dMzM0;
		result[1 * nD + d] = dMzT1;
		result[2 * nD + d] = dMzB1;
	}
}

double IRSPGR(double TI, double *p, double *c)
{
	double M0 = p[0], T1 = p[1], B1 = p[2];
	double flipAngle = c[0], TR = c[1];
	
	double irEfficiency = cos(B1 * M_PI) - 1;

	double fullRepTime = TI + TR;
	double eTI = exp(-TI / T1);
	double eFull = exp(-fullRepTime / T1);

	double irspgr = fabs(M0 * sin(B1 * flipAngle) *
					       (1. + irEfficiency * eTI + eFull));
	return irspgr;
}

void IRSPGR_Jacobian(double *data, int nD, double *p, double *c, double *result)
{
	double M0 = p[0], T1 = p[1], B1 = p[2];
	double alpha = c[0], TR = c[1], nReadout = c[2];
	
	for (int d = 0; d < nD; d++)
	{
		double TI = data[d];
		double irEff = cos(B1 * M_PI) - 1;
		
		double fullTR = TI + (nReadout * TR);
		double eTI = exp(-TI / T1);
		double eTR = exp(-fullTR / T1);
		
		double dMzM0 = sin(B1 * alpha) * (1. + eTR + irEff * eTI);
		double dMzT1 = (M0 * sin(B1 * alpha) / (T1 * T1)) *
					   (fullTR * eTR + TI * irEff * eTI);
		double b1 = M0 * alpha * cos(B1 * alpha) *
					   (1 + eTR + irEff * eTI);
		double b2 =    M0 * sin(B1 * alpha) *
					   (M_PI * sin(B1 * M_PI) * eTI);
		double dMzB1 = b1 - b2;
		result[0 * nD + d] = dMzM0;
		result[1 * nD + d] = dMzT1;
		result[2 * nD + d] = dMzB1;
	}
}

double SSFP(double flipAngle, double *p, double *c)
{
	double M0 = p[0], T2 = p[1], dO = p[2];
	double TR = c[0], T1 = c[1], B1 = c[2], offset = c[3];
	
	double eT1 = exp(-TR / T1);
	double eT2 = exp(-TR / T2);
	
	double phase = offset + dO * (TR / 1.e3) * 2. * M_PI;
	double sina = sin(B1 * flipAngle);
	double cosa = cos(B1 * flipAngle);
	double sinp = sin(phase);
	double cosp = cos(phase);
	
	double denom = ((1. - eT1 * cosa) * (1. - eT2 * cosp)) - 
				   (eT2 * (eT1 - cosa) * (eT2 - cosp));
	
	double Mx = ((1 - eT1) * eT2 * sina * (cosp - eT2)) / denom;
	double My = ((1.- eT1) * eT2 * sina * sinp) / denom;
	double ssfp = M0 * sqrt(Mx*Mx + My*My);
	return ssfp;
}

double nSSFP(double flipAngle, double *p, double *c)
{
	double M0 = 1., T2 = p[0], dO = p[1];
	double TR = c[0], T1 = c[1], B1 = c[2], offset = c[3];
	
	double eT1 = exp(-TR / T1);
	double eT2 = exp(-TR / T2);
	
	double phase = offset + dO * (TR / 1.e3) * 2. * M_PI;
	double sina = sin(B1 * flipAngle);
	double cosa = cos(B1 * flipAngle);
	double sinp = sin(phase);
	double cosp = cos(phase);
	
	double denom = ((1. - eT1 * cosa) * (1. - eT2 * cosp)) - 
				   (eT2 * (eT1 - cosa) * (eT2 - cosp));
	
	double Mx = ((1 - eT1) * eT2 * sina * (cosp - eT2)) / denom;
	double My = ((1.- eT1) * eT2 * sina * sinp) / denom;
	double ssfp = M0 * sqrt(Mx*Mx + My*My);
	return ssfp;
}

// Normalised, 2 component versions
/*	Full parameter vector is
	0 - T1_s
	1 - T1_f
	2 - T2_s
	3 - T2_f
	4 - f_s
	5 - tau_s
	6 - dw
	Relationships:
	f_s + f_f = 1. (Only two components)
	k_ = 1. / tau_ (Exchange is inverse of lifetime)
	f_s / tau_s = f_f / tau_f (Exchange equilibrium)
	Constants vector:
	0 - TR
	1 - B1
	2 - Phase cycle/offset
*/
double n2cSPGR(double alpha, double *p, double *c)
{
	double T1_s = p[0], T1_f = p[1],
		   f_s = p[4], tau_s = p[5],
		   TR = c[0], B1 = c[1];
	double f_f = 1. - f_s;
	double tau_f = f_f * tau_s / f_s;
	double k_s = 1. / tau_s, k_f = 1. / tau_f;
	
	double M0[2] = {f_s, f_f}, S[2];
	double A[4]  = {(-1./T1_f - k_f), k_s,
	                k_f, (-1./T1_s - k_s)};
	double eye[4] = { 1., 0.,
					  0., 1. };
	arrayExp(A, A, TR, 4);
	double sinterm[4];
	arraySub(sinterm, eye, A, 4);
	arrayScale(sinterm, sinterm, sin(B1 * alpha), 4);
	double costerm[4];
	arrayScale(costerm, A, cos(B1 * alpha), 4);
	arraySub(costerm, eye, costerm, 4);
	// Inverse
	int ipiv[2];
	clapack_dgetrf(CblasRowMajor, 2, 2, costerm, 2, ipiv);
	clapack_dgetri(CblasRowMajor, 2, costerm, 2, ipiv);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 1., sinterm, 2, costerm, 2, 0., A, 2);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 2, 2, 1., A, 2, M0, 1, 0., S, 1);
	double s = sqrt(S[0] * S[0] + S[1] * S[1]);
	return s;
}

double n2cSSFP(double alpha, double *p, double *c)
{
	double T1_s = p[0], T1_f = p[1], T2_s = p[2], T2_f = p[3],
		   f_s = p[4], tau_s = p[5], dw = p[6],
		   TR = c[0], B1 = c[1], rfPhase = c[2];
	double f_f = 1. - f_s;
	double tau_f = f_f * tau_s / f_s;
	double k_s = 1. / tau_s, k_f = 1. / tau_f;
	
	double iT2_f = -1./T2_f - k_f;
	double iT2_s = -1./T2_s - k_s;
	double iT1_f = -1./T1_f - k_f;
	double iT1_s = -1./T1_s - k_s;
	double phase = rfPhase + dw;

	int ipiv[6]; // General for all cblas ops
	double eye[36] = { 1., 0., 0., 0., 0., 0.,
	                   0., 1., 0., 0., 0., 0.,
					   0., 0., 1., 0., 0., 0.,
					   0., 0., 0., 1., 0., 0.,
					   0., 0., 0., 0., 1., 0.,
					   0., 0., 0., 0., 0., 1. };
	
	double A[36] = { iT2_f,  k_s,    phase, 0.,    0.,    0.,
					 k_f,    iT2_s,  0.,    phase, 0.,    0.,
					 -phase, 0.,     iT2_f, k_s,   0.,    0.,
					 0.,     -phase, k_f,   iT2_s, 0.,    0.,
					 0.,     0.,     0.,    0.,    iT1_f, k_s,
					 0.,     0.,     0.,    0.,    k_f,   iT1_s };
	double invA[36]; arrayCopy(invA, A, 36);
	clapack_dgetrf(CblasRowMajor, 6, 6, invA, 6, ipiv); // Inverse
	clapack_dgetri(CblasRowMajor, 6, invA, 6, ipiv);
	arrayExp(A, A, TR, 36);
	
	double ca = cos(B1 * alpha), sa = sin(B1 * alpha);
	double R[36] = { 1., 0., 0., 0., 0., 0.,
	                 0., 1., 0., 0., 0., 0.,
					 0., 0., ca, 0., sa, 0.,
					 0., 0., 0., ca, 0., sa,
					 0., 0.,-sa, 0., ca, 0.,
					 0., 0., 0.,-sa, 0., ca };
					 
	double temp1[36], temp2[36]; // First bracket
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6, 6, 6, 1., A, 6, R, 6, 0., temp1, 6);
	//arraySub(temp1, eye, temp1, 36);
	catlas_daxpby(36, 1., eye, 1, -1., temp1, 1);
	clapack_dgetrf(CblasRowMajor, 6, 6, temp1, 6, ipiv); // Inverse
	clapack_dgetri(CblasRowMajor, 6, temp1, 6, ipiv);	
	
	//arraySub(A, A, eye, 36); // Second bracket
	catlas_daxpby(36, -1., eye, 1, 1., A, 1);
	// Now multiply everything together
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6, 6, 6, 1., temp1, 6, A, 6, 0., temp2, 6);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6, 6, 6, 1., temp2, 6, invA, 6, 0., temp1, 6);
	
	double initC[6] = { 0., 0., 0., 0., f_f / (T1_f * tau_f), f_s / (T1_s * tau_s) };
	double finalC[6];
	
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 6, 1., temp1, 6, initC, 1, 0., finalC, 1);
	
	double s =  sqrt(finalC[0] * finalC[0] +
	                 finalC[1] * finalC[1] +
				     finalC[2] * finalC[2] +
				     finalC[3] * finalC[3]);
	return s;
}

void calcDESPOT1(double *flipAngles, double *spgrVals, int n,
				 double TR, double B1, double *M0, double *T1)
{
	// Linearise the data, then least-squares
	double X[n], Y[n], slope, inter;
	for (int i = 0; i < n; i++)
	{
		X[i] = spgrVals[i] / tan(flipAngles[i] * B1);
		Y[i] = spgrVals[i] / sin(flipAngles[i] * B1);
	}
	linearLeastSquares(X, Y, n, &slope, &inter);	
	*T1 = -TR / log(slope);
	*M0 = inter / (1. - slope);
}

void classicDESPOT2(double *flipAngles, double *ssfpVals, int n,
                    double TR, double T1, double B1, double *M0, double *T2)
{
	// As above, linearise, then least-squares
	double X[n], Y[n], slope, inter;
	for (int i = 0; i < n; i++)
	{
		X[i] = ssfpVals[i] / tan(flipAngles[i] * B1);
		Y[i] = ssfpVals[i] / sin(flipAngles[i] * B1);
	}
	linearLeastSquares(X, Y, n, &slope, &inter);
	double eT1 = exp(-TR / T1);
	*T2 = -TR / log((eT1 - slope) / (1. - slope * eT1));
	double eT2 = exp(-TR / (*T2));
	*M0 = inter * (1. - eT1 * eT2) / (1. - eT1);
}

void simplexDESPOT2(size_t nPhases, size_t *nD, double *phases, double **flipAngles, double **ssfp,
					double TR, double T1, double B1, double *M0, double *T2, double *dO)
{
	// Gather together all the data
	double p[3] = {*M0, *T2, *dO}, *c[nPhases], fRes = 0.;
	eval_type *f[nPhases];
	for (int i = 0; i < nPhases; i++)
	{
		c[i] = malloc(4 * sizeof(double));
		c[i][0] = TR; c[i][1] = T1; c[i][2] = B1;
		c[i][3] = phases[i];
		f[i] = SSFP;
	}
	
	simplex(p, 3, c, nPhases, flipAngles, ssfp, nD, f, NULL, &fRes);
	*M0 = p[0]; *T2 = p[1]; *dO = p[2];
}

void contractDESPOT2(size_t nPhases, size_t *nD, double *phases, double **flipAngles, double **ssfp,
					 double TR, double T1, double B1, double *M0, double *T2, double *dO)
{
	double loBounds[2] = {1., 150. / TR}; // A T2 of 0. causes nasty maths
	double hiBounds[2] = {300., 100. / TR};
	double *bounds[2] = {loBounds, hiBounds}; 
	double p[2] = {*T2, *dO}, *c[nPhases], fRes = 0.;
	eval_type *f[nPhases];
	for (int i = 0; i < nPhases; i++)
	{
		c[i] = malloc(4 * sizeof(double));
		c[i][0] = TR; c[i][1] = T1; c[i][2] = B1;
		c[i][3] = phases[i];
		f[i] = nSSFP;
		arrayScale(ssfp[i], ssfp[i], 1. / arrayMean(ssfp[i], nD[i]), nD[i]);
	}
	regionContraction(p, 2, c, nPhases, flipAngles, ssfp, nD, true, f,
					  bounds, 200, 20, 0.005, &fRes);
	*M0 = 1.; *T2 = p[0]; *dO = p[1];
}

/*
	Params:
	T1_s, T1_f, T2_s, T2_f,	f_s, tau_s, dw
	Consts:
	TR, B1, rfPhase
*/
void mcDESPOT(size_t nSPGR, double *spgrAlpha, double *spgr, double spgrTR,
			  size_t nPhases, size_t *nSSFPs, double *phases, double **ssfpAlphas, double **ssfp,
              double ssfpTR, double T1, double B1, double *p)
{
	double loBounds[7] = { 0.1, 0.1, 0.001, 0.001,   0., 0.025, 0. };
	double hiBounds[7] = { 2.5, 2.5, 0.150, 0.150, 0.45, 0.500,  1./ssfpTR };
	double *bounds[2] =  { loBounds, hiBounds };

	size_t nD[1 + nPhases];
	eval_type *f[1 + nPhases];
	double *alphas[1 + nPhases], *data[1 + nPhases];
	double *c[1 + nPhases], fRes = 0.;
	
	nD[0]     = nSPGR;
	f[0]      = n2cSPGR;
	alphas[0] = spgrAlpha;
	data[0]   = spgr;
	c[0] = malloc(3 * sizeof(double));
	c[0][0] = spgrTR; c[0][1] = B1; c[0][2] = 0.;
	arrayScale(spgr, spgr, 1. / arrayMean(spgr, nSPGR), nSPGR);
	for (int i = 0; i < nPhases; i++)
	{
		nD[i + 1]     = nSSFPs[i];
		f[i + 1]      = n2cSSFP;
		alphas[i + 1] = ssfpAlphas[i];
		data[i + 1]   = ssfp[i];
		c[i + 1] = malloc(3 * sizeof(double));
		c[i + 1][0] = ssfpTR; c[i + 1][1] = B1; c[i + 1][2] = phases[i];
		arrayScale(ssfp[i], ssfp[i], 1. / arrayMean(ssfp[i], nD[i]), nD[i]);
	}
	regionContraction(p, 7, c, 1 + nPhases, alphas, data, nD, true, f,
					  bounds, 2000, 50, 0.005, &fRes);
	for (int i = 0; i < 1 + nPhases; i++)
		free(c[i]);
}

double calcSPGR(double *angles, double *spgrVals, int n, double TR,
                double *M0, double *T1, double *B1)
{
	double par[3] = {*M0, *T1, *B1};
	double res = 0;
	levMar(par, 3, &TR, angles, spgrVals, n, SPGR, SPGR_Jacobian, &res);
	*M0 = par[0]; *T1 = par[1]; *B1 = par[2];
	return res;
}

double calcIR(double *TI, double *irVals, int nIR,
              double alpha, double TR,
			  double *M0, double *T1, double *B1)
{
	double par[3] = {*M0, *T1, *B1};
	double con[2] = {alpha, TR};
	double res = 0;
	levMar(par, 3, con, TI, irVals, nIR, IRSPGR, IRSPGR_Jacobian, &res);
	*M0 = par[0]; *T1 = par[1]; *B1 = par[2];
	return res;
}

double calcHIFI(double *flipAngles, double *spgrVals, int nSPGR, double spgrTR,
				double *TI, double *irVals, int nIR, double irFlipAngle, double irTR,
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
	double irConstants[2] = { irFlipAngle, irTR };
	double spgrRes[nSPGR], irRes[nIR];
	
	par[2] = B1_0;
	calcDESPOT1(flipAngles, spgrVals, nSPGR, spgrTR, par[2], &(par[0]), &(par[1]));
	double res1 = calcResiduals(par, spgrConstants, flipAngles, spgrVals, nSPGR, &SPGR, spgrRes, false) +
	              calcResiduals(par, irConstants, TI, irVals, nIR, &IRSPGR, irRes, false);
	par[2] = B1_3;
	calcDESPOT1(flipAngles, spgrVals, nSPGR, spgrTR, par[2], &(par[0]), &(par[1]));
	double res2 = calcResiduals(par, spgrConstants, flipAngles, spgrVals, nSPGR, SPGR, spgrRes, false) +
	              calcResiduals(par, irConstants, TI, irVals, nIR, IRSPGR, irRes, false);
	
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
	res1 = calcResiduals(par, spgrConstants, flipAngles, spgrVals, nSPGR, SPGR, spgrRes, false) +
	       calcResiduals(par, irConstants, TI, irVals, nIR, IRSPGR, irRes, false);
	par[2] = B1_2;
	calcDESPOT1(flipAngles, spgrVals, nSPGR, spgrTR, par[2], &(par[0]), &(par[1]));
	res2 = calcResiduals(par, spgrConstants, flipAngles, spgrVals, nSPGR, SPGR, spgrRes, false) +
	       calcResiduals(par, irConstants, TI, irVals, nIR, IRSPGR, irRes, false);
	
	while ( fabs(B1_3 - B1_0) > precision * (fabs(B1_1) + fabs(B1_2)))
	{
		if (res2 < res1)
		{
			B1_0 = B1_1; B1_1 = B1_2;
			B1_2 = R * B1_1 + C * B1_3;
			res1 = res2;
			par[2] = B1_2;
			calcDESPOT1(flipAngles, spgrVals, nSPGR, spgrTR, par[2], &(par[0]), &(par[1]));
			res2 = calcResiduals(par, spgrConstants, flipAngles, spgrVals, nSPGR, SPGR, spgrRes, false) +
	               calcResiduals(par, irConstants, TI, irVals, nIR, IRSPGR, irRes, false);
		}
		else
		{
			B1_3 = B1_2; B1_2 = B1_1;
			B1_1 = R * B1_2 + C * B1_0;
			res2 = res1;
			par[2] = B1_1;
			calcDESPOT1(flipAngles, spgrVals, nSPGR, spgrTR, par[2], &(par[0]), &(par[1]));
			res1 = calcResiduals(par, spgrConstants, flipAngles, spgrVals, nSPGR, SPGR, spgrRes, false) +
	               calcResiduals(par, irConstants, TI, irVals, nIR, IRSPGR, irRes, false);
		}
	}
	
	// Best value for B1
	*M0 = par[0]; *T1 = par[1];
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