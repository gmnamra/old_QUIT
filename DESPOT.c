/*
 *  DESPOT1.c
 *  MacRI
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright 2011 Tobias Wood. All rights reserved.
 *
 */

#include "DESPOT.h"
#include "mathMatrix.h"
#include "math2d.h"
#include "stdio.h"
#include "cblas.h"
#include "clapack.h"

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
	
	double phase = offset + dO * TR * 2. * M_PI;
	double sina = sin(B1 * flipAngle);
	double cosa = cos(B1 * flipAngle);
	double sinp = sin(phase);
	double cosp = cos(phase);
	
	double denom = ((1. - eT1 * cosa) * (1. - eT2 * cosp)) - 
				   (eT2 * (eT1 - cosa) * (eT2 - cosp));
	
	double My = ((1 - eT1) * eT2 * sina * (cosp - eT2)) / denom;
	double Mx = ((1.- eT1) * eT2 * sina * sinp) / denom;
	double ssfp = M0 * sqrt(Mx*Mx + My*My);
	return ssfp;
}

// Normalised, 2 component versions
/*	Full parameter vector is
	0 - T1_a
	1 - T1_b
	2 - T2_a
	3 - T2_b
	4 - f_a
	5 - tau_a
	6 - dw
	Relationships:
	f_b + f_a = 1. (Only two components)
	k_ = 1. / tau_ (Exchange is inverse of lifetime)
	f_b / tau_b = f_a / tau_a (Exchange equilibrium)
	Constants vector:
	0 - TR
	1 - B1
	2 - Phase cycle/offset
*/
double n2cSPGR(double alpha, double *p, double *c)
{
	double T1_a = p[0], T1_b = p[1],
		   f_a = p[4], tau_a = p[5],
		   TR = c[0], B1 = c[1];
	double f_b = 1. - f_a;
	double tau_b = f_b * tau_a / f_a;
	double k_b = 1. / tau_b, k_a = 1. / tau_a;
	
	double M0[2] = {f_a, f_b}, S[2];
	double A[4]  = {-(1./T1_a + k_a), k_b,
	                k_a, -(1./T1_b + k_b)};
	double eye[4] = { 1., 0.,
					  0., 1. };
	arrayScale(A, A, TR, 4);
	matrixExp(A, 2);
	
	double costerm[4];
	arrayScale(costerm, A, cos(B1 * alpha), 4);
	arraySub(costerm, eye, costerm, 4);
	// Inverse
	double invcterm[4];
	inv22(invcterm, costerm);
	
	double sinterm[4];
	arraySub(sinterm, eye, A, 4);
	arrayScale(sinterm, sinterm, sin(B1 * alpha), 4);
	matrixMult(A, invcterm, sinterm, 2, 2, 2);
	matrixMult(S, A, M0, 2, 2, 1);
	double s = S[0] + S[1];
	return s;
}

double n2cSSFP(double alpha, double *p, double *c)
{
	double T1_a = p[0], T1_b = p[1], T2_a = p[2], T2_b = p[3],
		   f_a = p[4], tau_a = p[5], dO = p[6],
		   TR = c[0], B1 = c[1], rfPhase = c[2];
	double f_b = 1. - f_a;
	double tau_b = f_b * tau_a / f_a;
	double k_a = 1. / tau_a, k_b = 1. / tau_b;
	
	double iT2_a = TR * -(1./T2_a + k_a);
	double iT2_b = TR * -(1./T2_b + k_b);
	double iT1_a = TR * -(1./T1_a + k_a);
	double iT1_b = TR * -(1./T1_b + k_b);

	k_b *= TR; k_a *= TR;
	double A[36] = { iT2_a, k_b,    0.,     0.,    0.,    0.,
					 k_a,   iT2_b,  0.,     0.,    0.,    0.,
					 0.,    0.,     iT2_a,  k_b,   0.,    0.,
					 0.,    0.,     k_a,    iT2_b, 0.,    0.,
					 0.,    0.,     0.,    0.,     iT1_a, k_b,
					 0.,    0.,     0.,    0.,     k_a,   iT1_b };

	double ca = cos(B1 * alpha), sa = sin(B1 * alpha);
	double R_rf[36] = { 1., 0., 0., 0., 0., 0.,
	                    0., 1., 0., 0., 0., 0.,
					    0., 0., ca, 0., sa, 0.,
					    0., 0., 0., ca, 0., sa,
					    0., 0.,-sa, 0., ca, 0.,
					    0., 0., 0.,-sa, 0., ca };
	
	double phase = rfPhase + (dO * TR * 2. * M_PI);
	double cp = cos(phase), sp = sin(phase);
	double R_ph[36] = { cp,  0., sp, 0., 0., 0.,
						0., cp, 0., sp, 0., 0.,
					   -sp,  0., cp, 0., 0., 0.,
					    0.,-sp, 0., cp, 0., 0.,
					    0., 0., 0., 0., 1., 0.,
					    0., 0., 0., 0., 0., 1. };
	double eye[36] = { 1., 0., 0., 0., 0., 0.,
	                   0., 1., 0., 0., 0., 0.,
					   0., 0., 1., 0., 0., 0.,
					   0., 0., 0., 1., 0., 0.,
					   0., 0., 0., 0., 1., 0.,
					   0., 0., 0., 0., 0., 1. };
					   
	double R[36], temp1[36], temp2[36], temp3[36]; // First bracket
	matrixExp(A, 6);
	matrixMult(R, R_ph, R_rf, 6, 6, 6);
	matrixMult(temp1, A, R, 6, 6, 6);
	arraySub(temp1, eye, temp1, 36);
	int ipiv[6]; // General for all cblas ops
	clapack_dgetrf(CblasRowMajor, 6, 6, temp1, 6, ipiv); // Inverse
	clapack_dgetri(CblasRowMajor, 6, temp1, 6, ipiv);	
	
	arraySub(temp2, A, eye, 36);

	// Now multiply everything together
	matrixMult(temp3, temp1, temp2, 6, 6, 6);
	double initC[6] = { 0., 0., 0., 0., f_a , f_b };
	double finalC[6];
	matrixMult(finalC, temp3, initC, 6, 6, 1);
	double s =  sqrt(pow(finalC[0] + finalC[1], 2.) +
					 pow(finalC[2] + finalC[3], 2.));
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
					  bounds, 200, 20, 100, 0.005, &fRes);
	*M0 = 1.; *T2 = p[0]; *dO = p[1];
}

/*
	Params:
	T1_m, T1_f, T2_m, T2_f,	f_m, tau_m, dw
	Consts:
	TR, B1, rfPhase
*/
void mcDESPOT(size_t nSPGR, double *spgrAlpha, double *spgr, double spgrTR,
			  size_t nPhases, size_t *nSSFPs, double *phases, double **ssfpAlphas, double **ssfp,
              double ssfpTR, double T1, double B1, double *p)
{
	double loBounds[7] = {  100.,  100.,   1.,   1.,  0.,  25.,        0. };
	double hiBounds[7] = { 2500., 2500., 150., 150.,  1., 500., 1./ssfpTR };
	double *bounds[2] =  { loBounds, hiBounds };

	size_t nD[1 + nPhases];
	eval_type *f[1 + nPhases];
	double *alphas[1 + nPhases], *data[1 + nPhases];
	double *c[1 + nPhases];
	
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
	size_t ctracts = regionContraction(p, 7, c, 1 + nPhases, alphas, data, nD, true, f,
					                   bounds, 5000, 50, 100, 0.005, &(p[7]));

	// Assume that the short T2 component is the Myelin
	// Hence swap parameters if p[2] is larger
	fprintf(stdout, "Finished after %ld contractions. ", ctracts);
	if (p[2] > p[3])
	{
		fprintf(stdout, "Swapped.\n");	
		double temp = p[2];
		p[2] = p[3]; p[3] = temp;
		temp = p[0]; p[0] = p[1]; p[1] = temp;
		p[4] = 1. - p[4];
		p[5] = p[4] * p[5] / (1. - p[4]);
	}
	else
		fprintf(stdout, "No swap.\n");
	ARR_D(bounds[0], 7);
	ARR_D(bounds[1], 7);
	ARR_D(p, 8);
	
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