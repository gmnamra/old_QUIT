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

	double fullTR = TI + TR;
	double eTI = exp(-TI / T1);
	double eFull = exp(-fullTR / T1);

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
	
	//fprintf(stdout, "MxMy: %f\t%f\n", Mx, My);
	double ssfp = M0 * sqrt(Mx*Mx + My*My);
	return ssfp;
}

double n1cSSFP(double alpha, double *p, double *c)
{
	double T2 = p[0], dO = p[1],
		   TR = c[0], T1 = c[1], B1 = c[2], rfPhase = c[3];
	
	double phase = rfPhase + (dO * TR * 2. * M_PI);
	double A[9] = { -TR / T2,    phase,       0.,
					  -phase, -TR / T2,       0.,
					      0.,       0., -TR / T1 };

	double ca = cos(B1 * alpha), sa = sin(B1 * alpha);
	double R_rf[9] = { 1., 0., 0.,
					   0., ca, sa,
					   0.,-sa, ca};
	double eye[9] = { 1., 0., 0.,
					  0., 1., 0.,
					  0., 0., 1. };
	double temp1[9], temp2[9], all[9];
	matrixExp(A, 3);
	matrixMult(temp1, A, R_rf, 3, 3, 3);
	arraySub(temp1, eye, temp1, 9);
	arraySub(temp2, eye, A, 9);
	// Now multiply everything together
	matrixSolve(all, temp1, temp2, 3, 3);
	double M0[3] = { 0., 0., 1. };
	double Mobs[3];
	matrixMult(Mobs, all, M0, 3, 3, 1);
	double s =  sqrt(Mobs[0]*Mobs[0] + Mobs[1]*Mobs[1]);
	return s;
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
	
	double M0[2] = {f_a, f_b}, S[2];
	double A[4]  = {-(TR/T1_a + TR/tau_a), TR/tau_b,
	                TR/tau_a, -(TR/T1_b + TR/tau_b)};
	double eye[4] = { 1., 0.,
					  0., 1. };
	//arrayScale(A, A, TR, 4);
	//MAT_D(A, 2, 2);
	matrixExp(A, 2);
	double costerm[4];
	arrayScale(costerm, A, cos(B1 * alpha), 4);
	arraySub(costerm, eye, costerm, 4);
	
	double sinterm[4];
	arraySub(sinterm, eye, A, 4);
	arrayScale(sinterm, sinterm, sin(B1 * alpha), 4);
	matrixSolve(A, costerm, sinterm, 2, 2);
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
	double iT2_a = -(1./T2_a + k_a);
	double iT2_b = -(1./T2_b + k_b);
	double iT1_a = -(1./T1_a + k_a);
	double iT1_b = -(1./T1_b + k_b);
	double phase = rfPhase + (dO * TR * 2. * M_PI);
	double A[36] = { iT2_a * TR, k_b * TR,    phase,     0.,    0.,    0.,
					 k_a  * TR,  iT2_b * TR,  0.,     phase,    0.,    0.,
					 -phase,             0.,     iT2_a * TR,  k_b * TR,   0.,    0.,
					 0.,    -phase,     k_a * TR,    iT2_b * TR, 0.,    0.,
					 0.,    0.,     0.,    0.,     iT1_a * TR, k_b * TR,
					 0.,    0.,     0.,    0.,     k_a * TR,   iT1_b * TR };

	double ca = cos(B1 * alpha), sa = sin(B1 * alpha);
	double R_rf[36] = { 1., 0., 0., 0., 0., 0.,
	                    0., 1., 0., 0., 0., 0.,
					    0., 0., ca, 0., sa, 0.,
					    0., 0., 0., ca, 0., sa,
					    0., 0.,-sa, 0., ca, 0.,
					    0., 0., 0.,-sa, 0., ca };
	double eye[36] = { 1., 0., 0., 0., 0., 0.,
	                   0., 1., 0., 0., 0., 0.,
					   0., 0., 1., 0., 0., 0.,
					   0., 0., 0., 1., 0., 0.,
					   0., 0., 0., 0., 1., 0.,
					   0., 0., 0., 0., 0., 1. };
					   
	double temp1[36], temp2[36], temp3[36]; // First bracket
	matrixExp(A, 6);
	matrixMult(temp1, A, R_rf, 6, 6, 6);
	arraySub(temp1, eye, temp1, 36);
	arraySub(temp2, eye, A, 36);
	
	// Now multiply everything together
	matrixSolve(temp3, temp1, temp2, 6, 6);
	double M0[6] = { 0., 0., 0., 0., f_a , f_b };
	double Mobs[6];
	matrixMult(Mobs, temp3, M0, 6, 6, 1);
	double s =  sqrt(pow(Mobs[0] + Mobs[1], 2.) +
					 pow(Mobs[2] + Mobs[3], 2.));
	return s;
}

void a2cSPGR(double *alpha, double *p, double *c, double *signal, size_t nA)
{
	double T1_a = p[0], T1_b = p[1],
		   f_a = p[4], tau_a = p[5],
		   TR = c[0], B1 = c[1];
	double f_b = 1. - f_a;
	double tau_b = f_b * tau_a / f_a;
	
	double M0[2] = {f_a, f_b}, Mobs[2];
	double A[4]  = {-(TR/T1_a + TR/tau_a), TR/tau_b,
	                TR/tau_a, -(TR/T1_b + TR/tau_b)};
	double eye[4] = { 1., 0.,
					  0., 1. };
	
	matrixExp(A, 2);
	double top[4], bottom[4], all[4];
	arraySub(top, eye, A, 4);
	for (size_t n = 0; n < nA; n++)
	{
		arrayScale(bottom, A, cos(B1 * alpha[n]), 4);
		arraySub(bottom, eye, bottom, 4);
	
		matrixSolve(all, bottom, top, 2, 2);
		arrayScale(all, all, sin(B1 * alpha[n]), 4);
		matrixMult(Mobs, all, M0, 2, 2, 1);
		signal[n] = Mobs[0] + Mobs[1];
	}
}

void a2cSSFP(double *alpha, double *p, double *c, double *signal, size_t nA)
{
	double T1_a = p[0], T1_b = p[1], T2_a = p[2], T2_b = p[3],
		   f_a = p[4], tau_a = p[5], dO = p[6],
		   TR = c[0], B1 = c[1], rfPhase = c[2];
	double f_b = 1. - f_a;
	double tau_b = f_b * tau_a / f_a;
	double eT2_a = -TR * (1./T2_a + 1./tau_a);
	double eT2_b = -TR * (1./T2_b + 1./tau_b);
	double eT1_a = -TR * (1./T1_a + 1./tau_a);
	double eT1_b = -TR * (1./T1_b + 1./tau_b);
	double k_a   = TR / tau_a;
	double k_b   = TR / tau_b;
	double phase = rfPhase + (dO * TR * 2. * M_PI);

	double eye[36] = { 1., 0., 0., 0., 0., 0.,
	                   0., 1., 0., 0., 0., 0.,
					   0., 0., 1., 0., 0., 0.,
					   0., 0., 0., 1., 0., 0.,
					   0., 0., 0., 0., 1., 0.,
					   0., 0., 0., 0., 0., 1. };
	double A[36] = {  eT2_a,    k_b, phase,    0.,    0.,    0.,
					    k_a,  eT2_b,    0., phase,    0.,    0.,
					 -phase,     0., eT2_a,   k_b,    0.,    0.,
					     0., -phase,   k_a, eT2_b,    0.,    0.,
					     0.,     0.,    0.,    0., eT1_a,   k_b,
					     0.,     0.,    0.,    0.,   k_a, eT1_b };
	double R_rf[36] = { 1., 0., 0., 0., 0., 0.,
	                    0., 1., 0., 0., 0., 0.,
					    0., 0., 0., 0., 0., 0.,
					    0., 0., 0., 0., 0., 0.,
					    0., 0., 0., 0., 0., 0.,
					    0., 0., 0., 0., 0., 0. };

	double top[36], bottom[36], all[36];
	double M0[6] = { 0., 0., 0., 0., f_a , f_b }, Mobs[6];
	
	matrixExp(A, 6);
	// Top of matrix divide
	arraySub(top, eye, A, 36);

	for (size_t n = 0; n < nA; n++)
	{
		double ca = cos(B1 * alpha[n]), sa = sin(B1 * alpha[n]);
		R_rf[14] = ca; R_rf[21] = ca; R_rf[28] = ca; R_rf[35] = ca;
		R_rf[16] = sa; R_rf[23] = sa; R_rf[26] = -sa; R_rf[33] = -sa;
		// Inverse bracket term (i.e. bottom of matrix divide)
		matrixMult(bottom, A, R_rf, 6, 6, 6);
		arraySub(bottom, eye, bottom, 36);
		
		// Matrix 'divide'
		matrixSolve(all, bottom, top, 6, 6);
		matrixMult(Mobs, all, M0, 6, 6, 1);
		signal[n] =  sqrt(pow(Mobs[0] + Mobs[1], 2.) +
					      pow(Mobs[2] + Mobs[3], 2.));
	}
}

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
*/
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
	
	matrixExp(A, 3);
	double top[9], bottom[9], all[9];
	arraySub(top, eye, A, 9);
	for (size_t n = 0; n < nA; n++)
	{
		arrayScale(bottom, A, cos(B1 * alpha[n]), 9);
		arraySub(bottom, eye, bottom, 9);
	
		matrixSolve(all, bottom, top, 3, 3);
		arrayScale(all, all, sin(B1 * alpha[n]), 9);
		matrixMult(Mobs, all, M0, 3, 3, 1);
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
	
	matrixExp(A, 9);
	// Top of matrix divide
	arraySub(top, eye, A, 81);

	for (size_t n = 0; n < nA; n++)
	{
		double ca = cos(B1 * alpha[n]), sa = sin(B1 * alpha[n]);
		R[30] = ca; R[40] = ca; R[50] = ca; R[60] = ca; R[70] = ca; R[80] = ca;
		R[33] = sa; R[43] = sa; R[53] = sa; R[57] = -sa; R[67] = -sa; R[77] = -sa;
		// Inverse bracket term (i.e. bottom of matrix divide)
		matrixMult(bottom, A, R, 9, 9, 9);
		arraySub(bottom, eye, bottom, 81);
		
		// Matrix 'divide'
		matrixSolve(all, bottom, top, 9, 9);
		matrixMult(Mobs, all, M0, 9, 9, 1);
		signal[n] =  sqrt(pow(Mobs[0] + Mobs[1] + Mobs[2], 2.) +
					      pow(Mobs[3] + Mobs[4] + Mobs[5], 2.));
	}
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
	/*double loBounds[2] = {1., 150. / TR}; // A T2 of 0. causes nasty maths
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
					  bounds, 200, 20, 100, 0.005, 0., &fRes);
	*M0 = 1.; *T2 = p[0]; *dO = p[1];*/
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
