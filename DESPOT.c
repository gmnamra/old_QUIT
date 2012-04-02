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

float SPGR(float flipAngle, float *p, float *c)
{
	float M0 = p[0], T1 = p[1], B1 = p[2], TR = c[0];
	float e1 = exp(-TR / T1);
	float spgr = M0 * (1. - e1) * sin(flipAngle * B1) /
				      (1. - e1 * cos(flipAngle * B1));
	return spgr;
}

void SPGR_Jacobian(float *angles, int nD, float *p, float *c, float *result)
{
	float M0 = p[0], T1 = p[1], B1 = p[2], TR = c[0];
	float eTR = exp(-TR / T1);
	for (int d = 0; d < nD; d++)
	{
		float alpha = angles[d];
		
		float denom = (1. - eTR * cos(B1 * alpha));
		
		float dMzM0 = (1 - eTR * sin(B1 * alpha)) / denom;
		float dMzT1 = (M0 * TR * sin(B1 * alpha) * eTR * (cos(B1 * alpha) - 1.)) /
		               (T1 * T1 * denom * denom);
		float dMzB1 = (M0 * B1 * eTR * (1. - eTR + cos(B1 * alpha))) / (denom * denom);
		result[0 * nD + d] = dMzM0;
		result[1 * nD + d] = dMzT1;
		result[2 * nD + d] = dMzB1;
	}
}

float IRSPGR(float TI, float *p, float *c)
{
	float M0 = p[0], T1 = p[1], B1 = p[2];
	float flipAngle = c[0], TR = c[1];
	
	float irEfficiency = cos(B1 * M_PI) - 1;

	float fullTR = TI + TR;
	float eTI = exp(-TI / T1);
	float eFull = exp(-fullTR / T1);

	float irspgr = fabs(M0 * sin(B1 * flipAngle) *
					       (1. + irEfficiency * eTI + eFull));
	return irspgr;
}

void IRSPGR_Jacobian(float *data, int nD, float *p, float *c, float *result)
{
	float M0 = p[0], T1 = p[1], B1 = p[2];
	float alpha = c[0], TR = c[1], nReadout = c[2];
	
	for (int d = 0; d < nD; d++)
	{
		float TI = data[d];
		float irEff = cos(B1 * M_PI) - 1;
		
		float fullTR = TI + (nReadout * TR);
		float eTI = exp(-TI / T1);
		float eTR = exp(-fullTR / T1);
		
		float dMzM0 = sin(B1 * alpha) * (1. + eTR + irEff * eTI);
		float dMzT1 = (M0 * sin(B1 * alpha) / (T1 * T1)) *
					   (fullTR * eTR + TI * irEff * eTI);
		float b1 = M0 * alpha * cos(B1 * alpha) *
					   (1 + eTR + irEff * eTI);
		float b2 =    M0 * sin(B1 * alpha) *
					   (M_PI * sin(B1 * M_PI) * eTI);
		float dMzB1 = b1 - b2;
		result[0 * nD + d] = dMzM0;
		result[1 * nD + d] = dMzT1;
		result[2 * nD + d] = dMzB1;
	}
}

float SSFP(float flipAngle, float *p, float *c)
{
	float M0 = p[0], T2 = p[1], B0 = p[2];
	float TR = c[0], T1 = c[1], B1 = c[2], offset = c[3];
	
	float eT1 = exp(-TR / T1);
	float eT2 = exp(-TR / T2);
	
	float phase = offset + B0 * (TR / 1.e3) * 2. * M_PI;
	float sina = sin(B1 * flipAngle);
	float cosa = cos(B1 * flipAngle);
	float sinp = sin(phase);
	float cosp = cos(phase);
	
	float denom = ((1. - eT1 * cosa) * (1. - eT2 * cosp)) - 
				   (eT2 * (eT1 - cosa) * (eT2 - cosp));
	
	float Mx = ((1 - eT1) * eT2 * sina * (cosp - eT2)) / denom;
	float My = ((1.- eT1) * eT2 * sina * sinp) / denom;
	float ssfp = M0 * sqrt(Mx*Mx + My*My);
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
float n2cSPGR(float alpha, float *p, float *c)
{
	float T1_a = p[0], T1_b = p[1],
		   f_a = p[4], tau_a = p[5],
		   TR = c[0], B1 = c[1];
	float f_b = 1. - f_a;
	float tau_b = f_b * tau_a / f_a;
	
	float M0[2] = {f_a, f_b}, S[2];
	float A[4]  = {-(TR/T1_a + TR/tau_a), TR/tau_b,
	                TR/tau_a, -(TR/T1_b + TR/tau_b)};
	float eye[4] = { 1., 0.,
					  0., 1. };
	matrixExp(A, 2);
	float costerm[4];
	arrayScale(costerm, A, cos(B1 * alpha), 4);
	arraySub(costerm, eye, costerm, 4);
	
	float sinterm[4];
	arraySub(sinterm, eye, A, 4);
	arrayScale(sinterm, sinterm, sin(B1 * alpha), 4);
	matrixSolve(A, costerm, sinterm, 2, 2);
	matrixMult(S, A, M0, 2, 2, 1);
	float s = S[0] + S[1];
	return s;
}

float n2cSSFP(float alpha, float *p, float *c)
{
	float T1_a = p[0], T1_b = p[1], T2_a = p[2], T2_b = p[3],
		   f_a = p[4], tau_a = p[5], dO = p[6],
		   TR = c[0], B1 = c[1], rfPhase = c[2];
	float f_b = 1. - f_a;
	float tau_b = f_b * tau_a / f_a;
	float k_a = 1. / tau_a, k_b = 1. / tau_b;
	float iT2_a = -(1./T2_a + k_a);
	float iT2_b = -(1./T2_b + k_b);
	float iT1_a = -(1./T1_a + k_a);
	float iT1_b = -(1./T1_b + k_b);
	float phase = rfPhase + (dO * TR * 2. * M_PI);
	float A[36] = { iT2_a * TR, k_b * TR,    phase,     0.,    0.,    0.,
					 k_a  * TR,  iT2_b * TR,  0.,     phase,    0.,    0.,
					 -phase,             0.,     iT2_a * TR,  k_b * TR,   0.,    0.,
					 0.,    -phase,     k_a * TR,    iT2_b * TR, 0.,    0.,
					 0.,    0.,     0.,    0.,     iT1_a * TR, k_b * TR,
					 0.,    0.,     0.,    0.,     k_a * TR,   iT1_b * TR };

	float ca = cos(B1 * alpha), sa = sin(B1 * alpha);
	float R_rf[36] = { 1., 0., 0., 0., 0., 0.,
	                    0., 1., 0., 0., 0., 0.,
					    0., 0., ca, 0., sa, 0.,
					    0., 0., 0., ca, 0., sa,
					    0., 0.,-sa, 0., ca, 0.,
					    0., 0., 0.,-sa, 0., ca };
	float eye[36] = { 1., 0., 0., 0., 0., 0.,
	                   0., 1., 0., 0., 0., 0.,
					   0., 0., 1., 0., 0., 0.,
					   0., 0., 0., 1., 0., 0.,
					   0., 0., 0., 0., 1., 0.,
					   0., 0., 0., 0., 0., 1. };
					   
	float temp1[36], temp2[36], temp3[36]; // First bracket
	matrixExp(A, 6);
	matrixMult(temp1, A, R_rf, 6, 6, 6);
	arraySub(temp1, eye, temp1, 36);
	arraySub(temp2, eye, A, 36);
	
	// Now multiply everything together
	matrixSolve(temp3, temp1, temp2, 6, 6);
	float M0[6] = { 0., 0., 0., 0., f_a , f_b };
	float Mobs[6];
	matrixMult(Mobs, temp3, M0, 6, 6, 1);
	float s =  sqrt(pow(Mobs[0] + Mobs[1], 2.) +
					 pow(Mobs[2] + Mobs[3], 2.));
	return s;
}

void a1cSSFP(float *alpha, float *p, float *c, float *signal, size_t nA)
{
	float T2 = p[0], B0 = p[1],
		   TR = c[0], T1 = c[1], B1 = c[2], rfPhase = c[3];
	
	float phase = rfPhase + (B0 * TR * 2. * M_PI);
	float A[9] = { -TR / T2,    phase,       0.,
					  -phase, -TR / T2,       0.,
					      0.,       0., -TR / T1 };

	float R_rf[9] = { 1., 0., 0.,
					   0., 0., 0.,
					   0., 0., 0.};
	float eye[9] = { 1., 0., 0.,
					  0., 1., 0.,
					  0., 0., 1. };
	float temp1[9], temp2[9], all[9];
	matrixExp(A, 3);
	arraySub(temp2, eye, A, 9);
	
	for (size_t n = 0; n < nA; n++)
	{
		float ca = cos(B1 * alpha[n]), sa = sin(B1 * alpha[n]);
		R_rf[4] = R_rf[8] = ca;
		R_rf[5] = sa;
		R_rf[7] = -sa;
		matrixMult(temp1, A, R_rf, 3, 3, 3);
		arraySub(temp1, eye, temp1, 9);
		matrixSolve(all, temp1, temp2, 3, 3);
		float M0[3] = { 0., 0., 1. };
		float Mobs[3];
		matrixMult(Mobs, all, M0, 3, 3, 1);
		signal[n] = sqrt(Mobs[0]*Mobs[0] + Mobs[1]*Mobs[1]);
	}
}

void a2cSPGR(float *alpha, float *p, float *c, float *signal, size_t nA)
{
	float T1_a = p[0], T1_b = p[1],
		   f_a = p[4], tau_a = p[5],
		   TR = c[0], B1 = c[1];
	float f_b = 1. - f_a;
	float tau_b = f_b * tau_a / f_a;
	
	float M0[2] = {f_a, f_b}, Mobs[2];
	float A[4]  = {-(TR/T1_a + TR/tau_a), TR/tau_b,
	                TR/tau_a, -(TR/T1_b + TR/tau_b)};
	float eye[4] = { 1., 0.,
					  0., 1. };
	
	matrixExp(A, 2);
	float top[4], bottom[4], all[4];
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

void a2cSSFP(float *alpha, float *p, float *c, float *signal, size_t nA)
{
	float T1_a = p[0], T1_b = p[1], T2_a = p[2], T2_b = p[3],
		   f_a = p[4], tau_a = p[5], dO = p[6],
		   TR = c[0], B1 = c[1], rfPhase = c[2];
	float f_b = 1. - f_a;
	float tau_b = f_b * tau_a / f_a;
	float eT2_a = -TR * (1./T2_a + 1./tau_a);
	float eT2_b = -TR * (1./T2_b + 1./tau_b);
	float eT1_a = -TR * (1./T1_a + 1./tau_a);
	float eT1_b = -TR * (1./T1_b + 1./tau_b);
	float k_a   = TR / tau_a;
	float k_b   = TR / tau_b;
	float phase = rfPhase + (dO * TR * 2. * M_PI);

	float eye[36] = { 1., 0., 0., 0., 0., 0.,
	                   0., 1., 0., 0., 0., 0.,
					   0., 0., 1., 0., 0., 0.,
					   0., 0., 0., 1., 0., 0.,
					   0., 0., 0., 0., 1., 0.,
					   0., 0., 0., 0., 0., 1. };
	float A[36] = {  eT2_a,    k_b, phase,    0.,    0.,    0.,
					    k_a,  eT2_b,    0., phase,    0.,    0.,
					 -phase,     0., eT2_a,   k_b,    0.,    0.,
					     0., -phase,   k_a, eT2_b,    0.,    0.,
					     0.,     0.,    0.,    0., eT1_a,   k_b,
					     0.,     0.,    0.,    0.,   k_a, eT1_b };
	float R_rf[36] = { 1., 0., 0., 0., 0., 0.,
	                    0., 1., 0., 0., 0., 0.,
					    0., 0., 0., 0., 0., 0.,
					    0., 0., 0., 0., 0., 0.,
					    0., 0., 0., 0., 0., 0.,
					    0., 0., 0., 0., 0., 0. };

	float top[36], bottom[36], all[36];
	float M0[6] = { 0., 0., 0., 0., f_a , f_b }, Mobs[6];
	
	matrixExp(A, 6);
	// Top of matrix divide
	arraySub(top, eye, A, 36);

	for (size_t n = 0; n < nA; n++)
	{
		float ca = cos(B1 * alpha[n]), sa = sin(B1 * alpha[n]);
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
void a3cSPGR(float *alpha, float *p, float *c, float *signal, size_t nA)
{
	float T1_a = p[0], T1_b = p[1], T1_c = p[2],
		   f_a = p[6], f_b = p[7],
		   tau_a = p[8], tau_b = p[9],
		   TR = c[0], B1 = c[1];
	float f_c = 1. - f_a - f_b;
	float tau_c = f_c * tau_b / f_a;
	
	float M0[3] = {f_a, f_b, f_c}, Mobs[3];
	float eT1_a = -TR * (1. / T1_a + 1 / tau_a);
	float eT1_b = -TR * (1. / T1_b + 1 / tau_b);
	float eT1_c = -TR * (1. / T1_c + 1 / tau_c);
	float k_a = TR / tau_a;
	float k_b = TR / tau_b;
	float k_c = TR / tau_c;
	float A[9]  = { eT1_a,   k_b,   k_c,
	                   k_a, eT1_b,   k_c,
					   k_a,   k_b, eT1_c }; 
	float eye[9] = { 1., 0., 0.,
					  0., 1., 0.,
					  0., 0., 1. };
	
	matrixExp(A, 3);
	float top[9], bottom[9], all[9];
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

void a3cSSFP(float *alpha, float *p, float *c, float *signal, size_t nA)
{
	float T1_a = p[0], T1_b = p[1], T1_c = p[2],
	       T2_a = p[3], T2_b = p[4], T2_c = p[5],
		   f_a = p[6], f_b = p[7],
		   tau_a = p[8], tau_b = p[9], dO = p[10],
		   TR = c[0], B1 = c[1], rfPhase = c[2];
		   
	float f_c = 1. - f_a - f_b;
	float tau_c = f_b * tau_a / f_a;
	float eT2_a = -TR * (1./T2_a + 1./tau_a);
	float eT2_b = -TR * (1./T2_b + 1./tau_b);
	float eT2_c = -TR * (1./T2_c + 1./tau_c);	
	float eT1_a = -TR * (1./T1_a + 1./tau_a);
	float eT1_b = -TR * (1./T1_b + 1./tau_b);
	float eT1_c = -TR * (1./T1_c + 1./tau_c);
	float k_a   = TR / tau_a;
	float k_b   = TR / tau_b;
	float k_c   = TR / tau_c;
	float phase = rfPhase + (dO * TR * 2. * M_PI);

	float eye[81] = { 1., 0., 0., 0., 0., 0., 0., 0., 0.,
	                   0., 1., 0., 0., 0., 0., 0., 0., 0.,
					   0., 0., 1., 0., 0., 0., 0., 0., 0.,
					   0., 0., 0., 1., 0., 0., 0., 0., 0.,
					   0., 0., 0., 0., 1., 0., 0., 0., 0.,
					   0., 0., 0., 0., 0., 1., 0., 0., 0.,
					   0., 0., 0., 0., 0., 0., 1., 0., 0.,
					   0., 0., 0., 0., 0., 0., 0., 1., 0.,
					   0., 0., 0., 0., 0., 0., 0., 0., 1. };
	float A[81] = {  eT2_a,    k_b,    k_c, phase,    0.,    0.,    0.,    0.,    0.,
					    k_a,  eT2_b,    k_c,    0., phase,    0.,    0.,    0.,    0.,
					    k_a,    k_b,  eT2_c,    0.,    0., phase,    0.,    0.,    0.,
					 -phase,     0.,     0., eT2_a,   k_b,   k_c,    0.,    0.,    0.,
					     0., -phase,     0.,   k_a, eT2_b,   k_c,    0.,    0.,    0.,
					     0.,     0., -phase,   k_a,   k_b, eT2_c,    0.,    0.,    0.,
					     0.,     0.,     0.,    0.,    0.,    0., eT1_a,   k_b,   k_c,
					     0.,     0.,     0.,    0.,    0.,    0.,   k_a, eT1_b,   k_c,
					     0.,     0.,     0.,    0.,    0.,    0.,   k_a,   k_b, eT1_c };
	float R[81] = { 1., 0., 0., 0., 0., 0., 0., 0., 0.,
				     0., 1., 0., 0., 0., 0., 0., 0., 0.,
					 0., 0., 1., 0., 0., 0., 0., 0., 0.,
					 0., 0., 0., 0., 0., 0., 0., 0., 0.,
					 0., 0., 0., 0., 0., 0., 0., 0., 0.,
					 0., 0., 0., 0., 0., 0., 0., 0., 0.,
					 0., 0., 0., 0., 0., 0., 0., 0., 0.,
					 0., 0., 0., 0., 0., 0., 0., 0., 0.,
					 0., 0., 0., 0., 0., 0., 0., 0., 0. };

	float top[81], bottom[81], all[81];
	float M0[9] = { 0., 0., 0., 0., 0., 0., f_a , f_b, f_c }, Mobs[9];
	
	matrixExp(A, 9);
	// Top of matrix divide
	arraySub(top, eye, A, 81);

	for (size_t n = 0; n < nA; n++)
	{
		float ca = cos(B1 * alpha[n]), sa = sin(B1 * alpha[n]);
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

void calcDESPOT1(float *flipAngles, float *spgrVals, int n,
				 float TR, float B1, float *M0, float *T1)
{
	// Linearise the data, then least-squares
	float X[n], Y[n], slope, inter;
	for (int i = 0; i < n; i++)
	{
		X[i] = spgrVals[i] / tan(flipAngles[i] * B1);
		Y[i] = spgrVals[i] / sin(flipAngles[i] * B1);
	}
	linearLeastSquares(X, Y, n, &slope, &inter);	
	*T1 = -TR / log(slope);
	*M0 = inter / (1. - slope);
}

void classicDESPOT2(float *flipAngles, float *ssfpVals, int n,
                    float TR, float T1, float B1, float *p)
{
	// As above, linearise, then least-squares
	// p[0] = M0, p[1] = T2
	float X[n], Y[n], slope, inter;
	for (int i = 0; i < n; i++)
	{
		X[i] = ssfpVals[i] / tan(flipAngles[i] * B1);
		Y[i] = ssfpVals[i] / sin(flipAngles[i] * B1);
	}
	linearLeastSquares(X, Y, n, &slope, &inter);
	float eT1 = exp(-TR / T1);
	p[1] = -TR / log((eT1 - slope) / (1. - slope * eT1));
	float eT2 = exp(-TR / p[1]);
	p[0] = inter * (1. - eT1 * eT2) / (1. - eT1);
}

void simplexDESPOT2(size_t nPhases, size_t *nD, float *phases, float **angles, float **ssfp,
					float TR, float T1, float B1, float *p)
{
	// Gather together all the data
	float *init[4], *c[nPhases], fRes = 0.;
	eval_type *f[nPhases];
	for (int i = 0; i < nPhases; i++)
	{
		c[i] = malloc(4 * sizeof(float));
		c[i][0] = TR; c[i][1] = T1; c[i][2] = B1; c[i][3] = phases[i];
		f[i] = SSFP;
	}
	
	// Use different phases and guesses at B0 to get initial starting location
	for (int i = 0; i < 4; i++)
		init[i] = malloc(3 * sizeof(float));
	classicDESPOT2(angles[0], ssfp[0], nD[0], TR, T1, B1, init[0]);
	classicDESPOT2(angles[0], ssfp[0], nD[0], TR, T1, B1, init[1]);
	classicDESPOT2(angles[1], ssfp[1], nD[1], TR, T1, B1, init[2]);
	classicDESPOT2(angles[1], ssfp[1], nD[1], TR, T1, B1, init[3]);
	init[0][2] = 0.; init[1][2] = 1. / (2. * TR);
	init[2][2] = 0.; init[3][2] = 1. / (2. * TR);
	
	simplex(p, 3, c, nPhases, angles, ssfp, nD, f, init, &fRes);
}

void contractDESPOT2(size_t nPhases, size_t *nD, float *phases, float **flipAngles, float **ssfp,
					 float TR, float T1, float B1, float *p)
{
	float loBounds[2] = {   1.,     0. };
	float hiBounds[2] = { 500., 1./TR };
	float *bounds[2] =  { loBounds, hiBounds };
	bool loConstraints[2] = { true, false };
	bool hiConstraints[2] = { true, false };
	bool *constraints[2] = { loConstraints, hiConstraints };
	float *c[nPhases], fRes = 0.;
	eval_array_type *f[nPhases];
	for (int i = 0; i < nPhases; i++)
	{
		c[i] = malloc(4 * sizeof(float));
		c[i][0] = TR; c[i][1] = T1; c[i][2] = B1; c[i][3] = phases[i];
		f[i] = a1cSSFP;
		arrayScale(ssfp[i], ssfp[i], 1. / arrayMean(ssfp[i], nD[i]), nD[i]);
	}
	regionContraction(p, 2, c, nPhases, flipAngles, ssfp, nD, true, f,
					  bounds, constraints, 4000, 25, 15, 0.05, 0.1, &fRes);
	p[1] = fmod(p[1], 1./TR); // Bring B0 back to one cycle 
	if (p[1] >  .5/TR)        // Move into range -pi/2 to +pi/2
		p[1] -= 1./TR;
	if (p[1] < -.5/TR)
		p[1] += 1./TR;
	p[1] *= 1.e3;             // Finally convert to Hz
}

float calcHIFI(float *flipAngles, float *spgrVals, int nSPGR, float spgrTR,
				float *TI, float *irVals, int nIR, float irFlipAngle, float irTR,
				float *M0, float *T1, float *B1)
{
	// Golden Section Search to find B1	
	// From www.mae.wvu.edu/~smirnov/nr/c10-1.pdf
	float R = 0.61803399; // Golden ratio - 1
	float C = 1 - R;
	float precision = 0.001;	
	
	// Set up initial bracket using some guesses
	float B1_0 = 0.3; float B1_3 = 1.8; float B1_1, B1_2;
	
	// Assemble parameters
	float par[3] = { *M0, *T1, B1_1 };
	float spgrConstants[1] = { spgrTR };
	float irConstants[2] = { irFlipAngle, irTR };
	float spgrRes[nSPGR], irRes[nIR];
	
	par[2] = B1_0;
	calcDESPOT1(flipAngles, spgrVals, nSPGR, spgrTR, par[2], &(par[0]), &(par[1]));
	float res1 = calcResiduals(par, spgrConstants, flipAngles, spgrVals, nSPGR, &SPGR, spgrRes, false) +
	              calcResiduals(par, irConstants, TI, irVals, nIR, &IRSPGR, irRes, false);
	par[2] = B1_3;
	calcDESPOT1(flipAngles, spgrVals, nSPGR, spgrTR, par[2], &(par[0]), &(par[1]));
	float res2 = calcResiduals(par, spgrConstants, flipAngles, spgrVals, nSPGR, SPGR, spgrRes, false) +
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
