/*
 *  DESPOT1.c
 *  MacRI
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright 2011 Tobias Wood. All rights reserved.
 *
 */

#include "DESPOT.h"
#include "mathsMatrix.h"
#include "maths2d.h"
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

void aSSFP(double *flipAngle, double *p, double *c, double *ssfp, size_t nA)
{
	double M0 = p[0], T2 = p[1];
	double TR = c[0], T1 = c[1], B0 = c[2], B1 = c[3], offset = c[4];
	
	double eT1 = exp(-TR / T1);
	double eT2 = exp(-TR / T2);
	
	double phase = offset + B0 * TR * 2. * M_PI;
	double sinp = sin(phase);
	double cosp = cos(phase);
	
	for (size_t i = 0; i < nA; i++)
	{
		double sina = sin(B1 * flipAngle[i]);
		double cosa = cos(B1 * flipAngle[i]);

	
		double denom = ((1. - eT1 * cosa) * (1. - eT2 * cosp)) - 
					   (eT2 * (eT1 - cosa) * (eT2 - cosp));
	
		double Mx = ((1 - eT1) * eT2 * sina * (cosp - eT2)) / denom;
		double My = ((1.- eT1) * eT2 * sina * sinp) / denom;
		ssfp[i] = M0 * sqrt(Mx*Mx + My*My);
	}
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
extern int LEV_DEBUG;
void a1cSSFP(gsl_vector *alpha, gsl_vector *p, void *constants, gsl_vector *signal)
{
	SSFP_constants *c = (SSFP_constants *)constants;
	
	Matrix3d A, R_rf = Matrix3d::Zero(), eye = Matrix3d::Identity(), temp1, temp2, all;
	Vector3d M0 = Vector3d::Zero(), Mobs;
	
	M0(2) = gsl_vector_get(p, 0);
	double T2 = gsl_vector_get(p, 1),
		   B0 = gsl_vector_get(p, 2),
		   TR = c->TR, T1 = c->T1, B1 = c->B1, rfPhase = c->rfPhase;		
	double phase = rfPhase + (B0 * TR * 2. * M_PI);
	A << -TR / T2,    phase,       0.,
		   -phase, -TR / T2,       0.,
			   0.,       0., -TR / T1;
	R_rf(0, 0) = 1.;
	MatrixExponential<Matrix3d> expA(A);
	expA.compute(A);
	temp2 = eye - A;
	LLT<Matrix3d> llt;
	for (size_t n = 0; n < alpha->size; n++)
	{
		double a = gsl_vector_get(alpha, n);
		double ca = cos(B1 * a), sa = sin(B1 * a);
		R_rf(1, 1) = R_rf(2, 2) = ca;
		R_rf(1, 2) = sa; R_rf(2, 1) = -sa;
		temp1 = eye - (A * R_rf);
		llt.compute(temp1);
		all = llt.solve(temp2);
		Mobs = all * M0;
		gsl_vector_set(signal, n, sqrt(Mobs(0) * Mobs(0) +
									   Mobs(1) * Mobs(1)));
	}
}

void a2cSPGR(gsl_vector *alpha, gsl_vector *p, void *constants, gsl_vector *signal)
{
	SPGR_constants *c = (SPGR_constants *)constants;
	double T1_a = gsl_vector_get(p, 0),
	       T1_b = gsl_vector_get(p, 1),
		   f_a  = gsl_vector_get(p, 4),
		   f_b = 1. - f_a,
		   tau_a = gsl_vector_get(p, 5),
		   tau_b = f_b * tau_a / f_a,
		   TR = c->TR, B1 = c->B1;
	
	gsl_matrix *A = gsl_matrix_alloc(2, 2);
	gsl_matrix *eye = gsl_matrix_alloc(2, 2);
	gsl_matrix *num = gsl_matrix_alloc(2, 2);
	gsl_matrix *den = gsl_matrix_alloc(2, 2);
	gsl_vector *M0   = gsl_vector_alloc(2);
	gsl_vector *Mobs = gsl_vector_alloc(2);
	matrix_eye(eye);
	double A_d[4]  = {-(TR/T1_a + TR/tau_a), TR/tau_b,
	                    TR/tau_a, -(TR/T1_b + TR/tau_b)};
	matrix_set_array(A, A_d);
	gsl_vector_set(M0, 0, c->M0 * f_a);
	gsl_vector_set(M0, 1, c->M0 * f_b); 
	
	matrix_exp(A);
	for (size_t n = 0; n < signal->size; n++)
	{
		double angle = gsl_vector_get(alpha, n);
		matrix_add_scale(num, eye, 1., A, -1.);
		gsl_matrix_memcpy(den, A);
		gsl_matrix_scale(den, cos(B1 * angle));
		matrix_add_scale(den, eye, 1., den, -1.);
	
		matrix_solve(den, num);
		gsl_matrix_scale(num, sin(B1 * angle));
		matrix_mulv(Mobs, num, M0);
		gsl_vector_set(signal, n, gsl_vector_get(Mobs, 0) + gsl_vector_get(Mobs, 1));
	}
	gsl_matrix_free(A);
	gsl_matrix_free(eye);
	gsl_matrix_free(num);
	gsl_matrix_free(den);
	gsl_vector_free(M0);
	gsl_vector_free(Mobs);
}

void a2cSSFP(gsl_vector *alpha, gsl_vector *p, void *constants, gsl_vector *signal)
{
	SSFP_constants *c = (SSFP_constants *)constants;
	double T1_a = gsl_vector_get(p, 0),
	       T1_b = gsl_vector_get(p, 1),
		   T2_a = gsl_vector_get(p, 2),
		   T2_b = gsl_vector_get(p, 3),
		   f_a  = gsl_vector_get(p, 4),
		   f_b  = 1. - f_a,
		   tau_a = gsl_vector_get(p, 5),
		   tau_b = f_b * tau_a / f_a,
		   B0 = gsl_vector_get(p, 6),
		   TR = c->TR, B1 = c->B1, rfPhase = c->rfPhase;
	double eT2_a = -TR * (1./T2_a + 1./tau_a),
	       eT2_b = -TR * (1./T2_b + 1./tau_b),
		   eT1_a = -TR * (1./T1_a + 1./tau_a),
		   eT1_b = -TR * (1./T1_b + 1./tau_b),
		   k_a   = TR / tau_a,
		   k_b   = TR / tau_b,
		   phase = rfPhase + (B0 * TR * 2. * M_PI);
	
	double A_d[36] = {  eT2_a,    k_b, phase,    0.,    0.,    0.,
					      k_a,  eT2_b,    0., phase,    0.,    0.,
					   -phase,     0., eT2_a,   k_b,    0.,    0.,
					       0., -phase,   k_a, eT2_b,    0.,    0.,
					       0.,     0.,    0.,    0., eT1_a,   k_b,
					       0.,     0.,    0.,    0.,   k_a, eT1_b };
	gsl_matrix *eye = gsl_matrix_alloc(6, 6); matrix_eye(eye);
	gsl_matrix *A   = gsl_matrix_alloc(6, 6); matrix_set_array(A, A_d);
	gsl_matrix *R_rf = gsl_matrix_calloc(6, 6);
	gsl_matrix_set(R_rf, 0, 0, 1.);
	gsl_matrix_set(R_rf, 1, 1, 1.);
	gsl_matrix *top = gsl_matrix_alloc(6, 6);
	gsl_matrix *bottom = gsl_matrix_alloc(6, 6);
	
	gsl_vector *M0 = gsl_vector_calloc(6);
	gsl_vector *Mobs = gsl_vector_alloc(6);
	gsl_vector_set(M0, 4, c->M0 * f_a);
	gsl_vector_set(M0, 5, c->M0 * f_b);
	
	matrix_exp(A);
	// Top of matrix divide
	for (size_t n = 0; n < signal->size; n++)
	{
		matrix_add_scale(top, eye, 1., A, -1.);
		double angle = gsl_vector_get(alpha, n);
		double ca = cos(B1 * angle), sa = sin(B1 * angle);
		gsl_matrix_set(R_rf, 2, 2, ca);
		gsl_matrix_set(R_rf, 3, 3, ca);
		gsl_matrix_set(R_rf, 4, 4, ca);
		gsl_matrix_set(R_rf, 5, 5, ca);
		gsl_matrix_set(R_rf, 2, 4, sa);
		gsl_matrix_set(R_rf, 3, 5, sa);
		gsl_matrix_set(R_rf, 4, 2, -sa);
		gsl_matrix_set(R_rf, 5, 3, -sa);
		// Inverse bracket term (i.e. bottom of matrix divide)
		matrix_mult(bottom, A, R_rf);
		matrix_add_scale(bottom, eye, 1., bottom, -1.);
		// Matrix 'divide'
		matrix_solve(bottom, top);
		matrix_mulv(Mobs, top, M0);
		gsl_vector_set(signal, n, sqrt(pow(gsl_vector_get(Mobs, 0) + gsl_vector_get(Mobs, 1), 2.) +
					                   pow(gsl_vector_get(Mobs, 2) + gsl_vector_get(Mobs, 3), 2.)));
	}
	
	gsl_matrix_free(eye);
	gsl_matrix_free(A);
	gsl_matrix_free(R_rf);
	gsl_matrix_free(top);
	gsl_matrix_free(bottom);
	gsl_vector_free(M0);
	gsl_vector_free(Mobs);
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

double calcDESPOT1(double *flipAngles, double *spgrVals, int n,
				   double TR, double B1, double *M0, double *T1)
{
	// Linearise the data, then least-squares
	double X[n], Y[n], slope, inter, res;
	for (int i = 0; i < n; i++)
	{
		X[i] = spgrVals[i] / tan(flipAngles[i] * B1);
		Y[i] = spgrVals[i] / sin(flipAngles[i] * B1);
	}
	linearLeastSquares(X, Y, n, &slope, &inter, &res);	
	*T1 = -TR / log(slope);
	*M0 = inter / (1. - slope);
	return res;
}

double classicDESPOT2(gsl_vector *flipAngles, gsl_vector *ssfpVals,
                      double TR, double T1, double B1, double *M0, double *T2)
{
	// As above, linearise, then least-squares
	// p[0] = M0, p[1] = T2
	int n = flipAngles->size;
	double X[n], Y[n], slope, inter, residual;
	for (int i = 0; i < n; i++)
	{
		X[i] = gsl_vector_get(ssfpVals, i) / tan(gsl_vector_get(flipAngles, i) * B1);
		Y[i] = gsl_vector_get(ssfpVals, i) / sin(gsl_vector_get(flipAngles, i) * B1);
	}
	linearLeastSquares(X, Y, n, &slope, &inter, &residual);
	double eT1 = exp(-TR / T1);
	*T2 = -TR / log((eT1 - slope) / (1. - slope * eT1));
	double eT2 = exp(-TR / *T2);
	*M0 = inter * (1. - eT1 * eT2) / (1. - eT1);
	return residual;
}
