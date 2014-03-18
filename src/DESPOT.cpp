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

// There's a bug in matrix log that leads to duplicate symbol definitions if
// this is included in the header file
#include <unsupported/Eigen/MatrixFunctions>

double clamp(double value, double low, double high)
{
	if (value < low)
		return low;
	if (value > high)
		return high;
	return value;
}

#ifdef AGILENT
bool ReadPP(const Nifti &nii, Agilent::ProcPar &pp) {
	const list<Nifti::Extension> &exts = nii.extensions();
	for (auto &e : exts) {
		if (e.code() == NIFTI_ECODE_COMMENT) {
			string s(e.data().begin(), e.data().end());
			stringstream ss(s);
			ss >> pp;
			return true;
		}
	}
	// If we got to here there are no procpar extensions, try the old method
	string path = nii.basePath() + ".procpar";
	ifstream pp_file(path);
	if (pp_file) {
		pp_file >> pp;
		return true;
	}
	return false;
}
#endif

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

	ArrayXd irspgr = ((M0 * sin(B1 * flip) * (1. + irEfficiency * eTI + eFull))).abs();
	return irspgr;
}

//******************************************************************************
#pragma mark Magnetisation Evolution Matrices, helper functions etc.
//******************************************************************************
typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 9, 9> Matrix9d;
typedef Matrix<double, 9, 1> Vector9d;
typedef Matrix<double, 3, Dynamic> MagVector;

// Sum a multi-component magnetisation vector
const MagVector SumMC(const MatrixXd &M_in) {
	MagVector M_out = M_in.topRows(3);
	for (MatrixXd::Index i = 3; i < M_in.rows(); i += 3) {
		M_out += M_in.block(i, 0, 3, M_in.cols());
	}
	return M_out;
}

// Turn a full XYZ magnetisation vector into a complex transverse magnetisation
const VectorXcd SigComplex(const MagVector &M_in) {
	VectorXcd cs(M_in.cols());
	cs.real() = M_in.topRows(1).transpose();
	cs.imag() = M_in.block(1, 0, 1, M_in.cols()).transpose();
	return cs;
}

const VectorXd SigMag(const MagVector &M_in) {
	VectorXd s = M_in.topRows(2).colwise().norm();
	return s;
}

// A 3x3 matrix rotation of alpha about X and beta around Z
// Corresponds to RF Flip angle of alpha, and phase-cycling of beta
const Matrix3d RF(const double &alpha, const double &beta)
{
	Matrix3d R;
	R = AngleAxisd(alpha, Vector3d::UnitX()) * AngleAxisd(beta, Vector3d::UnitZ());
	return R;
}

inline const Matrix3d Relax(const double &T1, const double &T2) {
	Matrix3d R;
	R << 1./T2,     0,     0,
	         0, 1./T2,     0,
			 0,     0, 1./T1;
	return R;
}

inline const Matrix3d InfinitesimalRF(const double &dalpha) {
	Matrix3d A;
	A << 0,      0,       0,
	     0,      0, -dalpha,
	     0, dalpha,       0;
	return A;
}

inline const Matrix3d OffResonance(const double &inHertz) {
	// Minus signs are this way round to make phase cycling go the right way
	Matrix3d O;
	double dw = inHertz * 2. * M_PI;
	O <<  0, dw, 0,
	    -dw,  0, 0,
		  0,  0, 0;
	return O;
}

inline const Matrix3d Spoiling() {
	// Destroy the x- and y- magnetization
	Matrix3d S = Matrix3d::Zero();
	S(2, 2) = 1.;
	return S;
}

inline const Matrix6d Exchange(const double &k_ab, const double &k_ba) {
	Matrix6d K = Matrix6d::Zero();
	K.block(0,0,3,3).diagonal().setConstant(k_ab);
	K.block(3,3,3,3).diagonal().setConstant(k_ba);
	K.block(3,0,3,3).diagonal().setConstant(-k_ab);
	K.block(0,3,3,3).diagonal().setConstant(-k_ba);
	return K;
}

// Calculate the exchange rates from the residence time and fractions
const void CalcExchange(const double tau_a, const double f_a, const double f_b, double &k_ab, double &k_ba) {
	double tau_b = f_b * tau_a / f_a;
	k_ab = 1./tau_a; k_ba = 1./tau_b;
	if ((f_a == 0.) || (f_b == 0.)) { // Only have 1 component, so no exchange
		k_ab = 0.;
		k_ba = 0.;
	}
}

//******************************************************************************
#pragma mark One Component Signals
//******************************************************************************
MagVector One_SPGR(const ArrayXd &flip, cdbl TR, cdbl PD, cdbl T1) {
	MagVector M(3, flip.size()); M.setZero();
	ArrayXd sa = flip.sin();
	ArrayXd ca = flip.cos();
	double expT1 = exp(-TR / T1);
	M.row(1) = PD * ((1. - expT1) * sa) / (1. - expT1*ca);
	return M;
}

MagVector One_SSFP(const ArrayXd &flip, cdbl TR, cdbl phase,
				   cdbl PD, cdbl T1, cdbl T2, cdbl f0) {
	Vector3d M0, Mobs;
	M0 << 0., 0., PD;
	Matrix3d L = (-(Relax(T1, T2) + OffResonance(f0))*TR).exp();
	const Vector3d RHS = (Matrix3d::Identity() - L) * M0;
	MagVector theory(3, flip.size());
	Matrix3d R_rf;
	for (int i = 0; i < flip.size(); i++) {
		const Matrix3d R_rf = RF(flip[i], phase);
		theory.col(i) = (Matrix3d::Identity() - (L * R_rf)).partialPivLu().solve(RHS);
	}
	return theory;
}

MagVector One_SSFP_Finite(const ArrayXd &flip, const bool spoil, cdbl TR, cdbl Trf, cdbl inTE, cdbl phase,
				          cdbl PD, cdbl T1, cdbl T2, cdbl f0) {
	const Matrix3d I = Matrix3d::Identity();
	const Matrix3d O = OffResonance(f0);
	Matrix3d P, R = Relax(T1, T2);
	double TE;
	if (spoil) {
		P = Spoiling();
		TE = inTE - Trf;
		assert(TE > 0.);
	} else {
		P = AngleAxisd(phase, Vector3d::UnitZ());
		TE = (TR - Trf) / 2; // Time AFTER the RF Pulse ends that echo is formed
	}
	
	const Matrix3d RpO = R + O;
	const Matrix3d E_e = (-TE * RpO).exp();
	const Matrix3d E = (-(TR - Trf) * RpO).exp();
	Vector3d m_inf; m_inf << 0, 0, PD;
		
	Matrix3d E_r;
	MagVector result(3, flip.size());
	for (int i = 0; i < flip.size(); i++) {
		const Matrix3d A = InfinitesimalRF(flip(i) / Trf);
		E_r.noalias() = (-Trf * (RpO + A)).exp();
		Vector3d m_rinf = (RpO + A).partialPivLu().solve(R * m_inf);
		Vector3d m_r = (I - E_r*P*E).partialPivLu().solve(E_r*P*(I-E)*m_inf + (I-E_r)*m_rinf);
		Vector3d m_e = E_e*(m_r - m_inf) + m_inf;
		result.col(i) = m_e;
	}
	return result;
}

//******************************************************************************
#pragma mark Two Component Signals
//******************************************************************************
MagVector Two_SPGR(const ArrayXd &flip, cdbl TR,
                   cdbl PD, cdbl T1_a, cdbl T1_b, cdbl tau_a, cdbl f_a) {
	Matrix2d A, eATR;
	Vector2d M0, Mobs;
	MagVector signal(3, flip.size()); signal.setZero();
	double k_ab, k_ba, f_b = 1. - f_a;
	CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
	M0 << f_a, f_b;
	A << -((1./T1_a) + k_ab),                    k_ba,
				       k_ab,      -((1./T1_b) + k_ba);
	eATR = (A*TR).exp();
	const Vector2d RHS = (Matrix2d::Identity() - eATR) * M0;
	for (int i = 0; i < flip.size(); i++) {
		double a = flip[i];
		Mobs = (Matrix2d::Identity() - eATR*cos(a)).partialPivLu().solve(RHS * sin(a));
		signal(1, i) = PD * Mobs.sum();
	}
	return signal;
}

MagVector Two_SSFP(const ArrayXd &flip, const double TR, const double phase,
                   cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b,
				   cdbl tau_a, cdbl f_a, cdbl f0) {
	MagVector signal(3, flip.size());
	Matrix6d R = Matrix6d::Zero();
	R.block(0,0,3,3) = Relax(T1_a, T2_a);
	R.block(3,3,3,3) = Relax(T1_b, T2_b);
	Matrix6d O = Matrix6d::Zero(); O.block(0,0,3,3) = O.block(3,3,3,3) = OffResonance(f0);
	double k_ab, k_ba, f_b = 1. - f_a;
	CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
	Matrix6d K = Exchange(k_ab, k_ba);
	Matrix6d L = (-(R+O+K)*TR).exp();
	Vector6d M0; M0 << 0., 0., PD * f_a, 0., 0., PD * f_b;
	const Vector6d eyemaM0 = (Matrix6d::Identity() - L) * M0;
	Matrix6d A = Matrix6d::Zero();
	for (int i = 0; i < flip.size(); i++) {
		const Matrix3d Ab = RF(flip[i], phase);
		A.block(0, 0, 3, 3) = Ab;
		A.block(3, 3, 3, 3) = Ab;
		Vector6d MTR = (Matrix6d::Identity() - L * A).partialPivLu().solve(eyemaM0);
		signal.col(i) = SumMC(MTR);
	}
	return signal;
}

MagVector Two_SSFP_Finite(const ArrayXd &flip, const bool spoil,
                          cdbl TR, cdbl Trf, cdbl inTE, cdbl phase,
						  cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b,
						  cdbl tau_a, cdbl f_a, cdbl f0) {
	const Matrix6d I = Matrix6d::Identity();
	Matrix6d R = Matrix6d::Zero(), C = Matrix6d::Zero();
	Matrix6d O = Matrix6d::Zero(); O.block(0,0,3,3) = O.block(3,3,3,3) = OffResonance(f0);
	Matrix3d C3;
	double TE;
	if (spoil) {
		TE = inTE - Trf;
		assert(TE > 0.);
		C3 = Spoiling();
		R.block(0,0,3,3) = Relax(T1_a, T2_a);
		R.block(3,3,3,3) = Relax(T1_b, T2_b);
	} else {
		TE = (TR - Trf) / 2; // Time AFTER the RF Pulse ends that echo is formed
		C3 = AngleAxisd(phase, Vector3d::UnitZ());
		R.block(0,0,3,3) = Relax(T1_a, T2_a);
		R.block(3,3,3,3) = Relax(T1_b, T2_b);
	}
	C.block(0,0,3,3) = C.block(3,3,3,3) = C3;
	Matrix6d RpO = R + O;
	double k_ab, k_ba, f_b = 1. - f_a;
	CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
	Matrix6d K = Exchange(k_ab, k_ba);
	Matrix6d RpOpK = RpO + K;
	Matrix6d l1;
	const Matrix6d le = (-(RpOpK)*TE).exp();
	const Matrix6d l2 = (-(RpOpK)*(TR-Trf)).exp();
	
	Vector6d m0, mp, me;
	m0 << 0, 0, f_a * PD, 0, 0, PD * f_b;
	const Vector6d Rm0 = R * m0;
	const Vector6d m2 = (RpO).partialPivLu().solve(Rm0);
	const Vector6d Cm2 = C * m2;
	
	MagVector theory(3, flip.size());
	Matrix6d A = Matrix6d::Zero();
	
	for (int i = 0; i < flip.size(); i++) {
		A.block(0,0,3,3) = A.block(3,3,3,3) = InfinitesimalRF(flip(i) / Trf);
		l1 = (-(RpOpK+A)*Trf).exp();
		Vector6d m1 = (RpO + A).partialPivLu().solve(Rm0);
		mp.noalias() = Cm2 + (I - l1*C*l2).partialPivLu().solve((I - l1)*(m1 - Cm2));
		me.noalias() = le*(mp - m2) + m2;				
		theory.col(i) = SumMC(me);
	}
	return theory;
}

//******************************************************************************
#pragma mark Three Component
//******************************************************************************
// Parameters are { T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c, f0 }
void splitParameters(const VectorXd &p, Ref<VectorXd> p_ab, Ref<VectorXd> p_c);
void splitParameters(const VectorXd &p, Ref<VectorXd> p_ab, Ref<VectorXd> p_c) {
	p_ab.segment(0, 4) = p.segment(0, 4);
	p_ab(4) = p(6); //tau_a
	p_ab(5) = p(7) / (1 - p(8)); // Adjust f_a so f_a + f_b = 1 for the 2c calculation
	p_ab(6) = p(9);
	p_c(0) = p(4); p_c(1) = p(5); p_c(2) = p(9);
}

MagVector Three_SPGR(const ArrayXd &flip, cdbl TR, cdbl PD,
                     cdbl T1_a, cdbl T1_b, cdbl T1_c, cdbl tau_a, cdbl f_a, cdbl f_c) {
	double f_ab = 1. - f_c;
	MagVector m_ab = Two_SPGR(flip, TR, PD * f_ab, T1_a, T1_b, tau_a, f_a / f_ab);
	MagVector m_c  = One_SPGR(flip, TR, PD * f_c, T1_c);
	MagVector r = m_ab + m_c;
	return r;
}

MagVector Three_SSFP(const ArrayXd &flip, cdbl TR, cdbl phase, cdbl PD,
                     cdbl T1_a, cdbl T2_a,
					 cdbl T1_b, cdbl T2_b,
					 cdbl T1_c, cdbl T2_c,
					 cdbl tau_a, cdbl f_a, cdbl f_c, cdbl f0) {
	double f_ab = 1. - f_c;
	MagVector m_ab = Two_SSFP(flip, TR, phase, PD * f_ab, T1_a, T2_a, T1_b, T2_b, tau_a, f_a / f_ab, f0);
	MagVector m_c  = One_SSFP(flip, TR, phase, PD * f_c, T1_c, T2_c, f0);
	MagVector r = m_ab + m_c;
	return r;
}

// Parameters are { T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c, f0 }
MagVector Three_SSFP_Finite(const ArrayXd &flip, const bool spoil,
                            cdbl TR, cdbl Trf, cdbl TE, cdbl ph, cdbl PD,
							cdbl T1_a, cdbl T2_a,
							cdbl T1_b, cdbl T2_b,
							cdbl T1_c, cdbl T2_c,
					        cdbl tau_a, cdbl f_a, cdbl f_c, cdbl f0) {
	double f_ab = 1. - f_c;
	MagVector m_ab = Two_SSFP_Finite(flip, spoil, TR, Trf, TE, ph, PD * f_ab, T1_a, T2_a, T1_b, T2_b, tau_a, f_a / f_ab, f0);
	MagVector m_c  = One_SSFP_Finite(flip, spoil, TR, Trf, TE, ph, PD * f_c, T1_c, T2_c, f0);
	MagVector r = m_ab + m_c;
	return r;
}

