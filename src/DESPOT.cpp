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

	VectorXd irspgr = ((M0 * sin(B1 * flip) * (1. + irEfficiency * eTI + eFull))).abs();
	return irspgr;
}

//******************************************************************************
#pragma Magnetisation Evolution Matrices, helper functions etc.
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
// Parameters are { T1, T2 } or T2*
//******************************************************************************
MagVector One_SPGR(const VectorXd &p, const ArrayXd &flip, const double TR, const double B1) {
	MagVector M(3, flip.size()); M.setZero();
	ArrayXd sa = (flip * B1).sin();
	ArrayXd ca = (flip * B1).cos();
	double expT1 = exp(-TR / p[0]);
	M.row(1) = ((1. - expT1) * sa) / (1. - expT1*ca);
	return M;
}

MagVector One_SPGR_Echo(const VectorXd &p, const ArrayXd &flip, const double TR, const double TE, const double B1) {
	MagVector M(3, flip.size()); M.setZero();
	ArrayXd sa = (flip * B1).sin();
	ArrayXd ca = (flip * B1).cos();
	double expT1 = exp(-TR / p[0]);
	double expT2 = exp(-TE / p[1]);
	M.row(1) = expT2 * ((1. - expT1) * sa) / (1. - expT1*ca);
	return M;
}

MagVector One_SSFP(const VectorXd &p, const ArrayXd &flip, const double TR, const double phase, const double B1, const double f0) {
	Vector3d M0, Mobs;
	M0 << 0., 0., 1.;
	Matrix3d L = (-(Relax(p[0], p[1]) + OffResonance(f0))*TR).exp();
	const Vector3d RHS = (Matrix3d::Identity() - L) * M0;
	MagVector theory(3, flip.size());
	Matrix3d R_rf;
	for (int i = 0; i < flip.size(); i++) {
		const Matrix3d R_rf = RF(B1 * flip[i], phase);
		theory.col(i) = (Matrix3d::Identity() - (L * R_rf)).partialPivLu().solve(RHS);
	}
	return theory;
}

MagVector One_SSFP_Echo(const VectorXd &p, const ArrayXd &flip, const double TR, const double phase, const double B1, const double f0) {
	Vector3d M0, Mobs;
	M0 << 0., 0., 1.;
	Matrix3d L = (-(Relax(p[0], p[1]) + OffResonance(f0))*TR).exp();
	Matrix3d E = (-(Relax(p[0], p[1]) + OffResonance(f0))*TR/2.).exp();
	const Vector3d RHS = (Matrix3d::Identity() - L) * M0;
	MagVector theory(3, flip.size());
	Matrix3d R_rf;
	for (int i = 0; i < flip.size(); i++) {
		const Matrix3d R_rf = RF(B1 * flip[i], phase);
		theory.col(i) = E * (Matrix3d::Identity() - (L * R_rf)).partialPivLu().solve(RHS);
	}
	return theory;
}

// Parameters are { T1, T2, delta_f }
MagVector One_SSFP_Finite(const VectorXd &p, const ArrayXd &flip, const bool spoil,
                          const double TR, const double Trf, const double inTE,
						  const double phase, const double B1, const double f0) {
	const Matrix3d I = Matrix3d::Identity();
	const Matrix3d O = OffResonance(f0);
	Matrix3d C, R;
	double TE;
	if (spoil) {
		C = Spoiling();
		TE = inTE - Trf;
		R = Relax(p[0], 1./(1./p[1]+p[2])); // For SPGR use T2*
	} else {
		C = AngleAxisd(phase, Vector3d::UnitZ());
		TE = (TR - Trf) / 2; // Time AFTER the RF Pulse ends that echo is formed
		R = Relax(p[0], p[1]); // For SSFP just use T2
	}
	
	Matrix3d l1, l2, le;
	le = (-(R + O) * TE).exp();
	l2 = (-(R + O) * (TR - Trf)).exp();
	Vector3d m0; m0 << 0, 0, 1.;
	Vector3d m2 = (R + O).partialPivLu().solve(R * m0);
	MagVector theory(3, flip.size());
	
	for (int i = 0; i < flip.size(); i++) {
		const Matrix3d A = InfinitesimalRF(B1 * flip(i) / Trf);
		l1 = (-(R + O + A)*(Trf)).exp();
		Vector3d m1 = (R + O + A).partialPivLu().solve(R * m0);
		Vector3d mp = C*m2 + (I - l1*C*l2).partialPivLu().solve((I - l1)*(m1 - C*m2));
		Vector3d me = le*(mp - m2) + m2;
		theory.col(i) = me;
	}
	return theory;
}

//******************************************************************************
#pragma mark Two Component Signals
//******************************************************************************
// Parameters are { T1_a, T2_a, T1_b, T2_b, tau_a, f_a }
MagVector Two_SPGR(const VectorXd &p, const ArrayXd &flip, const double TR, const double B1) {
	Matrix2d A, eATR;
	Vector2d M0, Mobs;
	MagVector signal(3, flip.size()); signal.setZero();
	double k_ab, k_ba, f_a = p[5], f_b = 1. - f_a;
	CalcExchange(p[4], f_a, f_b, k_ab, k_ba);
	M0 << f_a, f_b;
	A << -((1./p[0]) + k_ab),                    k_ba,
				        k_ab,      -((1./p[2]) + k_ba);
	eATR = (A*TR).exp();
	const Vector2d RHS = (Matrix2d::Identity() - eATR) * M0;
	for (int i = 0; i < flip.size(); i++) {
		double a = flip[i] * B1;
		Mobs = (Matrix2d::Identity() - eATR*cos(a)).partialPivLu().solve(RHS * sin(B1 * a));
		signal(1, i) = Mobs.sum();
	}
	return signal;
}

MagVector Two_SPGR_Echo(const VectorXd &p, const ArrayXd &flip, const double TR, const double TE, const double B1) {
	Matrix2d A, eATR, eATE;
	Vector2d M0, Mobs;
	MagVector signal(3, flip.size()); signal.setZero();
	double k_ab, k_ba, f_a = p[5], f_b = 1. - f_a;
	CalcExchange(p[4], f_a, f_b, k_ab, k_ba);
	M0 << f_a, f_b;
	A << -((1./p[0]) + k_ab),                    k_ba,
				        k_ab,      -((1./p[2]) + k_ba);
	eATR = (A*TR).exp();
	eATE = (A*TE).exp();
	const Vector2d RHS = (Matrix2d::Identity() - eATR) * M0;
	for (int i = 0; i < flip.size(); i++) {
		double a = flip[i] * B1;
		Mobs = eATE * (Matrix2d::Identity() - eATR*cos(a)).partialPivLu().solve(RHS * sin(B1 * a));
		signal(1, i) = Mobs.sum();
	}
	return signal;
}

MagVector Two_SSFP(const VectorXd &p, const ArrayXd &flip, const double TR, const double phase, const double B1, const double f0) {
	MagVector signal(3, flip.size());
	Vector6d M0; M0 << 0., 0., p[5], 0., 0., (1. - p[5]);
	Matrix6d R = Matrix6d::Zero();
	R.block(0,0,3,3) = Relax(p[0], p[1]);
	R.block(3,3,3,3) = Relax(p[2], p[3]);
	Matrix6d O = Matrix6d::Zero(); O.block(0,0,3,3) = O.block(3,3,3,3) = OffResonance(f0);
	double k_ab, k_ba;
	CalcExchange(p[4], p[5], (1 - p[5]), k_ab, k_ba);
	Matrix6d K = Exchange(k_ab, k_ba);
	Matrix6d L = (-(R+O+K)*TR).exp();
	const Vector6d eyemaM0 = (Matrix6d::Identity() - L) * M0;
	Matrix6d A = Matrix6d::Zero();
	for (int i = 0; i < flip.size(); i++) {
		const Matrix3d Ab = RF(B1 * flip[i], phase);
		A.block(0, 0, 3, 3) = Ab;
		A.block(3, 3, 3, 3) = Ab;
		Vector6d MTR = (Matrix6d::Identity() - L * A).partialPivLu().solve(eyemaM0);
		signal.col(i) = SumMC(MTR);
	}
	return signal;
}

MagVector Two_SSFP_Echo(const VectorXd &p, const ArrayXd &flip, const double TR, const double phase, const double B1, const double f0) {
	MagVector signal(3, flip.size());
	Vector6d M0; M0 << 0., 0., p[5], 0., 0., (1. - p[5]);
	Matrix6d R = Matrix6d::Zero();
	R.block(0,0,3,3) = Relax(p[0], p[1]);
	R.block(3,3,3,3) = Relax(p[2], p[3]);
	Matrix6d O = Matrix6d::Zero(); O.block(0,0,3,3) = O.block(3,3,3,3) = OffResonance(f0);
	double k_ab, k_ba;
	CalcExchange(p[4], p[5], (1 - p[5]), k_ab, k_ba);
	Matrix6d K = Exchange(k_ab, k_ba);
	Matrix6d L = (-(R+O+K)*TR).exp();
	Matrix6d E = (-(R+O+K)*TR/2.).exp();
	const Vector6d eyemaM0 = (Matrix6d::Identity() - L) * M0;
	Matrix6d A = Matrix6d::Zero();
	for (int i = 0; i < flip.size(); i++) {
		const Matrix3d Ab = RF(B1 * flip[i], phase);
		A.block(0, 0, 3, 3) = Ab;
		A.block(3, 3, 3, 3) = Ab;
		Vector6d MTE = E * (Matrix6d::Identity() - L * A).partialPivLu().solve(eyemaM0);
		signal.col(i) = SumMC(MTE);
	}
	return signal;
}

// Parameters are { T1_a, T2_a, T1_b, T2_b, tau_a, f_a, delta_f }
MagVector Two_SSFP_Finite(const VectorXd &p, const ArrayXd &flip, const bool spoil,
                          const double TR, const double Trf, const double inTE,
						  const double phase, const double B1, const double f0) {
	const Matrix6d I = Matrix6d::Identity();
	Matrix6d R = Matrix6d::Zero(), C = Matrix6d::Zero();
	Matrix6d O = Matrix6d::Zero(); O.block(0,0,3,3) = O.block(3,3,3,3) = OffResonance(f0);
	Matrix3d C3;
	double TE;
	if (spoil) {
		TE = inTE - Trf;
		C3 = Spoiling();
		R.block(0,0,3,3) = Relax(p[0], 1./(1./p[1] + p[6])); // For SPGR use T2*
		R.block(3,3,3,3) = Relax(p[2], 1./(1./p[3] + p[6]));
	} else {
		TE = (TR - Trf) / 2; // Time AFTER the RF Pulse ends that echo is formed
		C3 = AngleAxisd(phase, Vector3d::UnitZ());
		R.block(0,0,3,3) = Relax(p[0], p[1]); // For SSFP just use T2
		R.block(3,3,3,3) = Relax(p[2], p[3]);
	}
	C.block(0,0,3,3) = C.block(3,3,3,3) = C3;
	Matrix6d RpO = R + O;
	double k_ab, k_ba;
	CalcExchange(p[4], p[5], (1. - p[5]), k_ab, k_ba);
	Matrix6d K = Exchange(k_ab, k_ba);
	Matrix6d RpOpK = RpO + K;
	Matrix6d l1, temp;
	const Matrix6d le = (-(RpOpK)*TE).exp();
	const Matrix6d l2 = (-(RpOpK)*(TR-Trf)).exp();
	
	Vector6d m0, mp, me;
	m0 << 0, 0, p[5], 0, 0, (1. - p[5]);
	const Vector6d Rm0 = R * m0;
	const Vector6d m2 = (RpO).partialPivLu().solve(Rm0);
	const Vector6d Cm2 = C * m2;
	
	MagVector theory(3, flip.size());
	Matrix6d A = Matrix6d::Zero();
	
	for (int i = 0; i < flip.size(); i++) {
		A.block(0,0,3,3) = A.block(3,3,3,3) = InfinitesimalRF(B1 * flip(i) / Trf);
		l1 = (-(RpOpK+A)*Trf).exp();
		temp.noalias() = -(RpOpK + A)*(Trf);
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
// Parameters are { T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c }
void splitParameters(const VectorXd &p, Ref<VectorXd> p_ab, Ref<VectorXd> p_c);
void splitParameters(const VectorXd &p, Ref<VectorXd> p_ab, Ref<VectorXd> p_c) {
	p_ab.segment(0, 4) = p.segment(0, 4);
	p_ab(4) = p(6); //tau_a
	p_ab(5) = p(7) / (1 - p(8)); // Adjust f_a so f_a + f_b = 1 for the 2c calculation
	p_c(0) = p(4); p_c(1) = p(5);
}

MagVector Three_SPGR(const VectorXd &p, const ArrayXd &flip, const double TR, const double B1) {
	VectorXd p_ab(6), p_c(2);
	splitParameters(p, p_ab, p_c);
	MagVector m_ab = Two_SPGR(p_ab, flip, TR, B1);
	MagVector m_c  = One_SPGR(p_c, flip, TR, B1);
	MagVector r = (m_ab * (1. - p(8))) + (m_c * p(8));
	return r;
}

MagVector Three_SPGR_Echo(const VectorXd &p, const ArrayXd &flip, const double TR, const double TE, const double B1) {
	VectorXd p_ab(6), p_c(2);
	splitParameters(p, p_ab, p_c);
	MagVector m_ab = Two_SPGR_Echo(p_ab, flip, TR, TE, B1);
	MagVector m_c  = One_SPGR_Echo(p_c, flip, TR, TE, B1);
	MagVector r = (m_ab * (1. - p(8))) + (m_c * p(8));
	return r;
}

MagVector Three_SSFP(const VectorXd &p, const ArrayXd &flip, const double TR, const double phase, const double B1, const double f0) {
	VectorXd p_ab(6), p_c(2);
	splitParameters(p, p_ab, p_c);
	MagVector m_ab = Two_SSFP(p_ab, flip, TR, phase, B1, f0);
	MagVector m_c  = One_SSFP(p_c, flip, TR, phase, B1, f0);
	MagVector r = (m_ab * (1. - p(8))) + (m_c * p(8));
	return r;
}

MagVector Three_SSFP_Echo(const VectorXd &p, const ArrayXd &flip, const double TR, const double phase, const double B1, const double f0) {
	VectorXd p_ab(6), p_c(2);
	splitParameters(p, p_ab, p_c);
	MagVector m_ab = Two_SSFP_Echo(p_ab, flip, TR, phase, B1, f0);
	MagVector m_c  = One_SSFP_Echo(p_c, flip, TR, phase, B1, f0);
	MagVector r = (m_ab * (1. - p(8))) + (m_c * p(8));
	return r;
}

// Parameters are { T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c, delta_f }
MagVector Three_SSFP_Finite(const VectorXd &p, const ArrayXd &flip, const bool spoil,
                            const double TR, const double Trf, const double TE,
						    const double ph, const double B1, const double f0) {
	VectorXd p_ab(7), p_c(3);
	p_ab.segment(0, 4) = p.segment(0, 4);
	p_ab(4) = p(6); //tau_a
	p_ab(5) = p(7) / (1 - p(8)); // Adjust f_a so f_a + f_b = 1 for the 2c calculation
	p_ab(6) = p(9);
	p_c(0) = p(4); p_c(1) = p(5);
	p_c(2) = p(9);
	MagVector m_ab = Two_SSFP_Finite(p_ab, flip, spoil, TR, Trf, TE, ph, B1, f0);
	MagVector m_c  = One_SSFP_Finite(p_c,  flip, spoil, TR, Trf, TE, ph, B1, f0);
	MagVector r = (m_ab * (1. - p(8))) + (m_c * p(8));
	return r;
}

