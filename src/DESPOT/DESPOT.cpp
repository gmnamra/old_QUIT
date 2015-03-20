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
	if (value > low) {
		if (value < high) {
			return value;
		} else {
			return high;
		}
	} else {
		return low;
	}
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

inline const Matrix3d Relax(const double &T1, const double &T2) {
	Matrix3d R;
	R << 1./T2,     0,     0,
	         0, 1./T2,     0,
			 0,     0, 1./T1;
	return R;
}

inline const Matrix3d InfinitesimalRF(const double &dalpha) {
	Matrix3d A;
	A << 0,      0, -dalpha,
	     0,      0,       0,
	     dalpha, 0,       0;
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
const void CalcExchange(const double tau_a, const double f_a, double &f_b, double &k_ab, double &k_ba) {
	const double feps = numeric_limits<float>::epsilon(); // Because we read from float files
	f_b = 1.0 - f_a;
	const double tau_b = f_b * tau_a / f_a;
	k_ab = 1./tau_a; k_ba = 1./tau_b;
	if ((fabs(f_a - 1.) <= feps) || (fabs(f_b - 1.) <= feps)) {
		// Only have 1 component, so no exchange
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
	M.row(0) = PD * ((1. - expT1) * sa) / (1. - expT1*ca);
	return M;
}

MagVector One_SSFP(const ArrayXd &flip, cdbl TR, cdbl phase,
				   cdbl PD, cdbl T1, cdbl T2, cdbl f0) {
	double TE = TR / 2;
	const Vector3d m0(0., 0., PD);
	const Matrix3d E = (-Relax(T1, T2)*TR).exp();
	const Matrix3d E_TE = (-Relax(T1, T2)*TE).exp();
	const Matrix3d O(AngleAxisd(f0*2.*M_PI*TR, Vector3d::UnitZ()));
	const Matrix3d O_TE(AngleAxisd(f0*M_PI*TR, Vector3d::UnitZ()));
	const Matrix3d P(AngleAxisd(phase, Vector3d::UnitZ()));
	MagVector m_e(3, flip.size());
	for (int i = 0; i < flip.size(); i++) {
		const Matrix3d A(AngleAxisd(flip[i], Vector3d::UnitY()));
		const Vector3d m_minus = (Matrix3d::Identity() - P*O*E*A).partialPivLu().solve((1 - exp(-TR/T1)) * m0);
		m_e.col(i).noalias() = O_TE*E_TE*A*m_minus + (1 - exp(-TE/T1)) * m0;
	}
	return m_e;
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

MagVector One_SSFP_Ellipse(const ArrayXd &flip, cdbl TR, cdbl PD, cdbl T1, cdbl T2, cdbl f0) {
	double E1 = exp(-TR / T1);
    double E2 = exp(-TR / T2);

    double theta = M_PI * f0 * TR;
    ArrayXd M = PD * sqrt(E2) * (1 - E1)*sin(flip) / (1 - E1*E2*E2-(E1-E2*E2)*cos(flip));

    MagVector result(3, flip.size());
    result.row(0) = M * cos(theta);
    result.row(1) = M * sin(theta);

    return result;
}

MagVector MP_RAGE(cdbl flip, cdbl TR, const int N, const ArrayXd &TI, cdbl TD,
                  cdbl PD, cdbl T1) {
	const double M0 = PD;
	const double T1s = 1. / (1./T1 - log(cos(flip))/TR);
	const double M0s = M0 * (1. - exp(-TR/T1)) / (1 - exp(-TR/T1s));

	const double A1 = M0s*(1 - exp(-(N*TR)/T1s));
	const double A2 = M0*(1 - exp(-TD/T1));
	const ArrayXd A3 = M0*(1 - exp(-TI/T1));
	const double B1 = exp(-(N*TR)/T1s);
	const double B2 = exp(-TD/T1);
	const ArrayXd B3 = -exp(-TI/T1);

	const ArrayXd A = A3 + A2*B3 + A1*B2*B3;
	const ArrayXd B = B1*B2*B3;
	const ArrayXd M1 = A / (1. - B);

	MagVector M(3, TI.size()); M.setZero();
	M.row(0) = M1 * sin(flip);
	return M;
}

//******************************************************************************
#pragma mark Two Component Signals
//******************************************************************************
MagVector Two_SPGR(const ArrayXd &flip, cdbl TR,
                   cdbl PD, cdbl T1_a, cdbl T1_b, cdbl tau_a, cdbl f_a) {
	Matrix2d A, eATR;
	Vector2d M0, Mobs;
	MagVector signal(3, flip.size()); signal.setZero();
	double k_ab, k_ba, f_b;
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
	double k_ab, k_ba, f_b;
	CalcExchange(tau_a, f_a, f_b, k_ab, k_ba);
	Matrix6d K = Exchange(k_ab, k_ba);
	Matrix6d L = (-(R+O+K)*TR).exp();
	Vector6d M0; M0 << 0., 0., PD * f_a, 0., 0., PD * f_b;
	const Vector6d eyemaM0 = (Matrix6d::Identity() - L) * M0;
	Matrix6d A = Matrix6d::Zero();
	for (int i = 0; i < flip.size(); i++) {
		A.block(0, 0, 3, 3) = Matrix3d(AngleAxisd(flip[i], Vector3d::UnitY()) * AngleAxisd(phase, Vector3d::UnitZ()));
		A.block(3, 3, 3, 3).noalias() = A.block(0, 0, 3, 3);
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
	double k_ab, k_ba, f_b;
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

