/*
 *  DESPOT_Functors.h
 *
 *  Created by Tobias Wood on 16/08/2012.
 *  Copyright (c) 2012-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef DESPOT_Functors_h
#define DESPOT_Functors_h

#include <vector>
#include <array>
#include <map>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;
using namespace Eigen;

typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 9, 9> Matrix9d;
typedef Matrix<double, 9, 1> Vector9d;

struct DESPOTConstants {
	double TR, Trf, phase, B0, B1;
};

// A 3x3 matrix rotation of alpha about X and beta around Z
// Corresponds to RF Flip angle of alpha, and phase-cycling of beta
const Matrix3d RF(const double &alpha, const double &beta)
{
	Matrix3d R;
	R = AngleAxisd(alpha, Vector3d::UnitX()) * AngleAxisd(beta, Vector3d::UnitZ());
	return R;
}

// A 3x3 relaxation matrix for one species, including exchanging and off-resonance
const Matrix3d Relax(const double &R1, const double &R2,
					 const double &k, const double &dw)
{
	Matrix3d r = Matrix3d::Zero();
	
	r(0, 0) = r(1, 1) = -R2 - k;
	r(2, 2) = -R1 - k;
	r(1, 0) = dw;
	r(0, 1) = -dw;
	return r;
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
	Matrix3d O;
	double dw = inHertz * 2. * M_PI;
	O << 0, -dw, 0,
		dw,   0, 0,
		 0,   0, 0;
	return O;
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
// Parameters are { PD, T1, T2 }
//******************************************************************************
VectorXd One_SPGR(const VectorXd &flipAngles,
                  const DESPOTConstants& c, const VectorXd &p)
{
	VectorXd theory(flipAngles.size());
	ArrayXd sa = (flipAngles.array() * c.B1).sin();
	ArrayXd ca = (flipAngles.array() * c.B1).cos();
	double expT1 = exp(-c.TR / p[1]);
	theory = (p[0]*(1. - expT1) * sa) / (1. - expT1*ca);
	return theory;
}

VectorXd One_SSFP(const VectorXd &flipAngles,
                  const DESPOTConstants& c, const VectorXd &p)
{
	Matrix3d A = Relax(p[1], p[2]) + OffResonance(c.B0),
	         R_rf, eATR;
	Vector3d M0 = Vector3d::Zero(), Mobs;
	M0[2] = p[0]; // PD
	
	eATR = (A*c.TR).exp();
	const Vector3d RHS = (Matrix3d::Identity() - eATR) * M0;
	
	VectorXd theory(flipAngles.size());
	for (int i = 0; i < flipAngles.size(); i++) {
		const Matrix3d R_rf = RF(c.B1 * flipAngles[i], c.phase);
		Mobs = (Matrix3d::Identity() - (eATR * R_rf)).partialPivLu().solve(RHS);
		theory[i] = Mobs.head(2).norm();
	}
	return theory;
}

VectorXd One_SSFP_Finite(const VectorXd &flipAngles, const DESPOTConstants& c, const VectorXd &p)
{
	const Matrix3d I = Matrix3d::Identity();
	const Matrix3d R = Relax(p[1], p[2]);
	const Matrix3d O = OffResonance(c.B0);
	Matrix3d C; C = AngleAxisd(c.phase, Vector3d::UnitZ());
	double TE = (c.TR - c.Trf) / 2; // Time AFTER the RF Pulse ends that echo is formed
	
	Matrix3d l1, l2, le;
	l2 = le*le;
	le = ((-R + O) * TE).exp();
	l2 = ((-R + O) * (c.TR - c.Trf)).exp();
	
	Vector3d m0; m0 << 0, 0, p[0];
	Vector3d m2 = (R + O).partialPivLu().solve(R * m0);
	
	VectorXd theory(flipAngles.size());
	for (int i = 0; i < flipAngles.size(); i++) {
		const Matrix3d A = InfinitesimalRF(flipAngles[i] / c.Trf);
		l1 = (-(R + O + A)*(c.Trf)).exp();
		Vector3d m1 = (R + O + A).partialPivLu().solve(R * m0);
		Vector3d mp = C*m2 + (I - l1*C*l2).partialPivLu().solve((I - l1)*(m1 - C*m2));
		Vector3d me = le*(mp - m2) + m2;
		
		theory[i] = me.head(2).norm();
	}
	return theory;
}

//******************************************************************************
#pragma mark Two Component Signals
// Parameters are { PD, T1_a, T2_a, T1_b, T2_b, tau_a, f_a }
//******************************************************************************
VectorXd Two_SPGR(const VectorXd&flipAngles,
                  const DESPOTConstants &c, const VectorXd &p)
{
	Matrix2d A, eATR;
	Vector2d M0, Mobs;
	VectorXd signal(flipAngles.size());
	double k_ab, k_ba, f_a = p[6], f_b = 1. - f_a;
	CalcExchange(p[5], f_a, f_b, k_ab, k_ba);
	M0 << p[0] * f_a, p[0] * f_b;
	A << -((1./p[1]) + k_ab),                    k_ba,
				        k_ab,      -((1./p[3]) + k_ba);
	eATR = (A*c.TR).exp();
	const Vector2d RHS = (Matrix2d::Identity() - eATR) * M0;
	for (int i = 0; i < flipAngles.size(); i++) {
		double a = flipAngles[i];
		Mobs = (Matrix2d::Identity() - eATR*cos(c.B1 * a)).partialPivLu().solve(RHS * sin(c.B1 * a));
		signal[i] = Mobs.sum();
	}
	return signal;
}

VectorXd Two_SSFP(const VectorXd&flipAngles, const DESPOTConstants &c, const VectorXd &p)
{
	VectorXd signal(flipAngles.size());
	Matrix6d A = Matrix6d::Zero(), eATR, eATE, R = Matrix6d::Zero();
	Vector6d M0;
	M0 << 0., 0., p[0] * p[6], 0., 0., p[0] * (1. - p[6]);
	// Set up the 'A' matrix. It's quite complex.
	double TE = c.TR / 2., dw = c.B0 * 2. * M_PI, k_ab, k_ba;
	CalcExchange(p[5], p[6], (1 - p[6]), k_ab, k_ba);
	A.block(0, 0, 3, 3) = A.block(3, 3, 3, 3) = Relax(p[1], p[2]) + OffResonance(dw);
	A += Exchange(k_ab, k_ba);
	A(3, 0) = A(4, 1) = A(5, 2) = k_ab;
	A(0, 3) = A(1, 4) = A(2, 5) = k_ba;
	eATE.noalias() = (A*TE).exp();
	eATR.noalias() = eATE * eATE;
	const Vector6d eyemaM0 = (Matrix6d::Identity() - eATR) * M0;
	for (int i = 0; i < flipAngles.size(); i++) {
		const Matrix3d Rb = RF(c.B1 * flipAngles[i], c.phase);
		R.block(0, 0, 3, 3) = Rb;
		R.block(3, 3, 3, 3) = Rb;
		Vector6d MTR = (Matrix6d::Identity() - eATR * R).partialPivLu().solve(eyemaM0);
		Vector6d Mobs = eATE * R * MTR;
		signal[i] = sqrt(pow(Mobs[0] + Mobs[3], 2.) +
						 pow(Mobs[1] + Mobs[4], 2.));
	}
	return signal;
}

VectorXd Two_SSFP_Finite(const VectorXd &flipAngles, const DESPOTConstants& c, const VectorXd &p)
{
	const Matrix6d I = Matrix6d::Identity();
	Matrix6d R = Matrix6d::Zero();
	R.block(0,0,3,3) = Relax(p[1], p[2]);
	R.block(3,3,3,3) = Relax(p[3], p[4]);
	Matrix6d O = Matrix6d::Zero(); O.block(0,0,3,3) = O.block(3,3,3,3) = OffResonance(c.B0);
	Matrix3d C3; C3 = AngleAxisd(c.phase, Vector3d::UnitZ());
	Matrix6d C = Matrix6d::Zero(); C.block(0,0,3,3) = C.block(3,3,3,3) = C3;
	Matrix6d RpO = R + O;
	double k_ab, k_ba;
	CalcExchange(p[5], p[6], (1 - p[6]), k_ab, k_ba);
	Matrix6d K = Matrix6d::Zero();
	K.block(0,0,3,3).diagonal().setConstant(k_ab);
	K.block(3,3,3,3).diagonal().setConstant(k_ba);
	K.block(3,0,3,3).diagonal().setConstant(-k_ab);
	K.block(0,3,3,3).diagonal().setConstant(-k_ba);
	Matrix6d RpOpK = RpO + K;
	Matrix6d l1, le, temp;
	MatrixExponential<Matrix6d> b2(-(RpOpK) * (c.TR - c.Trf) / 2);
	b2.compute(le);
	const Matrix6d l2 = le*le;
	
	Vector6d m0, mp, me;
	m0 << 0, 0, p[0] * p[6], 0, 0, p[0] * (1-p[6]);
	const Vector6d Rm0 = R * m0;
	const Vector6d m2 = (RpO).partialPivLu().solve(Rm0);
	const Vector6d Cm2 = C * m2;
	
	VectorXd theory(flipAngles.size());
	Matrix6d A = Matrix6d::Zero();
	
	for (int i = 0; i < flipAngles.size(); i++) {
		A.block(0,0,3,3) = A.block(3,3,3,3) = InfinitesimalRF(flipAngles[i] / c.Trf);
		
		MatrixExponential<Matrix6d> b1(-(RpOpK + A)*(c.Trf));
		b1.compute(l1);
		temp.noalias() = -(RpOpK + A)*(c.Trf);
		Vector6d m1 = (RpO + A).partialPivLu().solve(Rm0);
		mp.noalias() = Cm2 + (I - l1*C*l2).partialPivLu().solve((I - l1)*(m1 - Cm2));
		me.noalias() = le*(mp - m2) + m2;
				
		theory[i] = (me.head(3) + me.tail(3)).head(2).norm();
	}
	return theory;
}

//******************************************************************************
#pragma mark Three Component
// Parameters are { PD, T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c }
//******************************************************************************
VectorXd Three_SPGR(const VectorXd&flipAngles, const DESPOTConstants &c, const VectorXd &p)
{
	Matrix3d A, eATR;
	Vector3d M0, Mobs;
	VectorXd signal(flipAngles.size());
	double k_ab, k_ba, f_a = p[8], f_c = p[9], f_b = 1. - f_a - f_c;
	CalcExchange(p[7], f_a, f_b, k_ab, k_ba);
	M0 << p[0] * f_a, p[0] * f_b, p[0] * f_c;
	A << -((1./p[1]) + k_ab),                     k_ba,         0.,
					    k_ab,      -((1./p[3]) + k_ba),         0.,
						  0.,                       0., -(1./p[5]);
	eATR = (A*c.TR).exp();
	const Vector3d RHS = (Matrix3d::Identity() - eATR) * M0;
	for (int i = 0; i < flipAngles.size(); i++) {
		double a = flipAngles[i];
		Mobs = (Matrix3d::Identity() - eATR*cos(c.B1 * a)).partialPivLu().solve(RHS * sin(c.B1 * a));
		signal[i] = Mobs.sum();
	}
	return signal;
}

VectorXd Three_SSFP(const VectorXd&flipAngles, const DESPOTConstants &c, const VectorXd &p)
{
	VectorXd signal(flipAngles.size());
	Matrix6d R = Matrix6d::Zero(), eATE6, eATR6, A6 = Matrix6d::Zero();
	Matrix3d eATE3, eATR3, A3 = Matrix3d::Zero();
	Vector6d M06;
	Vector3d M03;
	double TE = c.TR / 2., dw = c.B0 * 2. * M_PI, k_ab, k_ba,
	       f_a = p[8], f_c = p[9], f_b = 1. - (f_a + f_c);
	M06 << 0., 0., p[0] * f_a, 0., 0., p[0] * f_b;
	M03 << 0., 0., p[0] * f_c;
	CalcExchange(p[7], f_a, f_b, k_ab, k_ba);
	A6.block(0, 0, 3, 3) = A6.block(3, 3, 3, 3) = Relax(p[1], p[2]) + OffResonance(dw);
	A6 += Exchange(k_ab, k_ba);
	A3 = Relax(p[5], p[6]) + OffResonance(dw);
	eATE6 = (A6*TE).exp();
	eATE3 = (A3*TE).exp();
	eATR6.noalias() = eATE6 * eATE6;
	eATR3.noalias() = eATE3 * eATE3;
	const Vector6d eyemaM06 = (Matrix6d::Identity() - eATR6) * M06;
	const Vector3d eyemaM03 = (Matrix3d::Identity() - eATR3) * M03;
	for (int i = 0; i < flipAngles.size(); i++) {
		const Matrix3d Rb = RF(c.B1 * flipAngles[i], c.phase);
		R.block(0, 0, 3, 3) = Rb;
		R.block(3, 3, 3, 3) = Rb;
		Vector6d Mobs6;
		Vector3d Mobs3;
		Mobs6 = eATE6 * R * ((Matrix6d::Identity() - eATR6 * R).partialPivLu().solve(eyemaM06));
		Mobs3 = eATE3 * Rb* ((Matrix3d::Identity() - eATR3 * Rb).partialPivLu().solve(eyemaM03));
		signal[i] = sqrt(pow(Mobs6[0] + Mobs6[3] + Mobs3[0], 2.) +
						 pow(Mobs6[1] + Mobs6[4] + Mobs3[1], 2.));
	}
	return signal;
}

VectorXd Three_SSFP_Finite(const VectorXd &flipAngles, const DESPOTConstants& c, const VectorXd &p)
{
	// Parameters are { PD, T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c }
	const Matrix6d I6 = Matrix6d::Identity();
	const Matrix3d I3 = Matrix3d::Identity();
	Matrix6d R6 = Matrix6d::Zero();
	R6.block(0,0,3,3) = Relax(p[1], p[2]);
	R6.block(3,3,3,3) = Relax(p[3], p[4]);
	Matrix3d R = Relax(p[5], p[6]);
	
	Matrix6d O6 = Matrix6d::Zero();
	Matrix3d O = OffResonance(c.B0);
	O6.block(0,0,3,3) = O6.block(3,3,3,3) = O;
	
	Matrix6d C6 = Matrix6d::Zero();
	Matrix3d C; C = AngleAxisd(c.phase, Vector3d::UnitZ());
	C6.block(0,0,3,3) = C6.block(3,3,3,3) = C;
	
	double k_ab, k_ba, f_a = p[8], f_c = p[9], f_b = 1. - f_a - f_c;
	CalcExchange(p[7], f_a, f_b, k_ab, k_ba);
	Matrix6d K6 = Matrix6d::Zero();
	K6.block(0,0,3,3).diagonal().setConstant(k_ab);
	K6.block(3,3,3,3).diagonal().setConstant(k_ba);
	K6.block(3,0,3,3).diagonal().setConstant(-k_ab);
	K6.block(0,3,3,3).diagonal().setConstant(-k_ba);
	
	Matrix6d l16, l26, le6;
	Matrix3d l1, l2, le;
	MatrixExponential<Matrix6d> b26(-(R6 + O6 + K6) * (c.TR - c.Trf) / 2);
	b26.compute(le6);
	l26 = le6*le6;
	MatrixExponential<Matrix3d> b2(-(R + O) * (c.TR - c.Trf) / 2);
	b2.compute(le);
	l2 = le*le;
	
	Vector9d m0; m0 << 0., 0., p[0]*f_a, 0., 0., p[0]*f_b, 0., 0., p[0]*f_c;
	Vector9d m2;
	m2.head(6) = (R6 + O6).partialPivLu().solve(R6 * m0.head(6));
	m2.tail(3) = (R + O).inverse() * (R * m0.tail(3));
	Matrix6d A6 = Matrix6d::Zero();
	VectorXd theory(flipAngles.size());
	for (int i = 0; i < flipAngles.size(); i++) {
		Matrix3d A = InfinitesimalRF(flipAngles[i] / c.Trf);
		A6.block(0,0,3,3) = A6.block(3,3,3,3) = A;
		MatrixExponential<Matrix6d> b16(-(R6 + O6 + A6 + K6)*(c.Trf));
		MatrixExponential<Matrix3d> b1(-(R + O + A)*(c.Trf));
		b16.compute(l16);
		b1.compute(l1);
		Vector9d m1;
		m1.head(6) = (R6 + O6 + A6).partialPivLu().solve(R6 * m0.head(6));
		m1.tail(3) = (R + O + A).inverse() * (R * m0.tail(3));
		Vector9d mp;
		mp.head(6) = C6*m2.head(6) + (I6 - l16*C6*l26).partialPivLu().solve((I6 - l16)*(m1.head(6) - C6*m2.head(6)));
		mp.tail(3) = C*m2.tail(3) + (I3 - l1*C*l2).partialPivLu().solve((I3 - l1)*(m1.tail(3) - C*m2.tail(3)));
		Vector9d me;
		me.head(6) = le6*(mp.head(6) - m2.head(6)) + m2.head(6);
		me.tail(3) = le*(mp.tail(3) - m2.tail(3)) + m2.tail(3);
		theory[i] = (me.head(3) + me.segment(3,3) + me.tail(3)).head(2).norm();
	}
	return theory;
}

//******************************************************************************
#pragma mark Functor Base Class
//******************************************************************************
// From Nonlinear Tests in Eigen 
template<typename _Scalar, int NX=Dynamic, int NY=Dynamic>
class Functor
{
	public:
		typedef _Scalar Scalar;
		enum {
			InputsAtCompileTime = NX,
			ValuesAtCompileTime = NY
		};
		typedef Matrix<Scalar,InputsAtCompileTime,1> InputType;
		typedef Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
		typedef Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;
		
		const long m_inputs, m_values;
		
		Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
		Functor(long inputs, long values) : m_inputs(inputs), m_values(values) {}
		Functor(long values) : m_inputs(InputsAtCompileTime), m_values(values) {}
		
		virtual ~Functor() {};
		
		long inputs() const { return m_inputs; }
		long values() const { return m_values; }
		
		virtual int operator()(const VectorXd &params, ArrayXd &diffs) const = 0;
		
		virtual const ArrayXd theory(const VectorXd &params) const = 0;
		virtual const ArrayXd signals() const = 0;
};

//******************************************************************************
#pragma mark mcDESPOT Functor
//******************************************************************************
class mcDESPOT : public Functor<double> {
	public:
		enum SignalType {
			SignalSPGR = 0,
			SignalSSFP
		};
		
		enum B0Mode {
			B0_Map = 0,
			B0_Single,
			B0_Multi,
			B0_Bounded,
			B0_MultiBounded
		};
	
		static const int nP(const int &components) {
			switch (components) {
				case 1: return 3;
				case 2: return 7;
				case 3: return 10;
				default:
					cerr << "Cannot create a " << components << "-component mcDESPOT model." << endl;
					exit(EXIT_FAILURE);
			}
		}
		
		static const int nB0(const int &B0Mode, const size_t &nSignals) {
			switch (B0Mode) {
				case B0_Map: return 0; break;
				case B0_Single: return 1; break;
				case B0_Multi: return static_cast<int>(nSignals); break;
				case B0_Bounded: return 1; break;
				case B0_MultiBounded: return static_cast<int>(nSignals); break;
				default:
					cerr << "Invalid B0 mode." << endl;
					exit(EXIT_FAILURE);
			}
		}
		
		static const vector<string> names(const int components) {
			static const map<int, const vector<string> > _namesMap {
				{1, { "PD", "T1", "T2" } },
				{2, { "PD", "T1_a", "T2_a", "T1_b", "T2_b", "tau_a", "f_a"  } },
				{3, { "PD", "T1_a", "T2_a", "T1_b", "T2_b", "T1_c", "T2_c", "tau_a", "f_a", "f_c" } } };
			auto it = _namesMap.find(components);
			if (it == _namesMap.end()) {
				cerr << "Don't have file names for a " << components << "-component mcDESPOT model." << std::endl;
				exit(EXIT_FAILURE);
			} else
			return it->second;
		}
		
		static const ArrayXd &defaultLo(const int components, const int tesla) {
			static ArrayXd c1t3(3), c1t7(3), c2t3(7), c2t7(7), c3t3(10), c3t7(10);
			c1t3 << 0., 0.25, 0.01;
			c1t7 << 0., 0.25, 0.01;
			c2t3 << 0., 0.25, 0.01, 0.75, 0.01, 0.01, 0.001;
			c2t7 << 0., 0.25, 0.01, 0.75, 0.01, 0.01, 0.001;
			c3t3 << 0., 0.35, 0.002, 0.700, 0.075, 3.5, 0.175, 0.05, 0.001, 0.001;
			c3t7 << 0., 0.25, 0.01, 0.75, 0.02, 4.00, 0.15, 0., 0., 0.;
			
			switch (tesla) {
				case 3:
					switch (components) {
						case 1: return c1t3;
						case 2: return c2t3;
						case 3: return c3t3;
					}
				case 7:
					switch (components) {
						case 1: return c1t7;
						case 2: return c2t7;
						case 3: return c3t7;
					}
			}
			std::cerr << "Don't have defaults for a " << components << "-component model." << std::endl;
			exit(EXIT_FAILURE);
		}
		
		static const ArrayXd &defaultHi(const int components, const int tesla) {
			static ArrayXd c1t3(3), c1t7(3), c2t3(7), c2t7(7), c3t3(10), c3t7(10);
			c1t3 << 1.e7, 3.0, 0.25;
			c1t7 << 1.e7, 5.0, 0.10;
			c2t3 << 1.e7, 1.0, 0.05, 1.5, 0.05,  0.5, 0.95;
			c2t7 << 1.e7, 1.0, 0.02, 2.0, 0.05,  0.5, 0.95;
			c3t3 << 1.e7, 0.55, 0.016, 2.0, 0.145,  7.5, 0.5, 0.3, 0.3, 0.7;
			c3t7 << 1.e7, 1.0, 0.02, 2.0, 0.05, 20.0, 0.6, 0.5, 0.95, 0.95;
			
			switch (tesla) {
				case 3:
					switch (components) {
						case 1: return c1t3;
						case 2: return c2t3;
						case 3: return c3t3;
					}
				case 7:
					switch (components) {
						case 1: return c1t7;
						case 2: return c2t7;
						case 3: return c3t7;
					}
			}
			cerr << "Don't have defaults for a " << components << "-component model." << std::endl;
			exit(EXIT_FAILURE);
		}
	
	protected:
		const int _components, _B0Mode;
		long _nP, _nV, _nB0;
		const vector<SignalType> &_types;
		const vector<VectorXd> &_angles, &_signals;
		vector<DESPOTConstants> &_consts;
		const bool _normalise, _debug;
	
	public:
		mcDESPOT(const int components, const vector<SignalType> &types,
				 const vector<VectorXd> &angles, const vector<VectorXd> &signals,
				 vector<DESPOTConstants> &constants,
				 const int &B0Mode, const bool &normalise = false, const bool &debug = false) :
			_components(components), _types(types),
			_angles(angles), _signals(signals), _consts(constants),
			_normalise(normalise), _B0Mode(B0Mode), _debug(debug)
		{
			_nP = nP(components);
			_nB0 = nB0(B0Mode, signals.size());
			_nV = 0;
			for (int i = 0; i < angles.size(); i++) {
				if (angles[i].size() != signals[i].size()) {
					cerr << "Angles and signals size mis-match for signal " << i << endl;
					cerr << "Angles = " << angles[i].size() << " signal = " << signals[i].size() << endl;
					exit(EXIT_FAILURE);
				}
				_nV += angles[i].size();
			}
		}
		
		const bool constraint(const VectorXd &params) {
			// Negative PD or T1/T2 makes no sense
			if ((params[0] < 0.) || (params[1] <= 0.) || (params[2] <= 0.))
				return false;
			
			if (_components == 1) {
				return true;
			} else if (_components == 2) {
				// Check that T1_a, T2_a < T1_b, T2_b and that f_a makes sense
				if ((params[1] < params[3]) &&
					(params[2] < params[4]) &&
					(params[6] <= 1.0))
					return true;
				else
					return false;
			} else if (_components == 3) {
				// Check that T1/2_a < T1/2_b < T1/2_c and that f_a + f_c makes sense 
				if ((params[1] < params[3]) &&
					(params[2] < params[4]) &&
					(params[3] < params[5]) &&
					(params[4] < params[6]) &&
					((params[8] + params[9]) <= 1.0))
					return true;
				else
					return false;
			} else
				return true;
		}
				
		const long inputs() const { return _nP + _nB0; }
		const long values() const { return _nV; }
		
		const ArrayXd signals() const {
			ArrayXd v(values());
			int index = 0;
			for (int i = 0; i < _signals.size(); i++) {
				v.segment(index, _signals[i].size()) = _signals[i];
				//cout << "s:  " << _signals[i].transpose() << endl;
				index += _signals[i].size();
			}
			return v;
		}
		
		const ArrayXd theory(const VectorXd &params) const {
			ArrayXd t(values());
			int index = 0;
			for (int i = 0; i < _signals.size(); i++) {
				ArrayXd theory(_signals[i].size());
				if ((_B0Mode == B0_Single) || (_B0Mode == B0_Bounded))
					_consts[i].B0 = params[_nP];
				else if ((_B0Mode == B0_Multi) || (_B0Mode == B0_MultiBounded))
					_consts[i].B0 = params[_nP + i];
				if (_types[i] == SignalSPGR) {
					switch (_components) {
						case 1: theory = One_SPGR(_angles[i], _consts[i], params.head(_nP)); break;
						case 2: theory = Two_SPGR(_angles[i], _consts[i], params.head(_nP)); break;
						case 3: theory = Three_SPGR(_angles[i], _consts[i], params.head(_nP)); break;
					}
				} else if (_types[i] == SignalSSFP) {
					switch (_components) {
						case 1: theory = One_SSFP(_angles[i], _consts[i], params.head(_nP)); break;
						case 2: theory = Two_SSFP(_angles[i], _consts[i], params.head(_nP)); break;
						case 3: theory = Three_SSFP(_angles[i], _consts[i], params.head(_nP)); break;
					}
				}
				if (_debug) cout << "Params:     " << params.transpose() << endl;
				if (_debug) cout << "Theory:     " << theory.transpose() << endl;
				if (_normalise && (theory.square().sum() > 0.)) theory /= theory.mean();
				if (_debug) cout << "Normalised: " << theory.transpose() << endl;
				t.segment(index, _signals[i].size()) = theory;
				index += _signals[i].size();
			}
			return t;
		}
				
		int operator()(const VectorXd &params, ArrayXd &diffs) const {
			eigen_assert(diffs.size() == values());
			//cout << endl << "operator()" << endl;
			//cout << "p: " << params.transpose() << endl;
			ArrayXd t = theory(params);
			ArrayXd s = signals();
			diffs = t - s;
			//cout << "d:  " << diffs.transpose() << endl;
			//cout << "ds: " << diffs.square().transpose() << endl;
			//cout << "sum:" << diffs.square().sum() << endl;
			return 0;
		}
};

class mcFinite : public mcDESPOT {

	public:
		mcFinite(const int components, const vector<SignalType> &types,
				 const vector<VectorXd> &angles, const vector<VectorXd> &signals,
				 vector<DESPOTConstants> &constants,
				 const int &B0Mode, const bool &normalise = false) :
				mcDESPOT(components, types, angles, signals, constants, B0Mode, normalise)
		{
		}
		
		const ArrayXd theory(const VectorXd &params) const {
			ArrayXd t(values());
			int index = 0;
			for (int i = 0; i < _signals.size(); i++) {
				ArrayXd theory(_signals[i].size());
				if ((_B0Mode == B0_Single) || (_B0Mode == B0_Bounded))
					_consts[i].B0 = params[_nP];
				else if ((_B0Mode == B0_Multi) || (_B0Mode == B0_MultiBounded))
					_consts[i].B0 = params[_nP + i];
				if (_types[i] == SignalSPGR) {
					switch (_components) {
						case 1: theory = One_SPGR(_angles[i], _consts[i], params.head(_nP)); break;
						case 2: theory = Two_SPGR(_angles[i], _consts[i], params.head(_nP)); break;
						case 3: theory = Three_SPGR(_angles[i], _consts[i], params.head(_nP)); break;
					}
				} else if (_types[i] == SignalSSFP) {
					switch (_components) {
						case 1: theory = One_SSFP_Finite(_angles[i], _consts[i], params.head(_nP)); break;
						case 2: theory = Two_SSFP_Finite(_angles[i], _consts[i], params.head(_nP)); break;
						case 3:	theory = Three_SSFP_Finite(_angles[i], _consts[i], params.head(_nP)); break;
					}
				}
				if (_normalise && (theory.square().sum() > 0.)) theory /= theory.mean();
				t.segment(index, _signals[i].size()) = theory;
				index += _signals[i].size();
			}
			return t;
		}
};

//******************************************************************************
#pragma mark DESPOT2FM Functor
//******************************************************************************
class DESPOT2FM : public Functor<double> {
	private:
		long _nV;
		const VectorXd &_angles;
		const vector<VectorXd> &_signals;
		vector<DESPOTConstants> &_consts;
		const double _T1;
		const bool &_normalise, &_fitB0;
	
	public:
		static const vector<string> names() {
			static const vector<string> _names { { "FM_PD", "FM_T2", "FM_B0" } };
			return _names;
		}
		
		static const ArrayXd &defaultLo(const int tesla) {
			static ArrayXd t3(2), t7(2);
			t3 << 0., 0.010;
			t7 << 0., 0.005;
			
			switch (tesla) {
				case 3: return t3; break;
				case 7: return t7; break;
			}
			std::cerr << "Don't have defaults for " << tesla << "tesla." << std::endl;
			exit(EXIT_FAILURE);
		}
		
		static const ArrayXd defaultHi(const int tesla) {
			static ArrayXd t3(2), t7(2);
			t3 << 1.e7, 0.50;
			t7 << 1.e7, 0.25;
			
			switch (tesla) {
				case 3: return t3; break;
				case 7: return t7; break;
			}
			std::cerr << "Don't have defaults for " << tesla << "tesla." << std::endl;
			exit(EXIT_FAILURE);
		}
		
		const bool constraint(const VectorXd &params) {
			if (params[0] < 0.)
				return false;
			else
				return true;
		}
		
		const long inputs() { return _fitB0? 3 : 2; }
		const long values() const { return _nV; }
		
		DESPOT2FM(VectorXd &angles, const vector<VectorXd> &signals,
		          vector<DESPOTConstants> &consts,
				  const double T1, const bool &normalise = false, const bool &fitB0 = true) :
			_angles(angles), _signals(signals),
			_consts(consts),
			_T1(T1), _normalise(normalise), _fitB0(fitB0)
		{
			_nV = 0;
			
			for (int i = 0; i < signals.size(); i++) {
				if (angles.size() != signals[i].size()) {
					std::cerr << "Angles and signals size mis-match for signal " << i << std::endl;
					std::cerr << "Angles = " << angles.size() << " signal = " << signals[i].size() << std::endl;
					exit(EXIT_FAILURE);
				}
				_nV += signals[i].size();
			}
		}
		
		const ArrayXd signals() const {
			ArrayXd v(values());
			int index = 0;
			for (int i = 0; i < _signals.size(); i++) {
				v.segment(index, _signals[i].size()) = _signals[i];
				index += _signals[i].size();
			}
			return v;
		}
		
		// Params are { PD, T2, B0 }
		const ArrayXd theory(const VectorXd &params) const {
			VectorXd withT1(3);
			withT1 << params[0], _T1, params[1];
			
			ArrayXd t(values());
			int index = 0;
			for (int i = 0; i < _signals.size(); i++) {
				if (_fitB0)
					_consts[i].B0 = params[2];
				ArrayXd theory = One_SSFP(_angles, _consts[i], withT1);
				if (_normalise && (theory.square().sum() > 0.)) theory /= theory.mean();
				t.segment(index, _signals[i].size()) = theory;
				index += _signals[i].size();
			}
			return t;
		}
				
		int operator()(const VectorXd &params, ArrayXd &diffs) const {
			eigen_assert(diffs.size() == values());
			ArrayXd t = theory(params);
			ArrayXd s = signals();
			diffs = t - s;
			return 0;
		}
};

#endif
