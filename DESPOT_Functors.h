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
	double TR, Trf, TE, phase, B0, B1;
	bool spoil;
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

typedef Matrix<double, 3, Dynamic> MagVector;
// Sum a multi-component magnetisation vector
const MagVector SumMC(const MatrixXd &M_in) {
	MagVector M_out = M_in.topRows(3);
	for (size_t i = 3; i < M_in.rows(); i += 3) {
		M_out += M_in.block(i, 0, 3, M_in.cols());
	}
	return M_out;
}

const VectorXd SigMag(const MagVector &M_in) {
	VectorXd s = M_in.topRows(2).colwise().norm();
	return s;
}

//******************************************************************************
#pragma mark One Component Signals
// Parameters are { PD, T1, T2 }
//******************************************************************************
MagVector One_SPGR(const VectorXd &flipAngles,
                  const DESPOTConstants& c, const VectorXd &p)
{
	MagVector M(3, flipAngles.size()); M.setZero();
	ArrayXd sa = (flipAngles.array() * c.B1).sin();
	ArrayXd ca = (flipAngles.array() * c.B1).cos();
	double expT1 = exp(-c.TR / p[1]);
	M.row(1) = (p[0]*(1. - expT1) * sa) / (1. - expT1*ca);
	return M;
}

MagVector One_SSFP(const VectorXd &flipAngles,
                   const DESPOTConstants& c, const VectorXd &p)
{
	Vector3d M0 = Vector3d::Zero(), Mobs;
	M0[2] = p[0]; // PD
	Matrix3d L = (-(Relax(p[1], p[2]) + OffResonance(c.B0))*c.TR).exp();
	const Vector3d RHS = (Matrix3d::Identity() - L) * M0;
	MagVector theory(3, flipAngles.size());
	Matrix3d R_rf;
	for (int i = 0; i < flipAngles.size(); i++) {
		const Matrix3d R_rf = RF(c.B1 * flipAngles[i], c.phase);
		theory.col(i) = (Matrix3d::Identity() - (L * R_rf)).partialPivLu().solve(RHS);
	}
	return theory;
}

MagVector One_SSFP_Finite(const VectorXd &flipAngles,
                          const DESPOTConstants& c, const VectorXd &p)
{
	const Matrix3d I = Matrix3d::Identity();
	const Matrix3d R = Relax(p[1], p[2]);
	const Matrix3d O = OffResonance(c.B0);
	Matrix3d C;
	double TE;
	if (c.spoil) {
		C = Spoiling();
		TE = c.TE;
	} else {
		C = AngleAxisd(c.phase, Vector3d::UnitZ());
		TE = (c.TR - c.Trf) / 2; // Time AFTER the RF Pulse ends that echo is formed
	}
	
	Matrix3d l1, l2, le;
	le = (-(R + O) * TE).exp();
	l2 = (-(R + O) * (c.TR - c.Trf)).exp();
	Vector3d m0; m0 << 0, 0, p[0];
	Vector3d m2 = (R + O).partialPivLu().solve(R * m0);
	MagVector theory(3, flipAngles.size());
	for (int i = 0; i < flipAngles.size(); i++) {
		const Matrix3d A = InfinitesimalRF(c.B1 * flipAngles[i] / c.Trf);
		l1 = (-(R + O + A)*(c.Trf)).exp();
		Vector3d m1 = (R + O + A).partialPivLu().solve(R * m0);
		Vector3d mp = C*m2 + (I - l1*C*l2).partialPivLu().solve((I - l1)*(m1 - C*m2));
		Vector3d me = le*(mp - m2) + m2;
		theory.col(i) = me;
	}
	return theory;
}

//******************************************************************************
#pragma mark Two Component Signals
// Parameters are { PD, T1_a, T2_a, T1_b, T2_b, tau_a, f_a }
//******************************************************************************
MagVector Two_SPGR(const VectorXd &flipAngles,
                   const DESPOTConstants &c, const VectorXd &p)
{
	Matrix2d A, eATR;
	Vector2d M0, Mobs;
	MagVector signal(3, flipAngles.size()); signal.setZero();
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
		signal(1, i) = Mobs.sum();
	}
	return signal;
}

MagVector Two_SSFP(const VectorXd&flipAngles, const DESPOTConstants &c, const VectorXd &p)
{
	MagVector signal(3, flipAngles.size());
	Vector6d M0; M0 << 0., 0., p[0] * p[6], 0., 0., p[0] * (1. - p[6]);
	Matrix6d R = Matrix6d::Zero();
	R.block(0,0,3,3) = Relax(p[1], p[2]);
	R.block(3,3,3,3) = Relax(p[3], p[4]);
	Matrix6d O = Matrix6d::Zero(); O.block(0,0,3,3) = O.block(3,3,3,3) = OffResonance(c.B0);
	double k_ab, k_ba;
	CalcExchange(p[5], p[6], (1 - p[6]), k_ab, k_ba);
	Matrix6d K = Exchange(k_ab, k_ba);
	Matrix6d L = (-(R+O+K)*c.TR).exp();
	const Vector6d eyemaM0 = (Matrix6d::Identity() - L) * M0;
	Matrix6d A = Matrix6d::Zero();
	for (int i = 0; i < flipAngles.size(); i++) {
		const Matrix3d Ab = RF(c.B1 * flipAngles[i], c.phase);
		A.block(0, 0, 3, 3) = Ab;
		A.block(3, 3, 3, 3) = Ab;
		Vector6d MTR = (Matrix6d::Identity() - L * A).partialPivLu().solve(eyemaM0);
		signal.col(i) = SumMC(MTR);
	}
	return signal;
}

MagVector Two_SSFP_Finite(const VectorXd &flipAngles, const DESPOTConstants& c, const VectorXd &p)
{
	const Matrix6d I = Matrix6d::Identity();
	Matrix6d R = Matrix6d::Zero();
	R.block(0,0,3,3) = Relax(p[1], p[2]);
	R.block(3,3,3,3) = Relax(p[3], p[4]);
	Matrix6d O = Matrix6d::Zero(); O.block(0,0,3,3) = O.block(3,3,3,3) = OffResonance(c.B0);
	Matrix3d C3;
	double TE;
	if (c.spoil) {
		C3 = Spoiling();
		TE = c.TE;
	} else {
		C3 = AngleAxisd(c.phase, Vector3d::UnitZ());
		TE = (c.TR-c.Trf) / 2.;
	}
	Matrix6d C = Matrix6d::Zero(); C.block(0,0,3,3) = C.block(3,3,3,3) = C3;
	Matrix6d RpO = R + O;
	double k_ab, k_ba;
	CalcExchange(p[5], p[6], (1 - p[6]), k_ab, k_ba);
	Matrix6d K = Exchange(k_ab, k_ba);
	Matrix6d RpOpK = RpO + K;
	Matrix6d l1, temp;
	const Matrix6d le = (-(RpOpK)*TE).exp();
	const Matrix6d l2 = (-(RpOpK)*(c.TR-c.Trf)).exp();
	
	Vector6d m0, mp, me;
	m0 << 0, 0, p[0] * p[6], 0, 0, p[0] * (1-p[6]);
	const Vector6d Rm0 = R * m0;
	const Vector6d m2 = (RpO).partialPivLu().solve(Rm0);
	const Vector6d Cm2 = C * m2;
	
	MagVector theory(3, flipAngles.size());
	Matrix6d A = Matrix6d::Zero();
	
	for (int i = 0; i < flipAngles.size(); i++) {
		A.block(0,0,3,3) = A.block(3,3,3,3) = InfinitesimalRF(c.B1 * flipAngles[i] / c.Trf);
		l1 = (-(RpOpK+A)*c.Trf).exp();
		temp.noalias() = -(RpOpK + A)*(c.Trf);
		Vector6d m1 = (RpO + A).partialPivLu().solve(Rm0);
		mp.noalias() = Cm2 + (I - l1*C*l2).partialPivLu().solve((I - l1)*(m1 - Cm2));
		me.noalias() = le*(mp - m2) + m2;				
		theory.col(i) = SumMC(me);
	}
	return theory;
}

//******************************************************************************
#pragma mark Three Component
// Parameters are { PD, T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c }
//******************************************************************************
MagVector Three_SPGR(const VectorXd &flipAngles, const DESPOTConstants &c, const VectorXd &p)
{
	// Parameters are { PD, T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c }
	VectorXd p_1c(3), p_2c(7);
	p_1c(0) = p(0) * p(9); p_1c(1) = p(5); p_1c(2) = p(6);
	p_2c(0) = p(0) * (1. - p(9));
	p_2c.segment(1, 4) = p.segment(1, 4);
	p_2c(5) = p(7); p_2c(6) = p(8) / (1. - p(9)); // Adjust f_a so f_a + f_b = 1 for the 2c calculation
	
	MagVector m_ab = Two_SPGR(flipAngles, c, p_2c);
	MagVector m_c  = One_SPGR(flipAngles, c, p_1c);
	MagVector r = m_ab + m_c;
	return r;
}

MagVector Three_SSFP(const VectorXd &flipAngles, const DESPOTConstants &c, const VectorXd &p)
{
	// Parameters are { PD, T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c }
	VectorXd p_1c(3), p_2c(7);
	p_1c(0) = p(0) * p(9); p_1c(1) = p(5); p_1c(2) = p(6);
	p_2c(0) = p(0) * (1. - p(9));
	p_2c.segment(1, 4) = p.segment(1, 4);
	p_2c(5) = p(7); p_2c(6) = p(8) / (1. - p(9)); // Adjust f_a so f_a + f_b = 1 for the 2c calculation
	
	MagVector m_ab = Two_SSFP(flipAngles, c, p_2c);
	MagVector m_c  = One_SSFP(flipAngles, c, p_1c);
	MagVector r = m_ab + m_c;
	return r;
}

MagVector Three_SSFP_Finite(const VectorXd &flipAngles, const DESPOTConstants& c, const VectorXd &p)
{
	// Parameters are { PD, T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c }
	VectorXd p_1c(3), p_2c(7);
	p_1c(0) = p(0) * p(9); p_1c(1) = p(5); p_1c(2) = p(6);
	p_2c(0) = p(0) * (1. - p(9));
	p_2c.segment(1, 4) = p.segment(1, 4);
	p_2c(5) = p(7); p_2c(6) = p(8) / (1. - p(9)); // Adjust f_a so f_a + f_b = 1 for the 2c calculation
	
	MagVector m_ab = Two_SSFP_Finite(flipAngles, c, p_2c);
	MagVector m_c  = One_SSFP_Finite(flipAngles, c, p_1c);
	MagVector r = m_ab + m_c;
	return r;
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
			c3t3 << 1.e7, 0.55, 0.016, 2.0, 0.145,  7.5, 0.5, 0.3, 0.3, 0.95;
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
			if (_debug) cout << endl << "Params: " << params.transpose() << endl;
			for (int i = 0; i < _signals.size(); i++) {
				MagVector M(3, _signals[i].size());
				if ((_B0Mode == B0_Single) || (_B0Mode == B0_Bounded))
					_consts[i].B0 = params[_nP];
				else if ((_B0Mode == B0_Multi) || (_B0Mode == B0_MultiBounded))
					_consts[i].B0 = params[_nP + i];
				if (_types[i] == SignalSPGR) {
					switch (_components) {
						case 1: M = One_SPGR(_angles[i], _consts[i], params.head(_nP)); break;
						case 2: M = Two_SPGR(_angles[i], _consts[i], params.head(_nP)); break;
						case 3: M = Three_SPGR(_angles[i], _consts[i], params.head(_nP)); break;
					}
				} else if (_types[i] == SignalSSFP) {
					switch (_components) {
						case 1: M = One_SSFP(_angles[i], _consts[i], params.head(_nP)); break;
						case 2: M = Two_SSFP(_angles[i], _consts[i], params.head(_nP)); break;
						case 3: M = Three_SSFP(_angles[i], _consts[i], params.head(_nP)); break;
					}
				}
				ArrayXd theory = SigMag(M);
				if (_normalise && (theory.square().sum() > 0.)) theory /= theory.mean();
				t.segment(index, _signals[i].size()) = theory;
				index += _signals[i].size();
			}
			return t;
		}
				
		int operator()(const VectorXd &params, ArrayXd &diffs) const {
			eigen_assert(diffs.size() == values());
			if (_debug) cout << endl << __PRETTY_FUNCTION__ << endl;
			if (_debug) cout << "Params: " << params.transpose() << endl;
			ArrayXd t = theory(params);
			ArrayXd s = signals();
			diffs = t - s;
			if (_debug) cout << "Diffs:  " << diffs.transpose() << endl;
			if (_debug) cout << "Sum:    " << diffs.square().sum() << endl;
			return 0;
		}
};

#pragma mark mcFinite Functor
class mcFinite : public mcDESPOT {

	public:
		mcFinite(const int components, const vector<SignalType> &types,
				 const vector<VectorXd> &angles, const vector<VectorXd> &signals,
				 vector<DESPOTConstants> &constants,
				 const int &B0Mode, const bool &normalise = false, const bool &debug = false) :
				mcDESPOT(components, types, angles, signals, constants, B0Mode, normalise, debug)
		{
		}
		
		const ArrayXd theory(const VectorXd &params) const {
			ArrayXd t(values());
			int index = 0;
			if (_debug) cout << endl << "Params: " << params.transpose() << endl;
			for (int i = 0; i < _signals.size(); i++) {
				MagVector M(3, _signals[i].size());
				ArrayXd temp(_signals[i].size());
				if ((_B0Mode == B0_Single) || (_B0Mode == B0_Bounded))
					_consts[i].B0 = params[_nP];
				else if ((_B0Mode == B0_Multi) || (_B0Mode == B0_MultiBounded))
					_consts[i].B0 = params[_nP + i];
				switch (_components) {
					case 1: M = One_SSFP_Finite(_angles[i], _consts[i], params.head(_nP)); break;
					case 2: M = Two_SSFP_Finite(_angles[i], _consts[i], params.head(_nP)); break;
					case 3: M = Three_SSFP_Finite(_angles[i], _consts[i], params.head(_nP)); break;
				}
				ArrayXd theory = SigMag(M);
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
		static const vector<string> &names() {
			static const vector<string> _names = { { "FM_PD", "FM_T2", "FM_B0" } };
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
