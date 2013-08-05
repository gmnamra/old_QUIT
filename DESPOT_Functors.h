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
#include <exception>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;
using namespace Eigen;

typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 9, 9> Matrix9d;
typedef Matrix<double, 9, 1> Vector9d;

class DESPOTData {
	private:
		VectorXd _flip, _signal;
	public:
		// B0 is field-strength in T, f0_off is off-resonance in Hz
		double TR, Trf, TE, phase, f0_off, B1;
		bool spoil;
	
		DESPOTData() : TR(0.), Trf(0.), TE(0.), phase(0.), f0_off(0.), B1(0.), spoil(false), _flip(), _signal() {};
		DESPOTData(const size_t nData, bool inSpoil, double inTR, double inTrf, double inTE = 0., double inPhase = M_PI, double inf0_off = 0., double inB1 = 1.) :
			TR(inTR), Trf(inTrf), TE(inTE), phase(inPhase), f0_off(inf0_off), B1(inB1), spoil(inSpoil), _flip(), _signal()
		{
			_flip.resize(nData);
			_signal.resize(nData);
		};
		DESPOTData(const size_t nData, double inSpoil, double inTR, double inTrf, double inTE, double inPhase, double inf0_off, double inB1) = delete;
		
		const size_t size() const { return _flip.rows(); };
		void resize(const size_t n) {
			_flip.resize(n);
			_signal.resize(n);
		};
		
		const VectorXd &flip() const { return _flip; };
		void setFlip(const VectorXd &inFlip) {
			assert(inFlip.rows() == _flip.rows());
			_flip = inFlip;
		};
		const VectorXd &signal() const { return _signal; };
		void setSignal(const VectorXd &inSig) {
			assert(inSig.rows() == _signal.rows());
			_signal = inSig;
		};
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
MagVector One_SPGR(const DESPOTData &d, const VectorXd &p)
{
	MagVector M(3, d.flip().size()); M.setZero();
	ArrayXd sa = (d.flip().array() * d.B1).sin();
	ArrayXd ca = (d.flip().array() * d.B1).cos();
	double expT1 = exp(-d.TR / p[1]);
	M.row(1) = (p[0]*(1. - expT1) * sa) / (1. - expT1*ca);
	return M;
}

MagVector One_SSFP(const DESPOTData &d, const VectorXd &p)
{
	Vector3d M0 = Vector3d::Zero(), Mobs;
	M0[2] = p[0]; // PD
	Matrix3d L = (-(Relax(p[1], p[2]) + OffResonance(d.f0_off))*d.TR).exp();
	const Vector3d RHS = (Matrix3d::Identity() - L) * M0;
	MagVector theory(3, d.flip().size());
	Matrix3d R_rf;
	for (int i = 0; i < d.flip().size(); i++) {
		const Matrix3d R_rf = RF(d.B1 * d.flip()[i], d.phase);
		theory.col(i) = (Matrix3d::Identity() - (L * R_rf)).partialPivLu().solve(RHS);
	}
	return theory;
}

MagVector One_SSFP_Finite(const DESPOTData &d, const VectorXd &p)
{
	const Matrix3d I = Matrix3d::Identity();
	const Matrix3d O = OffResonance(d.f0_off);
	Matrix3d C, R;
	double TE;
	if (d.spoil) {
		C = Spoiling();
		TE = d.TE - d.Trf;
		R = Relax(p[1], 1./(1./p[2]+d.f0_off/10.)); // For SPGR use T2*
	} else {
		C = AngleAxisd(d.phase, Vector3d::UnitZ());
		TE = (d.TR - d.Trf) / 2; // Time AFTER the RF Pulse ends that echo is formed
		R = Relax(p[1], p[2]); // For SSFP just use T2
	}
	
	Matrix3d l1, l2, le;
	le = (-(R + O) * TE).exp();
	l2 = (-(R + O) * (d.TR - d.Trf)).exp();
	Vector3d m0; m0 << 0, 0, p[0];
	Vector3d m2 = (R + O).partialPivLu().solve(R * m0);
	MagVector theory(3, d.flip().size());
	for (int i = 0; i < d.flip().size(); i++) {
		const Matrix3d A = InfinitesimalRF(d.B1 * d.flip()[i] / d.Trf);
		l1 = (-(R + O + A)*(d.Trf)).exp();
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
MagVector Two_SPGR(const DESPOTData &d, const VectorXd &p)
{
	Matrix2d A, eATR;
	Vector2d M0, Mobs;
	MagVector signal(3, d.flip().size()); signal.setZero();
	double k_ab, k_ba, f_a = p[6], f_b = 1. - f_a;
	CalcExchange(p[5], f_a, f_b, k_ab, k_ba);
	M0 << p[0] * f_a, p[0] * f_b;
	A << -((1./p[1]) + k_ab),                    k_ba,
				        k_ab,      -((1./p[3]) + k_ba);
	eATR = (A*d.TR).exp();
	const Vector2d RHS = (Matrix2d::Identity() - eATR) * M0;
	for (int i = 0; i < d.flip().size(); i++) {
		double a = d.flip()[i];
		Mobs = (Matrix2d::Identity() - eATR*cos(d.B1 * a)).partialPivLu().solve(RHS * sin(d.B1 * a));
		signal(1, i) = Mobs.sum();
	}
	return signal;
}

MagVector Two_SSFP(const DESPOTData &d, const VectorXd &p)
{
	MagVector signal(3, d.flip().size());
	Vector6d M0; M0 << 0., 0., p[0] * p[6], 0., 0., p[0] * (1. - p[6]);
	Matrix6d R = Matrix6d::Zero();
	R.block(0,0,3,3) = Relax(p[1], p[2]);
	R.block(3,3,3,3) = Relax(p[3], p[4]);
	Matrix6d O = Matrix6d::Zero(); O.block(0,0,3,3) = O.block(3,3,3,3) = OffResonance(d.f0_off);
	double k_ab, k_ba;
	CalcExchange(p[5], p[6], (1 - p[6]), k_ab, k_ba);
	Matrix6d K = Exchange(k_ab, k_ba);
	Matrix6d L = (-(R+O+K)*d.TR).exp();
	const Vector6d eyemaM0 = (Matrix6d::Identity() - L) * M0;
	Matrix6d A = Matrix6d::Zero();
	for (int i = 0; i < d.flip().size(); i++) {
		const Matrix3d Ab = RF(d.B1 * d.flip()[i], d.phase);
		A.block(0, 0, 3, 3) = Ab;
		A.block(3, 3, 3, 3) = Ab;
		Vector6d MTR = (Matrix6d::Identity() - L * A).partialPivLu().solve(eyemaM0);
		signal.col(i) = SumMC(MTR);
	}
	return signal;
}

MagVector Two_SSFP_Finite(const DESPOTData &d, const VectorXd &p)
{
	const Matrix6d I = Matrix6d::Identity();
	Matrix6d R = Matrix6d::Zero(), C = Matrix6d::Zero();
	Matrix6d O = Matrix6d::Zero(); O.block(0,0,3,3) = O.block(3,3,3,3) = OffResonance(d.f0_off);
	Matrix3d C3;
	double TE;
	if (d.spoil) {
		TE = d.TE - d.Trf;
		C3 = Spoiling();
		R.block(0,0,3,3) = Relax(p[1], 1./(1./p[2] + d.f0_off/10.)); // For SPGR use T2*
		R.block(3,3,3,3) = Relax(p[3], 1./(1./p[4] + d.f0_off/10.));
	} else {
		TE = (d.TR - d.Trf) / 2; // Time AFTER the RF Pulse ends that echo is formed
		C3 = AngleAxisd(d.phase, Vector3d::UnitZ());
		R.block(0,0,3,3) = Relax(p[1], p[2]); // For SSFP just use T2
		R.block(3,3,3,3) = Relax(p[3], p[4]);
	}
	C.block(0,0,3,3) = C.block(3,3,3,3) = C3;
	Matrix6d RpO = R + O;
	double k_ab, k_ba;
	CalcExchange(p[5], p[6], (1 - p[6]), k_ab, k_ba);
	Matrix6d K = Exchange(k_ab, k_ba);
	Matrix6d RpOpK = RpO + K;
	Matrix6d l1, temp;
	const Matrix6d le = (-(RpOpK)*TE).exp();
	const Matrix6d l2 = (-(RpOpK)*(d.TR-d.Trf)).exp();
	
	Vector6d m0, mp, me;
	m0 << 0, 0, p[0] * p[6], 0, 0, p[0] * (1-p[6]);
	const Vector6d Rm0 = R * m0;
	const Vector6d m2 = (RpO).partialPivLu().solve(Rm0);
	const Vector6d Cm2 = C * m2;
	
	MagVector theory(3, d.flip().size());
	Matrix6d A = Matrix6d::Zero();
	
	for (int i = 0; i < d.flip().size(); i++) {
		A.block(0,0,3,3) = A.block(3,3,3,3) = InfinitesimalRF(d.B1 * d.flip()[i] / d.Trf);
		l1 = (-(RpOpK+A)*d.Trf).exp();
		temp.noalias() = -(RpOpK + A)*(d.Trf);
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
MagVector Three_SPGR(const DESPOTData &d, const VectorXd &p)
{
	// Parameters are { PD, T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c }
	VectorXd p_1c(3), p_2c(7);
	p_1c(0) = p(0) * p(9); p_1c(1) = p(5); p_1c(2) = p(6);
	p_2c(0) = p(0) * (1. - p(9));
	p_2c.segment(1, 4) = p.segment(1, 4);
	p_2c(5) = p(7); p_2c(6) = p(8) / (1. - p(9)); // Adjust f_a so f_a + f_b = 1 for the 2c calculation
	
	MagVector m_ab = Two_SPGR(d, p_2c);
	MagVector m_c  = One_SPGR(d, p_1c);
	MagVector r = m_ab + m_c;
	return r;
}

MagVector Three_SSFP(const DESPOTData &d, const VectorXd &p)
{
	// Parameters are { PD, T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c }
	VectorXd p_1c(3), p_2c(7);
	p_1c(0) = p(0) * p(9); p_1c(1) = p(5); p_1c(2) = p(6);
	p_2c(0) = p(0) * (1. - p(9));
	p_2c.segment(1, 4) = p.segment(1, 4);
	p_2c(5) = p(7); p_2c(6) = p(8) / (1. - p(9)); // Adjust f_a so f_a + f_b = 1 for the 2c calculation
	
	MagVector m_ab = Two_SSFP(d, p_2c);
	MagVector m_c  = One_SSFP(d, p_1c);
	MagVector r = m_ab + m_c;
	return r;
}

MagVector Three_SSFP_Finite(const DESPOTData &d, const VectorXd &p)
{
	// Parameters are { PD, T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c }
	VectorXd p_1c(3), p_2c(7);
	p_1c(0) = p(0) * p(9); p_1c(1) = p(5); p_1c(2) = p(6);
	p_2c(0) = p(0) * (1. - p(9));
	p_2c.segment(1, 4) = p.segment(1, 4);
	p_2c(5) = p(7); p_2c(6) = p(8) / (1. - p(9)); // Adjust f_a so f_a + f_b = 1 for the 2c calculation
	
	MagVector m_ab = Two_SSFP_Finite(d, p_2c);
	MagVector m_c  = One_SSFP_Finite(d, p_1c);
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
		enum class Components {
			One, Two, Three
		};
		static string to_string(const Components& c) {
			switch (c) {
				case Components::One: return "one";
				case Components::Two: return "two";
				case Components::Three: return "three";
			}
		};
		enum class FieldStrength {
			Three, Seven, Unknown
		};
		static string to_string(const FieldStrength& f) {
			switch (f) {
				case FieldStrength::Three: return "3";
				case FieldStrength::Seven: return "7";
				case FieldStrength::Unknown: return "User";
			}
		}
		enum class OffResMode {
			Map = 0, Single, Multi, Bounded, MultiBounded
		};
		
		static const int nP(const Components &c) {
			switch (c) {
				case Components::One: return 3;
				case Components::Two: return 7;
				case Components::Three: return 10;
			}
		}
		
		static const int nOffRes(const OffResMode &m, const size_t &nSignals) {
			switch (m) {
				case OffResMode::Map: return 0;
				case OffResMode::Single: return 1;
				case OffResMode::Multi: return static_cast<int>(nSignals);
				case OffResMode::Bounded: return 1;
				case OffResMode::MultiBounded: return static_cast<int>(nSignals);
			}
		}
		
		static const vector<string> names(const Components &c) {
			switch (c) {
				case Components::One: return { "PD", "T1", "T2" };
				case Components::Two: return { "PD", "T1_a", "T2_a", "T1_b", "T2_b", "tau_a", "f_a"  };
				case Components::Three: return { "PD", "T1_a", "T2_a", "T1_b", "T2_b", "T1_c", "T2_c", "tau_a", "f_a", "f_c" };
			}
		}
		
		static const ArrayXd &defaultLo(const Components &c, const FieldStrength &tesla) {
			static ArrayXd lo(nP(c));
			switch (tesla) {
				case FieldStrength::Three:
					switch (c) {
						case Components::One: lo << 0., 0.25, 0.01; break;
						case Components::Two: lo << 0., 0.25, 0.01, 0.75, 0.01, 0.01, 0.001; break;
						case Components::Three: lo << 0., 0.35, 0.002, 0.700, 0.075, 3.5, 0.175, 0.05, 0.001, 0.001; break;
					}
				case FieldStrength::Seven:
					switch (c) {
						case Components::One: lo << 0., 0.25, 0.01; break;
						case Components::Two: lo << 0., 0.25, 0.01, 0.75, 0.01, 0.01, 0.001; break;
						case Components::Three: lo << 0., 0.25, 0.01, 0.75, 0.02, 4.00, 0.15, 0., 0., 0.; break;
					}
				case FieldStrength::Unknown: lo.setZero(); break;
			}
			return lo;
		}
		
		static const ArrayXd &defaultHi(const Components &c, const FieldStrength &tesla) {
			static ArrayXd hi(nP(c));
			switch (tesla) {
				case FieldStrength::Three:
					switch (c) {
						case Components::One: hi << 1.e7, 3.0, 0.25; break;
						case Components::Two: hi << 1.e7, 1.0, 0.05, 1.5, 0.05,  0.5, 0.95; break;
						case Components::Three: hi << 1.e7, 0.55, 0.016, 2.0, 0.145,  7.5, 0.5, 0.3, 0.3, 0.95; break;
						default: throw(invalid_argument("Bad number of components."));
					}
				case FieldStrength::Seven:
					switch (c) {
						case Components::One: hi << 1.e7, 5.0, 0.10; break;
						case Components::Two: hi << 1.e7, 1.0, 0.02, 2.0, 0.05,  0.5, 0.95; break;
						case Components::Three: hi << 1.e7, 1.0, 0.02, 2.0, 0.05, 20.0, 0.6, 0.5, 0.95, 0.95; break;
						default: throw(invalid_argument("Bad number of components."));
					}
				case FieldStrength::Unknown: hi.setZero(); break;
			}
			return hi;
		}
	
	protected:
		const Components m_components;
		const OffResMode m_offRes;
		size_t m_nP, m_nV, m_nOffRes;
		vector<DESPOTData> &m_data;
		const bool m_normalise, m_finite, m_debug;
	
	public:
		mcDESPOT(const Components &c, vector<DESPOTData> &data,
				 const OffResMode &offRes, const bool &normalise = false, const bool &finite = false, const bool &debug = false) :
			m_components(c), m_data(data),
			m_normalise(normalise), m_finite(finite), m_offRes(offRes), m_debug(debug)
		{
			m_nP = nP(c);
			m_nOffRes = nOffRes(offRes, m_data.size());
			m_nV = 0;
			for (int i = 0; i < data.size(); i++) {
				if (data[i].flip().size() != data[i].signal().size()) {
					cerr << "Angles and signals size mis-match for signal " << i << endl;
					cerr << "Angles = " << data[i].flip().size() << " signal = " << data[i].signal().size() << endl;
					exit(EXIT_FAILURE);
				}
				m_nV += m_data[i].flip().size();
			}
		}
		
		const bool constraint(const VectorXd &params) {
			// Negative PD or T1/T2 makes no sense
			if ((params[0] < 0.) || (params[1] <= 0.) || (params[2] <= 0.))
				return false;
			
			if (m_components == Components::One) {
				return true;
			} else if (m_components == Components::Two) {
				// Check that T1_a, T2_a < T1_b, T2_b and that f_a makes sense
				if ((params[1] < params[3]) &&
					(params[2] < params[4]) &&
					(params[6] <= 1.0))
					return true;
				else
					return false;
			} else if (m_components == Components::Three) {
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
				
		const long inputs() const { return m_nP + m_nOffRes; }
		const long values() const { return m_nV; }
		const ArrayXd signals() const {
			ArrayXd v(values());
			int index = 0;
			for (int i = 0; i < m_data.size(); i++) {
				v.segment(index, m_data[i].signal().size()) = m_data[i].signal();
				//cout << "s:  " << _signals[i].transpose() << endl;
				index += m_data[i].signal().size();
			}
			return v;
		}
		
		const ArrayXd theory(const VectorXd &params) const {
			ArrayXd t(values());
			int index = 0;
			if (m_debug) cout << endl << "Params: " << params.transpose() << endl;
			for (int i = 0; i < m_data.size(); i++) {
				MagVector M(3, m_data[i].flip().size());
				if ((m_offRes == OffResMode::Single) || (m_offRes == OffResMode::Bounded))
					m_data[i].f0_off = params[m_nP];
				else if ((m_offRes == OffResMode::Multi) || (m_offRes == OffResMode::MultiBounded))
					m_data[i].f0_off = params[m_nP + i];
				if (m_finite) {
					switch (m_components) {
						case Components::One: M = One_SSFP_Finite(m_data[i], params.head(m_nP)); break;
						case Components::Two: M = Two_SSFP_Finite(m_data[i], params.head(m_nP)); break;
						case Components::Three: M = Three_SSFP_Finite(m_data[i], params.head(m_nP)); break;
					}
				} else {
					if (m_data[i].spoil == true) {
						switch (m_components) {
							case Components::One: M = One_SPGR(m_data[i], params.head(m_nP)); break;
							case Components::Two: M = Two_SPGR(m_data[i], params.head(m_nP)); break;
							case Components::Three: M = Three_SPGR(m_data[i], params.head(m_nP)); break;
						}
					} else {
						switch (m_components) {
							case Components::One: M = One_SSFP(m_data[i], params.head(m_nP)); break;
							case Components::Two: M = Two_SSFP(m_data[i], params.head(m_nP)); break;
							case Components::Three: M = Three_SSFP(m_data[i], params.head(m_nP)); break;
						}
					}
				}
				ArrayXd theory = SigMag(M);
				if (m_normalise && (theory.square().sum() > 0.)) theory /= theory.mean();
				t.segment(index, m_data[i].signal().size()) = theory;
				index += m_data[i].signal().size();
			}
			return t;
		}
				
		int operator()(const VectorXd &params, ArrayXd &diffs) const {
			eigen_assert(diffs.size() == values());
			if (m_debug) cout << endl << __PRETTY_FUNCTION__ << endl;
			if (m_debug) cout << "Params: " << params.transpose() << endl;
			ArrayXd t = theory(params);
			ArrayXd s = signals();
			diffs = t - s;
			if (m_debug) cout << "Diffs:  " << diffs.transpose() << endl;
			if (m_debug) cout << "Sum:    " << diffs.square().sum() << endl;
			return 0;
		}
};

//******************************************************************************
#pragma mark DESPOT2FM Functor
//******************************************************************************
class DESPOT2FM : public Functor<double> {
	private:
		long _nV;
		vector<DESPOTData> &_data;
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
		
		DESPOT2FM(vector<DESPOTData> &data, const double T1, const bool &normalise = false, const bool &fitB0 = true) :
			_data(data), _T1(T1), _normalise(normalise), _fitB0(fitB0)
		{
			_nV = 0;
			
			for (int i = 0; i < _data.size(); i++) {
				if (_data[i].flip().size() != _data[i].signal().size()) {
					cerr << "Angles and signals size mis-match for signal " << i << endl;
					cerr << "Angles = " << _data[i].flip().size() << " signal = " << _data[i].signal().size() << endl;
					exit(EXIT_FAILURE);
				}
				_nV += _data[i].flip().size();
			}
		}
		
		const ArrayXd signals() const {
			ArrayXd v(values());
			int index = 0;
			for (int i = 0; i < _data.size(); i++) {
				v.segment(index, _data[i].signal().size()) = _data[i].signal();
				index += _data[i].signal().size();
			}
			return v;
		}
		
		// Params are { PD, T2, B0 }
		const ArrayXd theory(const VectorXd &params) const {
			VectorXd withT1(3);
			withT1 << params[0], _T1, params[1];
			
			ArrayXd t(values());
			int index = 0;
			for (int i = 0; i < _data.size(); i++) {
				if (_fitB0)
					_data[i].f0_off = params[2];
				ArrayXd theory = SigMag(One_SSFP(_data[i], withT1));
				if (_normalise && (theory.square().sum() > 0.)) theory /= theory.mean();
				t.segment(index, _data[i].flip().size()) = theory;
				index += _data[i].flip().size();
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
