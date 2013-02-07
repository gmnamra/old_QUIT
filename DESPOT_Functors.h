/*
 *  DESPOT_Functors.h
 *
 *  Created by Tobias Wood on 16/08/2012.
 *  Copyright (c) 2012 Tobias Wood. All rights reserved.
 */

#ifndef DESPOT_Functors_h
#define DESPOT_Functors_h

#include <vector>
#include <array>
#include <map>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;
using namespace Eigen;

struct DESPOTConstants {
	double TR, phase, B0, B1;
};

// A 3x3 matrix rotation of alpha about X and beta around Z
// Corresponds to RF Flip angle of alpha, and phase-cycling of beta
// Phase-cycling is negative as it rotates the reference frame instead of the
// the magnetic vector, hence is seen as "going backwards".
const Matrix3d RF(const double &alpha, const double &beta)
{
	double ca = cos(alpha), sa = sin(alpha), cb = cos(-beta), sb = sin(-beta);
	Matrix3d R;
	R << cb, ca*sb, sa*sb,
	    -sb, ca*cb, sa*cb,
		  0,   -sa,    ca;
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
// Parameters are { B0, PD, T1, T2 }
//******************************************************************************
VectorXd One_SPGR(const VectorXd &flipAngles,
                  const DESPOTConstants& c, const VectorXd &p)
{
	VectorXd theory(flipAngles.size());
	VectorXd sa = (flipAngles.array() * c.B1).sin();
	VectorXd ca = (flipAngles.array() * c.B1).cos();
	double expT1 = exp(-c.TR / p[2]);
	theory = p[1] * (1. - expT1) * sa.array() / (1. - expT1*ca.array());
	
	return theory;
}

VectorXd One_SSFP(const VectorXd &flipAngles,
                  const DESPOTConstants& c, const VectorXd &p)
{
	Matrix3d A = Relax(p[2], p[3], 0., 2. * M_PI * p[0]),
	         R_rf, eATR;
	Vector3d M0 = Vector3d::Zero(), Mobs;
	M0[2] = p[1]; // PD

	MatrixExponential<Matrix3d> mexp(A * c.TR);
	mexp.compute(eATR);
	const Vector3d RHS = (Matrix3d::Identity() - eATR) * M0;
	
	VectorXd theory(flipAngles.size());
	for (int i = 0; i < flipAngles.size(); i++) {
		const Matrix3d R_rf = RF(c.B1 * flipAngles[i], c.phase);
		Mobs = (Matrix3d::Identity() - (eATR * R_rf)).partialPivLu().solve(RHS);
		theory[i] = Mobs.head(2).norm();
	}
	return theory;
}

//******************************************************************************
#pragma mark Two Component Signals
// Parameters are { B0, PD, T1_a, T2_a, T1_b, T2_b, tau_a, f_a }
//******************************************************************************
typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;

VectorXd Two_SPGR(const VectorXd&flipAngles,
                  const DESPOTConstants &c, const VectorXd &p)
{
	Matrix2d A, eATR;
	Vector2d M0, Mobs;
	VectorXd signal(flipAngles.size());
	double k_ab, k_ba, f_a = p[7], f_b = 1. - f_a;
	CalcExchange(p[6], f_a, f_b, k_ab, k_ba);
	M0 << p[1] * f_a, p[1] * f_b;
	A << -((1./p[2]) + k_ab),                    k_ba,
				        k_ab,      -((1./p[4]) + k_ba);
	MatrixExponential<Matrix2d> expA(A*c.TR);
	expA.compute(eATR);
	const Vector2d RHS = (Matrix2d::Identity() - eATR) * M0;
	for (int i = 0; i < flipAngles.size(); i++) {
		double a = flipAngles[i];
		Mobs = (Matrix2d::Identity() - eATR*cos(c.B1 * a)).partialPivLu().solve(RHS * sin(c.B1 * a));
		signal[i] = Mobs.sum();
	}
	return signal;
}

VectorXd Two_SSFP(const VectorXd&flipAngles, const DESPOTConstants &c, const VectorXd &p, const bool normalise)
{
	VectorXd signal(flipAngles.size());
	Matrix6d A = Matrix6d::Zero(), eATR, eATE, R = Matrix6d::Zero();
	Vector6d M0;
	M0 << 0., 0., p[1] * p[7], 0., 0., p[1] * (1. - p[7]);
	// Set up the 'A' matrix. It's quite complex.
	double TE = c.TR / 2., dw = p[0] * 2. * M_PI, k_ab, k_ba;
	CalcExchange(p[6], p[7], (1 - p[7]), k_ab, k_ba);
	A.block(0, 0, 3, 3) = Relax(1./p[2], 1./p[3], k_ab, dw);
	A.block(3, 3, 3, 3) = Relax(1./p[4], 1./p[5], k_ba, dw);
	A(3, 0) = A(4, 1) = A(5, 2) = k_ab;
	A(0, 3) = A(1, 4) = A(2, 5) = k_ba;
	MatrixExponential<Matrix6d> exp(A*TE);
	exp.compute(eATE);
	eATR.noalias() = eATE * eATE;
	const Vector6d eyemaM0 = (Matrix6d::Identity() - eATR) * M0;
	for (int i = 0; i < flipAngles.size(); i++) {
		const Matrix3d Rb = RF(c.B1 * flipAngles[i], c.phase);
		R.block(0, 0, 3, 3) = Rb;
		R.block(3, 3, 3, 3) = Rb;
		Vector6d MTR = (Matrix6d::Identity() - eATR * R).partialPivLu().solve(eyemaM0);
		if (!normalise) {
			Vector6d Mobs = eATE * R * MTR;
			signal[i] = sqrt(pow(Mobs[0] + Mobs[3], 2.) +
							 pow(Mobs[1] + Mobs[4], 2.));
		} else {
			signal[i] = sqrt(pow(MTR[0] + MTR[3], 2.) +
		                     pow(MTR[1] + MTR[4], 2.));
		}
	}
	return signal;
}

//******************************************************************************
#pragma mark Three Component
// Parameters are { B0, PD, T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c }
//******************************************************************************
typedef Matrix<double, 9, 9> Matrix9d;
typedef Matrix<double, 9, 1> Vector9d;

VectorXd Three_SPGR(const VectorXd&flipAngles, const DESPOTConstants &c, const VectorXd &p)
{
	Matrix3d A, eATR;
	Vector3d M0, Mobs;
	VectorXd signal(flipAngles.size());
	double k_ab, k_ba, f_a = p[9], f_c = p[10], f_b = 1. - f_a - f_c;
	CalcExchange(p[8], f_a, f_b, k_ab, k_ba);
	M0 << p[1] * f_a, p[1] * f_b, p[1] * f_c;
	A << -((1./p[2]) + k_ab),                     k_ba,         0.,
					    k_ab,      -((1./p[4]) + k_ba),         0.,
						  0.,                       0., -(1./p[6]);
	MatrixExponential<Matrix3d> expA(A*c.TR);
	expA.compute(eATR);
	const Vector3d RHS = (Matrix3d::Identity() - eATR) * M0;
	for (int i = 0; i < flipAngles.size(); i++) {
		double a = flipAngles[i];
		Mobs = (Matrix3d::Identity() - eATR*cos(c.B1 * a)).partialPivLu().solve(RHS * sin(c.B1 * a));
		signal[i] = Mobs.sum();
	}
	return signal;
}

VectorXd Three_SSFP(const VectorXd&flipAngles, const DESPOTConstants &c, const VectorXd &p,
                    const bool normalise)
{
	VectorXd signal(flipAngles.size());
	Matrix9d A = Matrix9d::Zero(), eATR, eATE, R = Matrix9d::Zero();
	Vector9d M0;
	double TE = c.TR / 2., dw = p[0] * 2. * M_PI, k_ab, k_ba,
	       f_a = p[9], f_c = p[10], f_b = 1. - f_a - f_c;
	M0 << 0., 0., p[1] * f_a, 0., 0., p[1] * f_b, 0., 0., p[1] * f_c;
	// Set up the 'A' matrix. It's quite complex.

	CalcExchange(p[8], f_a, f_b, k_ab, k_ba);
	A.block(0, 0, 3, 3) = Relax(1./p[2], 1./p[3], k_ab, dw);
	A.block(3, 3, 3, 3) = Relax(1./p[4], 1./p[5], k_ba, dw);
	A.block(6, 6, 3, 3) = Relax(1./p[6], 1./p[7], 0., dw);
	A(3, 0) = A(4, 1) = A(5, 2) = k_ab;
	A(0, 3) = A(1, 4) = A(2, 5) = k_ba;
	MatrixExponential<Matrix9d> exp(A*TE);
	exp.compute(eATE);
	eATR.noalias() = eATE * eATE;
	const Vector9d eyemaM0 = (Matrix9d::Identity() - eATR) * M0;
	for (int i = 0; i < flipAngles.size(); i++) {
		const Matrix3d Rb = RF(c.B1 * flipAngles[i], c.phase);
		R.block(0, 0, 3, 3) = Rb;
		R.block(3, 3, 3, 3) = Rb;
		R.block(6, 6, 3, 3) = Rb;
		Vector9d MTR = (Matrix9d::Identity() - eATR * R).partialPivLu().solve(eyemaM0);
		if (!normalise) {
			Vector9d Mobs = eATE * R * MTR;
			signal[i] = sqrt(pow(Mobs[0] + Mobs[3] + Mobs[6], 2.) +
							 pow(Mobs[1] + Mobs[4] + Mobs[7], 2.));
		} else {
			signal[i] = sqrt(pow(MTR[0] + MTR[3] + MTR[6], 2.) +
		                     pow(MTR[1] + MTR[4] + MTR[7], 2.));
		}
	}
	return signal;
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
		
		virtual int operator()(const VectorXd &params, VectorXd &diffs) const = 0;
		
		virtual const VectorXd theory(const VectorXd &params) const = 0;
		virtual const VectorXd signals() const = 0;
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
	
	private:
		const int _components;
		long _nP, _nV;
		const vector<SignalType> &_types;
		const vector<VectorXd> &_angles, &_signals;
		const vector<DESPOTConstants> &_consts;
		const bool &_normalise, &_fitB0;
	
	public:
		static const int nP(const int components) {
			switch (components) {
				case 1: return 4;
				case 2: return 8;
				case 3: return 11;
				default:
					std::cerr << "Cannot create a " << components << "-component mcDESPOT model." << std::endl;
					exit(EXIT_FAILURE);
			}
		}
		
		static const vector<const string> &names(const int components) {
			static const map<int, const vector<const string> > _namesMap {
				{1, { "1c_B0", "1c_PD", "1c_T1", "1c_T2" } },
				{2, { "2c_B0", "2c_PD", "2c_T1_a", "2c_T2_a", "2c_T1_b", "2c_T2_b", "2c_tau_a", "2c_f_a"  } },
				{3, { "3c_B0", "3c_PD", "3c_T1_a", "3c_T2_a", "3c_T1_b", "3c_T2_b", "3c_T1_c", "3c_T2_c", "3c_tau_a", "3c_f_a", "3c_f_c" } } };
			map<int, const vector<const string> >::const_iterator it = _namesMap.find(components);
			if (it == _namesMap.end()) {
				std::cerr << "Don't have file names for a " << components << "-component mcDESPOT model." << std::endl;
				exit(EXIT_FAILURE);
			} else
			return it->second;
		}
		
		static const ArrayXd &defaultLo(const int components, const int tesla) {
			static ArrayXd c1t3(4), c1t7(4), c2t3(8), c2t7(8), c3t3(11), c3t7(11);
			c1t3 << -150., 0., 0.25, 0.01;
			c1t7 << -150., 0., 0.25, 0.01;
			c2t3 << -150., 0., 0.25, 0.01, 0.75, 0.01, 0.01, 0.0;
			c2t7 << -150., 0., 0.25, 0.01, 0.75, 0.01, 0.01, 0.0;
			c3t3 << -150., 0., 0.25, 0.01, 0.75, 0.02, 2.00, 0.15, 0., 0., 0.;
			c3t7 << -150., 0., 0.25, 0.01, 0.75, 0.02, 4.00, 0.15, 0., 0., 0.;
			
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
		
		static const ArrayXd defaultHi(const int components, const int tesla) {
			static ArrayXd c1t3(4), c1t7(4), c2t3(8), c2t7(8), c3t3(11), c3t7(11);
			c1t3 << 150., 1.e7, 3.0, 0.25;
			c1t7 << 150., 1.e7, 5.0, 0.10;
			c2t3 << 150., 1.e7, 1.0, 0.05, 1.5, 0.05,  0.5, 0.95;
			c2t7 << 150., 1.e7, 1.0, 0.02, 2.0, 0.05,  0.5, 0.95;
			c3t3 << 150., 1.e7, 1.0, 0.05, 1.5, 0.05,  7.5, 0.6, 0.5, 0.49, 0.95;
			c3t7 << 150., 1.e7, 1.0, 0.02, 2.0, 0.05, 20.0, 0.6, 0.5, 0.95, 0.95;
			
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
		
		
		const bool constraint(const VectorXd &params) {
			// Negative PD or T1/T2 makes no sense
			if ((params[1] < 0.) || (params[2] <= 0.) || (params[3] <= 0.))
				return false;
			
			if (_components == 1) {
				return true;
			} else if (_components == 2) {
				// Check that T1_a, T2_a < T1_b, T2_b and that f_a makes sense
				if ((params[4] > 0.) && (params[5] > 0.) &&
					(params[2] < params[4]) &&
					(params[3] < params[5]) &&
					(params[7] <= 1.0))
					return true;
				else
					return false;
			} else if (_components == 3) {
				if ((params[4] > 0.) && (params[5] > 0.) && (params[6] > 0.) && (params[7] > 0.) &&
					(params[2] < params[4]) &&
					(params[3] < params[5]) &&
					(params[4] < params[6]) &&
					(params[5] < params[7]) &&
					((params[9] + params[10]) <= 0.95))
					return true;
				else
					return false;
			} else
				return true;
		}
				
		const long inputs() const { return _nP; }
		const long values() const { return _nV; }
		
		mcDESPOT(const int components, const vector<SignalType> &types,
				 const vector<VectorXd> &angles, const vector<VectorXd> &signals,
		         const vector<DESPOTConstants> &constants,
				 const bool &normalise = false, const bool &fitB0 = true) :
			_components(components), _types(types),
			_angles(angles), _signals(signals), _consts(constants),
			_normalise(normalise), _fitB0(fitB0)
		{
			_nP = nP(components);
			_nV = 0;
			
			for (int i = 0; i < angles.size(); i++) {
				if (angles[i].size() != signals[i].size()) {
					std::cerr << "Angles and signals size mis-match for signal " << i << std::endl;
					std::cerr << "Angles = " << angles[i].size() << " signal = " << signals[i].size() << std::endl;
					exit(EXIT_FAILURE);
				}
				_nV += angles[i].size();
			}
		}
		
		const VectorXd signals() const {
			VectorXd v(values());
			int index = 0;
			for (int i = 0; i < _signals.size(); i++) {
				v.segment(index, _signals[i].size()) = _signals[i];
				index += _signals[i].size();
			}
			return v;
		}
		
		const VectorXd theory(const VectorXd &params) const {
			VectorXd localP = params;
			VectorXd t(values());
			int index = 0;
			for (int i = 0; i < _signals.size(); i++) {
				VectorXd theory(_signals[i].size());
				if (!_fitB0)
					localP[0] = _consts[i].B0;
				if (_types[i] == SignalSPGR) {
					switch (_components) {
						case 1: theory = One_SPGR(_angles[i], _consts[i], localP); break;
						case 2: theory = Two_SPGR(_angles[i], _consts[i], localP); break;
						case 3: theory = Three_SPGR(_angles[i], _consts[i], localP); break;
					}
				} else if (_types[i] == SignalSSFP) {
					switch (_components) {
						case 1: theory = One_SSFP(_angles[i], _consts[i], localP); break;
						case 2: theory = Two_SSFP(_angles[i], _consts[i], localP, _normalise); break;
						case 3: theory = Three_SSFP(_angles[i], _consts[i], localP, _normalise); break;
					}
				}
				if (_normalise && (theory.norm() > 0.)) theory /= theory.mean();
				t.segment(index, _signals[i].size()) = theory;
				index += _signals[i].size();
			}
			return t;
		}
				
		int operator()(const VectorXd &params, VectorXd &diffs) const {
			eigen_assert(diffs.size() == values());
			VectorXd t = theory(params);
			VectorXd s = signals();
			diffs = t - s;
			return 0;
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
		const vector<DESPOTConstants> &_consts;
		const double _T1;
		const bool &_normalise, &_fitB0;
	
	public:
		static const vector<const string> &names() {
			static const vector<const string> _names { { "FM_B0", "FM_PD", "FM_T2" } };
			return _names;
		}
		
		static const ArrayXd &defaultLo(const int tesla) {
			static ArrayXd t3(3), t7(3);
			t3 << 0.,   0.010, -150.;
			t7 << 0.,   0.005, -150.;
			
			switch (tesla) {
				case 3: return t3; break;
				case 7: return t7; break;
			}
			std::cerr << "Don't have defaults for " << tesla << "tesla." << std::endl;
			exit(EXIT_FAILURE);
		}
		
		static const ArrayXd defaultHi(const int tesla) {
			static ArrayXd t3(3), t7(3);
			t3 << 1.e6,   1.00, 150.;
			t7 << 1.e6,   0.25, 150.;
			
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
		
		static const long inputs() { return 3; }
		const long values() const { return _nV; }
		
		DESPOT2FM(VectorXd &angles, const vector<VectorXd> &signals,
		          const vector<DESPOTConstants> &consts,
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
		
		const VectorXd signals() const {
			VectorXd v(values());
			int index = 0;
			for (int i = 0; i < _signals.size(); i++) {
				v.segment(index, _signals[i].size()) = _signals[i];
				index += _signals[i].size();
			}
			return v;
		}
		
		// Params are { B0, PD, T2 }
		const VectorXd theory(const VectorXd &params) const {
			VectorXd withT1(4);
			withT1 << params[0], params[1], _T1, params[2];
			
			VectorXd t(values());
			int index = 0;
			for (int i = 0; i < _signals.size(); i++) {
				VectorXd theory(_signals[i].size());
				if (!_fitB0)
					withT1[0] = _consts[i].B0;
				theory = One_SSFP(_angles, _consts[i], withT1);
				if (_normalise && (theory.norm() > 0.)) theory /= theory.mean();
				t.segment(index, _signals[i].size()) = theory;
				index += _signals[i].size();
			}
			return t;
		}
				
		int operator()(const VectorXd &params, VectorXd &diffs) const {
			eigen_assert(diffs.size() == values());
			VectorXd t = theory(params);
			VectorXd s = signals();
			diffs = t - s;
			return 0;
		}
};

#endif
