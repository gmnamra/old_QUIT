//
//  DESPOT_Functors.h
//  DESPOT
//
//  Created by Tobias Wood on 16/08/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#ifndef DESPOT_Functors_h
#define DESPOT_Functors_h

#include <vector>
#include <array>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "NiftiImage.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
#pragma mark Functor Base Classes
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
		
		virtual VectorXd signals() = 0;
};

// DESPOT Base Functor
template<int nPt=Dynamic>
class DESPOT_Functor : public Functor<double, nPt>
{
	protected:
		const VectorXd &_spgrAngles, &_spgrB1,
					   &_ssfpAngles, &_ssfpB0, &_ssfpB1;
		const vector<VectorXd> &_spgrSignals, &_ssfpSignals;
		const vector<double> &_ssfpPhases;
		double _spgrTR, _ssfpTR;
		const bool _normalise, _fitB0;
			
	public:
		static const int nP = nPt;
		static const array<string, nPt> names;
		static const array<double, nPt> lo3Bounds, hi3Bounds, lo7Bounds, hi7Bounds;
		
		
		const bool constraint(const VectorXd &params) const { return true; }

		DESPOT_Functor(const VectorXd &spgrAngles, const vector<VectorXd> &spgrSignals,
					   const VectorXd &spgrB1, const double &spgrTR,
		               const VectorXd &ssfpAngles, const vector<double> &ssfpPhases, const vector<VectorXd> &ssfpSignals,
					   const VectorXd &ssfpB0, const VectorXd &ssfpB1, const double &ssfpTR,
					   const bool &normalise = false, const bool &fitB0 = true) :
					   Functor<double, nP>(spgrAngles.size() * spgrSignals.size() + ssfpAngles.size() * ssfpSignals.size()),
				       _spgrAngles(spgrAngles), _spgrSignals(spgrSignals), _spgrB1(spgrB1),
				       _ssfpAngles(ssfpAngles), _ssfpPhases(ssfpPhases), _ssfpSignals(ssfpSignals),
					   _ssfpB0(ssfpB0), _ssfpB1(ssfpB1),
					   _spgrTR(spgrTR), _ssfpTR(ssfpTR),
					   _normalise(normalise), _fitB0(fitB0)
				      {
						//std::cout << "SPGR B1 " << _spgrB1.transpose() << " SSFP B0 " << _ssfpB0.transpose() << " B1 " << _ssfpB1.transpose() << std::endl;
					  }
		
		VectorXd signals()
		{
			VectorXd v(this->values());
			int index = 0;
			for (int i = 0; i < _spgrSignals.size(); i++)
			{
				v.segment(index, _spgrAngles.size()) = _spgrSignals[i];
				index += _spgrAngles.size();
			}
			for (int i = 0; i < _ssfpSignals.size(); i++)
			{
				v.segment(index, _ssfpAngles.size()) = _ssfpSignals[i];
				index += _ssfpAngles.size();
			}
			return v;
		}
};
		

//******************************************************************************
#pragma mark One Component Signals
//******************************************************************************
VectorXd One_SPGR(const VectorXd&flipAngles,
                  const double TR, const double M0, const double T1,
				  const double B1)
{
	VectorXd theory(flipAngles.size());
	VectorXd sa = (flipAngles.array() * B1).sin();
	VectorXd ca = (flipAngles.array() * B1).cos();
	double expT1 = exp(-TR / T1);
	theory = M0 * (1. - expT1) * sa.array() / (1. - expT1*ca.array());
	
	return theory;
}

VectorXd One_SSFP(const VectorXd &flipAngles, const double rfPhase,
                  const double TR, const double M0,
				  const double T1, const double T2,
				  const double B0, const double B1)
{
	Matrix3d A = Matrix3d::Zero(), R_rf = Matrix3d::Zero(),
	                eye = Matrix3d::Identity(), expA, eyemA;
	Vector3d M = Vector3d::Zero(), Mobs;
	M[2] = M0;
	R_rf(0, 0) = 1.;
	double phase = rfPhase / TR + (B0 * 2. * M_PI);
	A(0, 0) = A(1, 1) = -1 / T2;
	A(0, 1) =  phase;
	A(1, 0) = -phase;
	A(2, 2) = -1 / T1;
	MatrixExponential<Matrix3d> expmA(A*TR);
	expmA.compute(expA);
	eyemA.noalias() = eye - expA;
	
	VectorXd theory(flipAngles.size());
	for (int i = 0; i < flipAngles.size(); i++)
	{
		double a = flipAngles[i];
		double ca = cos(B1 * a), sa = sin(B1 * a);
		R_rf(1, 1) = R_rf(2, 2) = ca;
		R_rf(1, 2) = sa; R_rf(2, 1) = -sa;
		Mobs = (eye - (expA * R_rf)).partialPivLu().solve(eyemA) * M;				
		theory[i] = Mobs.head(2).norm();
	}
	return theory;
}
//******************************************************************************
#pragma mark OneComponent Functor
//******************************************************************************
class OneComponent : public DESPOT_Functor<4>
{
	public:
		OneComponent(const VectorXd &spgrAngles, const vector<VectorXd> &spgrSignals,
					 const VectorXd &spgrB1, const double &spgrTR,
					 const VectorXd &ssfpAngles, const vector<double> &ssfpPhases, const vector<VectorXd> &ssfpSignals,
					 const VectorXd &ssfpB0, const VectorXd &ssfpB1, const double &ssfpTR,
					 const bool &normalise = false, const bool &fitB0 = true)
					: DESPOT_Functor(spgrAngles, spgrSignals, spgrB1, spgrTR,
					                 ssfpAngles, ssfpPhases, ssfpSignals, ssfpB0, ssfpB1, ssfpTR,
					                 normalise, true)
					{}
		static const array<string, nP> names;
		static const array<double, nP> lo3Bounds, hi3Bounds, lo7Bounds, hi7Bounds;
		
		int operator()(const VectorXd &params, VectorXd &diffs) const
		{
			double PD = params[0], T1 = params[1], T2 = params[2], B0 = params[3];
			eigen_assert(diffs.size() == values());
			int index = 0;
			for (int i = 0; i < _spgrSignals.size(); i++)
			{
				VectorXd theory = One_SPGR(_spgrAngles, _spgrTR, PD, T1, _spgrB1[i]);
				//std::cout << "SPGR theory " << theory.transpose() << std::endl;
				if (_normalise) theory /= theory.mean();
				//std::cout << "SPGR normalise " << theory.transpose() << std::endl;
				//std::cout << "SPGR signal " << _spgrSignals[i].transpose() << std::endl;
				diffs.segment(index, _spgrAngles.size()) = theory - _spgrSignals[i];
				index += _spgrAngles.size();
			}

			for (int i = 0; i < _ssfpSignals.size(); i++)
			{
				VectorXd theory = One_SSFP(_ssfpAngles, _ssfpPhases[i], _ssfpTR,
										   PD, T1, T2, B0, _ssfpB1[i]);
				//std::cout << "SSFP theory " << theory.transpose() << std::endl;
				if (_normalise) theory /= theory.mean();
				//std::cout << "SSFP normalise " << theory.transpose() << std::endl;
				//std::cout << "SSFP signal " << _ssfpSignals[i].transpose() << std::endl;
				diffs.segment(index, _ssfpAngles.size()) = theory - _ssfpSignals[i];
				index += _ssfpAngles.size();
			}
			//std::cout << "diffs " << diffs.transpose() << std::endl;
			return 0;
		}
};
const array<string, OneComponent::nP> OneComponent::names{ { string("M0"), "T1", "T2", "B0" } };
const array<double, OneComponent::nP> OneComponent::lo3Bounds{ { 0.,   0.1, 0.01, -150. } };
const array<double, OneComponent::nP> OneComponent::hi3Bounds{ { 1.e7, 3.0, 1.0, 150. } };
const array<double, OneComponent::nP> OneComponent::lo7Bounds{ { 0.,   0.1, 0.005, -150. } };
const array<double, OneComponent::nP> OneComponent::hi7Bounds{ { 1.e7, 5.0, 0.05, 150. } };

//******************************************************************************
#pragma mark OneComponentSSFP
//******************************************************************************

class OneComponentSSFP : public Functor<double>
{
	public:
		const VectorXd &_flipAngles, &_rfPhases;
		const vector<VectorXd> &_signals;
		const double _TR, _T1, _M0, _B0, _B1;
		const bool _normalise;
		
		static const int nP = 1;
		static const char *names[];
		
		OneComponentSSFP(const VectorXd &flipAngles, const VectorXd &rfPhases,
		                 const vector<VectorXd> &signals,
				         const double TR, const double M0, const double T1,
						 const double B0, const double B1, const bool normalise = false) :
				         Functor<double>(OneComponentSSFP::nP, flipAngles.size() * rfPhases.size()),
				         _flipAngles(flipAngles),
				         _rfPhases(rfPhases),
				         _signals(signals),
				         _TR(TR), _T1(T1), _M0(M0), _B0(B0), _B1(B1),
						 _normalise(normalise)
				         {}
			
		int operator()(const VectorXd &params, VectorXd &diffs) const
		{
			eigen_assert(diffs.size() == values());
			int index = 0;
			for (int p = 0; p < _rfPhases.size(); p++)
			{
				VectorXd temp = One_SSFP(_flipAngles, _rfPhases(p), _TR, _M0, _T1, params[0], _B0, _B1);
				if (_normalise)
					temp /= temp.mean();
				diffs.segment(index, _flipAngles.size()) = temp - _signals[p];
				index += _flipAngles.size();
			}
			return 0;
		}
};
const char *OneComponentSSFP::names[] = { "T2" };

//******************************************************************************
#pragma mark Two Component Signals
//******************************************************************************
typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;

VectorXd Two_SPGR(const VectorXd&flipAngles,
                  const double &TR, const double &PD, const double &B1,
				  const double &T1_a, const double &T1_b,
				  const double &f_a, const double &f_b,
				  const double &k_ab, const double &k_ba)
{
	Matrix2d A;
	const Matrix2d eye2 = Matrix2d::Identity();
	Vector2d M0, Mobs;
	VectorXd signal(flipAngles.size());
	
	M0 << PD * f_a, PD * f_b;
	A << -(TR/T1_a + TR*k_ab),             TR*k_ba,
					  TR*k_ab, -(TR/T1_b + TR*k_ba);
	//cout << A << endl;
	MatrixExponential<Matrix2d> expA(A);
	expA.compute(A);
	const Matrix2d eyema = eye2 - A;
	for (int i = 0; i < flipAngles.size(); i++)
	{
		double a = flipAngles[i];
		Mobs = (eye2 - A*cos(B1 * a)).partialPivLu().solve(eyema * sin(B1 * a)) * M0;
		signal[i] = Mobs.sum();
	}
	
	return signal;
}

VectorXd Two_SSFP(const VectorXd&flipAngles, const double &rfPhase,
                  const double &TR, const double &PD,
				  const double &B0, const double &B1,
				  const double &T1_a, const double &T1_b,
				  const double &T2_a, const double &T2_b,
				  const double &f_a, const double &f_b,
				  const double &k_ab, const double &k_ba)
{
	Matrix6d A, expA, R_rf, eye_mAR;
	const Matrix6d eye6 = Matrix6d::Identity();
	Vector6d M0, Mobs;
	PartialPivLU<Matrix6d> solver;
	A.setZero(); R_rf.setZero();
	R_rf(0, 0) = R_rf(1, 1) = 1.;
	M0 << 0., 0., 0., 0., PD * f_a, PD * f_b;

	VectorXd signal(flipAngles.size());
	double phase = rfPhase + (B0 * TR * 2. * M_PI);
	// Can get away with this because the block structure of the
	// matrix ensures that the zero blocks are always zero after
	// the matrix exponential.
	A(0, 0) = A(2, 2) = -TR * (1./T2_a + k_ab);
	A(1, 1) = A(3, 3) = -TR * (1./T2_b + k_ba);
	A(0, 1) = A(2, 3) = A(4, 5) = TR * k_ba;
	A(0, 2) = A(1, 3) = phase;
	A(1, 0) = A(3, 2) = A(5, 4) = TR * k_ab;
	A(2, 0) = A(3, 1) = -phase; 
	A(4, 4) = -TR * (1./T1_a + k_ab);
	A(5, 5) = -TR * (1./T1_b + k_ba);
	A(0, 3) = A(1, 2) = A(2, 1) = A(3, 0) = 0.;
	MatrixExponential<Matrix6d> exp(A);
	exp.compute(expA);
	const Matrix6d eyema = eye6 - expA;
	for (int i = 0; i < flipAngles.size(); i++)
	{
		double a = flipAngles[i];
		double ca = cos(B1 * a), sa = sin(B1 * a);
		R_rf(2, 2) = R_rf(3, 3) = R_rf(4, 4) = R_rf(5, 5) = ca;
		R_rf(2, 4) = R_rf(3, 5) = sa;
		R_rf(4, 2) = R_rf(5, 3) = -sa;
		eye_mAR.noalias() = eye6 - (expA * R_rf);
		solver.compute(eye_mAR);
		Mobs.noalias() = solver.solve(eyema) * M0;
		signal[i] = sqrt(pow(Mobs[0] + Mobs[1], 2.) +
						  pow(Mobs[2] + Mobs[3], 2.));
	}
	return signal;
}

//******************************************************************************
#pragma mark Two Component Functor
//******************************************************************************
class TwoComponent : public DESPOT_Functor<8>
{
	public:
		TwoComponent(const VectorXd &spgrAngles, const vector<VectorXd> &spgrSignals,
					 const VectorXd &spgrB1, const double &spgrTR,
					 const VectorXd &ssfpAngles, const vector<double> &ssfpPhases, const vector<VectorXd> &ssfpSignals,
					 const VectorXd &ssfpB0, const VectorXd &ssfpB1, const double &ssfpTR,
					 const bool &normalise = false, const bool &fitB0 = true)
					: DESPOT_Functor(spgrAngles, spgrSignals, spgrB1, spgrTR,
					                 ssfpAngles, ssfpPhases, ssfpSignals, ssfpB0, ssfpB1, ssfpTR,
					                 normalise, true)
					{}
		static const array<string, nP> names;
		static const array<double, nP> lo3Bounds, hi3Bounds, lo7Bounds, hi7Bounds;
		
		static bool constraint(const VectorXd &params)
		{
			if ((params[1] < params[2]) &&
				(params[3] < params[4]) &&
				(params[5] <= 1.0))
				return true;
			else
				return false;

		}
	
		int operator()(const VectorXd &params, VectorXd &diffs) const
		{
			eigen_assert(diffs.size() == values());
			double PD   = params[0],
			       T1_a = params[1],   T1_b = params[2],
				   T2_a = params[3],   T2_b = params[4],
			       f_a  = params[5],    f_b = 1. - f_a,
			       tau_a = params[6], tau_b = f_b * tau_a / f_a,
				   k_ab = 1. / tau_a,  k_ba = 1. / tau_b,
				   B0 = params[7];
			// Only have 1 component, so no exchange
			if ((f_a == 0.) || (f_b == 0.))
			{
				k_ab = 0.;
				k_ba = 0.;
			}
			int index = 0;
			for (int i = 0; i < _spgrSignals.size(); i++)
			{
				VectorXd theory = Two_SPGR(_spgrAngles, _spgrTR, PD, _spgrB1[i],
				                           T1_a, T1_b, f_a, f_b, k_ab, k_ba);
				if (_normalise) theory /= theory.mean();
				diffs.segment(index, _spgrAngles.size()) = theory - _spgrSignals[i];
				index += _spgrAngles.size();
			}

			for (int i = 0; i < _ssfpSignals.size(); i++)
			{
				VectorXd theory = Two_SSFP(_ssfpAngles, _ssfpPhases[i], _ssfpTR,
				                           PD, B0, _ssfpB1[i],
				                           T1_a, T1_b, T2_a, T2_b,
										   f_a, f_b, k_ab, k_ba);
				if (_normalise) theory /= theory.mean();
				diffs.segment(index, _ssfpAngles.size()) = theory - _ssfpSignals[i];
				index += _ssfpAngles.size();
			}
			return 0;
		}
};
const array<string, TwoComponent::nP> TwoComponent::names{ { "M0", "T1_a", "T1_b", "T2_a", "T2_b", "f_a", "tau_a", "B0" } };
const array<double, TwoComponent::nP> TwoComponent::lo3Bounds{ { 0., 0.1, 0.8, 0.001, 0.01, 0.0, 0.05, 0. } };
const array<double, TwoComponent::nP> TwoComponent::hi3Bounds{ { 1.e7, 1.0, 3.0, 0.050, 0.25, 1.0, 2.00, 0. } };
const array<double, TwoComponent::nP> TwoComponent::lo7Bounds{ { 0., 0.1, 0.8, 0.001, 0.01, 0.0, 0.05, 0. } };
const array<double, TwoComponent::nP> TwoComponent::hi7Bounds{ { 1.e7, 1.0, 3.0, 0.050, 0.25, 1.0, 2.00, 0. } };

//******************************************************************************
#pragma mark Three Component
//******************************************************************************
typedef Matrix<double, 9, 9> Matrix9d;
typedef Matrix<double, 9, 1> Vector9d;

VectorXd Three_SPGR(const VectorXd&flipAngles,
                    const double &TR, const double &PD, const double &B1,
				    const double &T1_a, const double &T1_b, const double &T1_c,
				    const double &f_a, const double &f_b, const double &f_c,
				    const double &k_ab, const double &k_ba)
{
	VectorXd signal(flipAngles.size());
	Matrix3d A;
	const Matrix3d eye3 = Matrix3d::Identity();
	Vector3d M0, Mobs;
	VectorXd signals(flipAngles.size());
	M0 << PD * f_a, PD * f_b, PD * f_c;
	A << -(1./T1_a + k_ab),               k_ba, 0,
					  k_ab,  -(1./T1_b + k_ba), 0,
						 0,                  0, -1./T1_c;
	A *= TR;
	MatrixExponential<Matrix3d> expA(A);
	expA.compute(A);
	const Matrix3d eyema = eye3 - A;
	for (int i = 0; i < flipAngles.size(); i++)
	{
		double a = flipAngles[i];
		Mobs = (eye3 - A*cos(B1 * a)).partialPivLu().solve(eyema * sin(B1 * a)) * M0;
		signal[i] = Mobs.sum();
	}
	return signal;
}

VectorXd Three_SSFP(const VectorXd&flipAngles, const double &rfPhase,
                    const double &TR, const double &PD,
				    const double &B0, const double &B1,
				    const double &T1_a, const double &T1_b, const double &T1_c,
				    const double &T2_a, const double &T2_b, const double &T2_c,
				    const double &f_a, const double &f_b, const double &f_c,
				    const double &k_ab, const double &k_ba)
{
	VectorXd signal(flipAngles.size());
	Matrix9d A = Matrix9d::Zero(), expA, R_rf = Matrix9d::Zero(), eye_mAR;
	const Matrix9d eye9 = Matrix9d::Identity();
	Vector9d M0, Mobs;
	PartialPivLU<Matrix9d> solver;
	// Set up the 'A' matrix. It's quite complex.
	A(0, 0) = A(3, 3) = -TR * (1./T2_a + k_ab);
	A(1, 1) = A(4, 4) = -TR * (1./T2_b + k_ba);
	A(2, 2) = A(5, 5) = -TR * (1./T2_c);
	A(0, 1) = A(3, 4) = A(6, 7) = TR * k_ba;
	A(1, 0) = A(4, 3) = A(7, 6) = TR * k_ab;
	A(6, 6) = -TR * (1./T1_a + k_ab);
	A(7, 7) = -TR * (1./T1_b + k_ba);
	A(8, 8) = -TR * (1./T1_c);
	R_rf(0, 0) = R_rf(1, 1) = R_rf(2, 2) = 1.;
	M0 << 0., 0., 0., 0., 0., 0., PD * f_a, PD * f_b, PD * f_c;
	
	double phase = rfPhase + (B0 * TR * 2. * M_PI);
	A(0, 3) = A(1, 4) = A(2, 5) = phase;
	A(3, 0) = A(4, 1) = A(5, 2) = -phase;
	MatrixExponential<Matrix9d> exp(A);
	exp.compute(expA);
	const Matrix9d eyema = eye9 - expA;
	for (int i = 0; i < flipAngles.size(); i++)
	{
		double a = flipAngles[i];
		double ca = cos(B1 * a), sa = sin(B1 * a);
		R_rf(3, 3) = R_rf(4, 4) = R_rf(5, 5) =
		R_rf(6, 6) = R_rf(7, 7) = R_rf(8, 8) =  ca;
		R_rf(3, 6) = R_rf(4, 7) = R_rf(5, 8) =  sa;
		R_rf(6, 3) = R_rf(7, 4) = R_rf(8, 5) = -sa;
		eye_mAR.noalias() = eye9 - (expA * R_rf);
		solver.compute(eye_mAR);
		Mobs.noalias() = solver.solve(eyema) * M0;
		signal[i] = sqrt(pow(Mobs[0] + Mobs[1] + Mobs[2], 2.) +
						 pow(Mobs[3] + Mobs[4] + Mobs[5], 2.));
	}
	return signal;
}

class ThreeComponent : public DESPOT_Functor<11>
{
	public:
		ThreeComponent(const VectorXd &spgrAngles, const vector<VectorXd> &spgrSignals,
					   const VectorXd &spgrB1, const double &spgrTR,
					   const VectorXd &ssfpAngles, const vector<double> &ssfpPhases, const vector<VectorXd> &ssfpSignals,
					   const VectorXd &ssfpB0, const VectorXd &ssfpB1, const double &ssfpTR,
					   const bool &normalise = false, const bool &fitB0 = true)
					  : DESPOT_Functor(spgrAngles, spgrSignals, spgrB1, spgrTR,
					                 ssfpAngles, ssfpPhases, ssfpSignals, ssfpB0, ssfpB1, ssfpTR,
					                 normalise, true)
					  {}
		static const array<string, nP> names;
		static const array<double, nP> lo3Bounds, hi3Bounds, lo7Bounds, hi7Bounds;
		
		static bool constraint(const VectorXd &params)
		{
			if ((params[1] < params[2]) &&
			    (params[2] < params[3]) &&
				(params[4] < params[5]) &&
				(params[5] < params[6]) &&
			    ((params[7] + params[8]) <= 0.95))
				return true;
			else
				return false;
		}
		
		int operator()(const VectorXd &params, VectorXd &diffs) const
		{
			eigen_assert(diffs.size() == values());
			double PD   = params[0],
			       T1_a = params[1], T1_b = params[2], T1_c = params[3],
				   T2_a = params[4], T2_b = params[5], T2_c = params[6],
			       f_a  = params[7], f_c  = params[8], f_b  = 1. - f_a - f_c,
			       tau_a = params[9], tau_b = f_b * tau_a / f_a,
				   k_ab = 1. / tau_a, k_ba = 1. / tau_b,
				   B0 = params[10];
			// Only have 1 component, so no exchange
			if ((f_a == 0.) || (f_b == 0.))
			{
				k_ab = 0.;
				k_ba = 0.;
			}
			int index = 0;
			for (int i = 0; i < _spgrSignals.size(); i++)
			{
				VectorXd theory = Three_SPGR(_spgrAngles, _spgrTR, PD, _spgrB1[i],
				                           T1_a, T1_b, T1_c, f_a, f_b, f_c, k_ab, k_ba);
				if (_normalise) theory /= theory.mean();
				diffs.segment(index, _spgrAngles.size()) = theory - _spgrSignals[i];
				index += _spgrAngles.size();
			}

			for (int i = 0; i < _ssfpSignals.size(); i++)
			{
				VectorXd theory = Three_SSFP(_ssfpAngles, _ssfpPhases[i], _ssfpTR,
				                           PD, B0, _ssfpB1[i],
				                           T1_a, T1_b, T1_c, T2_a, T2_b, T2_c,
										   f_a, f_b, f_c, k_ab, k_ba);
				if (_normalise) theory /= theory.mean();
				diffs.segment(index, _ssfpAngles.size()) = theory - _ssfpSignals[i];
				index += _ssfpAngles.size();
			}
			return 0;
		}
};
const array<string, ThreeComponent::nP> ThreeComponent::names{ { "M0", "3c_T1_a", "3c_T1_b", "3c_T1_c", "3c_T2_a", "3c_T2_b", "3c_T2_c", "3c_f_a", "3c_f_c", "3c_tau_a", "B0" } };
const array<double, ThreeComponent::nP> ThreeComponent::lo3Bounds{ { 0., 0.250, 0.250, 1.500, 0.000, 0.000, 0.150, 0.00, 0.00, 0.025, 0. } };
const array<double, ThreeComponent::nP> ThreeComponent::hi3Bounds{ { 1.e7, 0.750, 3.500, 7.500, 0.150, 0.250, 1.000, 0.49, 0.75, 1.500, 0. } };
const array<double, ThreeComponent::nP> ThreeComponent::lo7Bounds{ { 0.,   0.250, 0.750,  4.000, 0.010, 0.020, 0.150, 0.00, 0.00, 0.0, 0. } };
const array<double, ThreeComponent::nP> ThreeComponent::hi7Bounds{ { 1.e7, 0.750, 3.000, 20.000, 0.020, 0.050, 0.600, 0.95, 0.95, 0.5, 0. } };

//******************************************************************************
#pragma mark Utility functions
//******************************************************************************
template<typename Functor_t>
void write_results(const string outPrefix, double **paramsData,
				   double *residualData, NiftiImage &hdr)
{
	string outPath;
	hdr.setnt(1);
	hdr.setDatatype(NIFTI_TYPE_FLOAT32);
	for (int p = 0; p < Functor_t::nP; p++)
	{
		outPath = outPrefix + "_" + Functor_t::names[p] + ".nii.gz";
		cout << "Writing parameter file: " << outPath << endl;
		hdr.open(outPath, 'w');
		hdr.writeVolume(0, paramsData[p]);
		hdr.close();
	}
	
	if (residualData)
	{
		outPath = outPrefix + "_residual.nii.gz";
		cout << "Writing residual file: " << outPath << endl;
		hdr.open(outPath, 'w');
		hdr.writeVolume(0, residualData);
		hdr.close();
	}
}

#endif
