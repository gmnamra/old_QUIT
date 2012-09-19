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
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "nifti3_io.h"

using namespace Eigen;

//******************************************************************************
#pragma mark Functors
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
		
		virtual ~Functor() {};
		
		long inputs() const { return m_inputs; }
		long values() const { return m_values; }
		
		virtual int operator()(const VectorXd &params, VectorXd &diffs) const = 0;
};
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
#pragma mark OneComponent
//******************************************************************************
class OneComponent : public Functor<double>
{
	public:
		const VectorXd &_spgrAngles, &_ssfpAngles, &_rfPhases, &_spgrSignals;
		const std::vector<VectorXd> &_ssfpSignals;
		double _spgrTR, _ssfpTR, _B0, _B1;
		
		static const int nP = 3;
		static const char *names[];
		static const double lo3Bounds[], hi3Bounds[], lo7Bounds[], hi7Bounds[];
		
		static bool f_constraint(const VectorXd &params)
		{
			return true;
		}
		
		OneComponent(const VectorXd &spgrAngles, const VectorXd &spgrSignals,
		             const VectorXd &ssfpAngles, const VectorXd &rfPhases,
					 const std::vector<VectorXd> &ssfpSignals,
					 const double spgrTR, const double ssfpTR, 
					 const double B0, const double B1) :
					 Functor<double>(OneComponent::nP, spgrAngles.size() + ssfpAngles.size() * rfPhases.size()),
				     _spgrAngles(spgrAngles), _spgrSignals(spgrSignals),
				     _ssfpAngles(ssfpAngles), _rfPhases(rfPhases),
					 _ssfpSignals(ssfpSignals),
					 _spgrTR(spgrTR), _ssfpTR(ssfpTR), _B0(B0), _B1(B1)
				     {
						
					 }
	
		int operator()(const VectorXd &params, VectorXd &diffs) const
		{
			eigen_assert(diffs.size() == values());
			
			int index = 0;
			//std::cout << "****************************************************" << std::endl;
			//std::cout << "SPGR First" << std::endl;
			VectorXd theory = One_SPGR(_spgrAngles, _spgrTR, params[0], params[1], _B1);
			//std::cout << "SPGR Theory: " << theory.transpose() << " Sig: " << _spgrSignals.transpose() << std::endl;
			diffs.head(_spgrSignals.size()) = theory - _spgrSignals;
			index += _spgrSignals.size();
			
			//std::cout << "****************************************************" << std::endl;
			//std::cout << "Now SSFP" << std::endl;
			for (int p = 0; p < _rfPhases.size(); p++)
			{
				theory = One_SSFP(_ssfpAngles, _rfPhases(p), _ssfpTR, params[0], params[1], params[2], _B0, _B1);
				diffs.segment(index, _ssfpAngles.size()) = theory - _ssfpSignals[p];
				//std::cout << "SSFP Theory: " << theory.transpose() << " Sig: " << _ssfpSignals[p].transpose() << std::endl;
				index += _ssfpAngles.size();
			}
			//std::cout << "Diffs: " << diffs.transpose() << std::endl;
			return 0;
		}
};
const char *OneComponent::names[] = { "M0", "T1", "T2" };
const double OneComponent::lo3Bounds[] = { 0., 0.1, 0.01 };
const double OneComponent::hi3Bounds[] = { 1.e7, 3.0, 1.0 };
const double OneComponent::lo7Bounds[] = { 0., 0., 0.005 };
const double OneComponent::hi7Bounds[] = { 1.e7, 3., 0.1 };

//******************************************************************************
#pragma mark OneComponentSSFP
//******************************************************************************

class OneComponentSSFP : public Functor<double>
{
	public:
		const VectorXd &_flipAngles, &_rfPhases;
		const std::vector<VectorXd> &_signals;
		const double _TR, _T1, _M0, _B0, _B1;
		
		static const int nP = 1;
		static const char *names[];
		
		OneComponentSSFP(const VectorXd &flipAngles, const VectorXd &rfPhases,
		                 const std::vector<VectorXd> &signals,
				         const double TR, const double M0, const double T1,
						 const double B0, const double B1) :
				         Functor<double>(OneComponentSSFP::nP, flipAngles.size() * rfPhases.size()),
				         _flipAngles(flipAngles),
				         _rfPhases(rfPhases),
				         _signals(signals),
				         _TR(TR), _T1(T1), _M0(M0), _B0(B0), _B1(B1)
				         {}
			
		int operator()(const VectorXd &params, VectorXd &diffs) const
		{
			eigen_assert(diffs.size() == values());
			int index = 0;
			for (int p = 0; p < _rfPhases.size(); p++)
			{
				VectorXd temp = One_SSFP(_flipAngles, _rfPhases(p), _TR, _M0, _T1, params[0], _B0, _B1);
				diffs.segment(index, _flipAngles.size()) = temp - _signals[p];
				index += _flipAngles.size();
			}
			return 0;
		}
};
const char *OneComponentSSFP::names[] = { "T2" };

//******************************************************************************
#pragma mark Two Component
//******************************************************************************

typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;
class TwoComponent : public Functor<double>
{
	public:
		const VectorXd &_spgrAngles, &_ssfpAngles, &_rfPhases, &_spgrSignals;
		const std::vector<VectorXd> &_ssfpSignals;
		const double _spgrTR, _ssfpTR, _B0, _B1;
		
		static const int nP = 7;
		static const char *names[];
		static const double lo3Bounds[], hi3Bounds[], lo7Bounds[], hi7Bounds[];
		
		static bool f_constraint(const VectorXd &params)
		{
			return true;
		}
		
		TwoComponent(const VectorXd &spgrAngles, const VectorXd &spgrSignals,
		             const VectorXd &ssfpAngles, const VectorXd &rfPhases,
					 const std::vector<VectorXd> &ssfpSignals,
					 const double spgrTR, const double ssfpTR, 
					 const double B0, const double B1) :
					 Functor<double>(TwoComponent::nP, spgrAngles.size() + ssfpAngles.size() * rfPhases.size()),
				     _spgrAngles(spgrAngles), _spgrSignals(spgrSignals),
				     _ssfpAngles(ssfpAngles), _rfPhases(rfPhases),
					 _ssfpSignals(ssfpSignals),
					 _spgrTR(spgrTR), _ssfpTR(ssfpTR), _B0(B0), _B1(B1)
				     {
						
					 }
	
		int operator()(const VectorXd &params, VectorXd &diffs) const
		{
			double PD   = params[0],
			       T1_a = params[1],
			       T1_b = params[2],
				   T2_a = params[3],
				   T2_b = params[4],
			       f_a  = params[5],
				   f_b = 1. - f_a,
			       tau_a = params[6],
			       tau_b = f_b * tau_a / f_a,
				   k_ab = 1. / tau_a,
				   k_ba = 1. / tau_b;
			// Only have 1 component, so no exchange
			if ((f_a == 0.) || (f_b == 0.))
			{
				k_ab = 0.;
				k_ba = 0.;
			}
			
			eigen_assert(diffs.size() == values());
			
			int index = 0;
			//std::cout << "****************************************************" << std::endl;
			//std::cout << "SPGR First" << std::endl;
			{
				Matrix2d A;
				const Matrix2d eye2 = Matrix2d::Identity();
				Vector2d M0, Mobs;
				VectorXd signals(_spgrSignals.size());
				
				M0 << PD * f_a, PD * f_b;
				A << -(_spgrTR/T1_a + _spgrTR*k_ab),                  _spgrTR*k_ba,
									   _spgrTR*k_ab, -(_spgrTR/T1_b + _spgrTR*k_ba);
				//std::cout << A << std::endl;
				MatrixExponential<Matrix2d> expA(A);
				expA.compute(A);
				const Matrix2d eyema = eye2 - A;
				for (int i = 0; i < _spgrAngles.size(); i++)
				{
					double a = _spgrAngles[i];
					Mobs = (eye2 - A*cos(_B1 * a)).partialPivLu().solve(eyema * sin(_B1 * a)) * M0;
					signals[i] = Mobs.sum();
				}
				diffs.head(_spgrSignals.size()) = signals - _spgrSignals;
				index += _spgrSignals.size();
			}
			
			//std::cout << "****************************************************" << std::endl;
			//std::cout << "Now SSFP" << std::endl;
			{
				Matrix6d A, expA, R_rf, eye_mAR;
				const Matrix6d eye6 = Matrix6d::Identity();
				Vector6d M0, Mobs;
				PartialPivLU<Matrix6d> solver;
				A.setZero(); R_rf.setZero();
				R_rf(0, 0) = R_rf(1, 1) = 1.;
				M0 << 0., 0., 0., 0., PD * f_a, PD * f_b;
				
				for (int p = 0; p < _rfPhases.size(); p++)
				{
					VectorXd signals(_ssfpAngles.size());
					double phase = _rfPhases[p] + (_B0 * _ssfpTR * 2. * M_PI);
					// Can get away with this because the block structure of the
					// matrix ensures that the zero blocks are always zero after
					// the matrix exponential.
					A(0, 0) = A(2, 2) = -_ssfpTR * (1./T2_a + k_ab);
					A(1, 1) = A(3, 3) = -_ssfpTR * (1./T2_b + k_ba);
					A(0, 1) = A(2, 3) = A(4, 5) = _ssfpTR * k_ba;
					A(0, 2) = A(1, 3) = phase;
					A(1, 0) = A(3, 2) = A(5, 4) = _ssfpTR * k_ab;
					A(2, 0) = A(3, 1) = -phase; 
					A(4, 4) = -_ssfpTR * (1./T1_a + k_ab);
					A(5, 5) = -_ssfpTR * (1./T1_b + k_ba);
					A(0, 3) = A(1, 2) = A(2, 1) = A(3, 0) = 0.;
					MatrixExponential<Matrix6d> exp(A);
					exp.compute(expA);
					const Matrix6d eyema = eye6 - expA;
					for (int i = 0; i < _ssfpAngles.size(); i++)
					{
						double a = _ssfpAngles[i];
						double ca = cos(_B1 * a), sa = sin(_B1 * a);
						R_rf(2, 2) = R_rf(3, 3) = R_rf(4, 4) = R_rf(5, 5) = ca;
						R_rf(2, 4) = R_rf(3, 5) = sa;
						R_rf(4, 2) = R_rf(5, 3) = -sa;
						eye_mAR.noalias() = eye6 - (expA * R_rf);
						solver.compute(eye_mAR);
						Mobs.noalias() = solver.solve(eyema) * M0;
						signals[i] = sqrt(pow(Mobs[0] + Mobs[1], 2.) +
										  pow(Mobs[2] + Mobs[3], 2.));
					}
					diffs.segment(index, _ssfpAngles.size()) = signals - _ssfpSignals[p];
					index += _ssfpAngles.size();
				}
			}
			return 0;
		}
};
const char *TwoComponent::names[] = { "M0", "T1_a", "T1_b", "T2_a", "T2_b", "f_a", "tau_a" };
const double TwoComponent::lo3Bounds[] = { 0., 0.1, 0.8, 0.001, 0.01, 0.0, 0.05, 0. };
const double TwoComponent::hi3Bounds[] = { 1.e7, 1.0, 3.0, 0.050, 0.25, 1.0, 2.00, 0. };
const double TwoComponent::lo7Bounds[] = { 0., 0.1, 0.8, 0.001, 0.01, 0.0, 0.05, 0. };
const double TwoComponent::hi7Bounds[] = { 1.e7, 1.0, 3.0, 0.050, 0.25, 1.0, 2.00, 0. };

//******************************************************************************
#pragma mark Three Component
//******************************************************************************

typedef Matrix<double, 9, 9> Matrix9d;
typedef Matrix<double, 9, 1> Vector9d;
class ThreeComponent : public Functor<double>
{
	private:
		const double _spgrTR, _ssfpTR, _B0, _B1;
		const VectorXd &_spgrAngles, &_ssfpAngles, &_rfPhases, &_spgrSignals;
		const std::vector<VectorXd> &_ssfpSignals;
		
	public:
		static const int nP = 10;
		static const char *names[];
		static const double lo3Bounds[], hi3Bounds[], lo7Bounds[], hi7Bounds[];
		
		static bool f_constraint(const VectorXd &params)
		{
			if ((params[7] + params[8]) > 0.95)
				return false;
			else
				return true;
		}
		
		
		ThreeComponent(const VectorXd &spgrAngles, const VectorXd &spgrSignals,
		               const VectorXd &ssfpAngles, const VectorXd &rfPhases,
					   const std::vector<VectorXd> &ssfpSignals,
					   const double spgrTR, const double ssfpTR, 
					   const double B0, const double B1) :
					   Functor<double>(ThreeComponent::nP, spgrAngles.size() + ssfpAngles.size() * rfPhases.size()),
				       _spgrAngles(spgrAngles), _spgrSignals(spgrSignals),
				       _ssfpAngles(ssfpAngles), _rfPhases(rfPhases),
					   _ssfpSignals(ssfpSignals),
					   _spgrTR(spgrTR), _ssfpTR(ssfpTR), _B0(B0), _B1(B1)
				      {
					  	eigen_assert(_rfPhases.size() == _ssfpSignals.size());
					  }
	
		int operator()(const VectorXd &params, VectorXd &diffs) const
		{
			double PD   = params[0],
			       T1_a = params[1],
			       T1_b = params[2],
				   T1_c = params[3],
				   T2_a = params[4],
				   T2_b = params[5],
				   T2_c = params[6],
			       f_a  = params[7],
				   f_c  = params[8],
				   f_b  = 1. - f_a - f_c,
			       tau_a = params[9],
			       tau_b = f_b * tau_a / f_a,
				   k_ab = 1. / tau_a,
				   k_ba = 1. / tau_b;
			// Only have 1 component, so no exchange
			if ((f_a == 0.) || (f_b == 0.))
			{
				k_ab = 0.;
				k_ba = 0.;
			}
			
			eigen_assert(diffs.size() == values());
			
			int index = 0;
			//std::cout << "****************************************************" << std::endl;
			//std::cout << "SPGR First" << std::endl;
			{
				Matrix3d A;
				const Matrix3d eye3 = Matrix3d::Identity();
				Vector3d M0, Mobs;
				VectorXd signals(_spgrSignals.size());
				M0 << PD * f_a, PD * f_b, PD * f_c;
				A << -(1./T1_a + k_ab),               k_ba, 0,
								  k_ab,  -(1./T1_b + k_ba), 0,
									 0,                  0, -1./T1_c;
				//std::cout << A << std::endl;
				MatrixExponential<Matrix3d> expA(_spgrTR * A);
				expA.compute(A);
				const Matrix3d eyema = eye3 - A;
				for (int i = 0; i < _spgrAngles.size(); i++)
				{
					double a = _spgrAngles[i];
					Mobs = (eye3 - A*cos(_B1 * a)).partialPivLu().solve(eyema * sin(_B1 * a)) * M0;
					signals[i] = Mobs.sum();
				}
				diffs.head(_spgrSignals.size()) = signals - _spgrSignals;
				index += _spgrSignals.size();
			}
			
			//std::cout << "****************************************************" << std::endl;
			//std::cout << "Now SSFP" << std::endl;
			{
				Matrix9d A, expA, R_rf, eye_mAR;
				const Matrix9d eye9 = Matrix9d::Identity();
				Vector9d M0, Mobs;
				PartialPivLU<Matrix9d> solver;
				R_rf.setZero();
				// Set up the 'A' matrix. It's quite complex.
				A(0, 0) = A(3, 3) = -_ssfpTR * (1./T2_a + k_ab);
				A(1, 1) = A(4, 4) = -_ssfpTR * (1./T2_b + k_ba);
				A(2, 2) = A(5, 5) = -_ssfpTR * (1./T2_c);
				A(0, 1) = A(3, 4) = A(6, 7) = _ssfpTR * k_ba;
				A(1, 0) = A(4, 3) = A(7, 6) = _ssfpTR * k_ab;
				A(6, 6) = -_ssfpTR * (1./T1_a + k_ab);
				A(7, 7) = -_ssfpTR * (1./T1_b + k_ba);
				A(8, 8) = -_ssfpTR * (1./T1_c);
				R_rf(0, 0) = R_rf(1, 1) = R_rf(2, 2) = 1.;
				A(0, 2) = A(0, 4) = A(0, 5) =
				A(1, 2) = A(1, 3) = A(1, 5) =
				A(2, 0) = A(2, 1) = A(2, 3) = A(2, 4) =
				A(3, 1) = A(3, 2) = A(3, 5) = 
				A(4, 0) = A(4, 2) = A(4, 5) =
				A(5, 0) = A(5, 1) = A(5, 3) = A(5, 4) =
				A(6, 8) = A(7, 8) = A(8, 6) = A(8, 7) = 0.;
				A.topRightCorner(6, 3).setZero();
				A.bottomLeftCorner(3, 6).setZero();
				M0 << 0., 0., 0., 0., 0., 0., PD * f_a, PD * f_b, PD * f_c;
				
				for (int p = 0; p < _rfPhases.size(); p++)
				{
					VectorXd signals(_ssfpAngles.size());
					double phase = _rfPhases[p] + (_B0 * _ssfpTR * 2. * M_PI);
					A(0, 3) = A(1, 4) = A(2, 5) = phase;
					A(3, 0) = A(4, 1) = A(5, 2) = -phase;
					MatrixExponential<Matrix9d> exp(A);
					exp.compute(expA);
					const Matrix9d eyema = eye9 - expA;
					for (int i = 0; i < _ssfpAngles.size(); i++)
					{
						double a = _ssfpAngles[i];
						double ca = cos(_B1 * a), sa = sin(_B1 * a);
						R_rf(3, 3) = R_rf(4, 4) = R_rf(5, 5) =
						R_rf(6, 6) = R_rf(7, 7) = R_rf(8, 8) =  ca;
						R_rf(3, 6) = R_rf(4, 7) = R_rf(5, 8) =  sa;
						R_rf(6, 3) = R_rf(7, 4) = R_rf(8, 5) = -sa;
						eye_mAR.noalias() = eye9 - (expA * R_rf);
						solver.compute(eye_mAR);
						Mobs.noalias() = solver.solve(eyema) * M0;
						signals[i] = sqrt(pow(Mobs[0] + Mobs[1] + Mobs[2], 2.) +
										  pow(Mobs[3] + Mobs[4] + Mobs[5], 2.));
					}
					diffs.segment(index, _ssfpSignals[p].size()) = signals - _ssfpSignals[p];
					index += _ssfpSignals[p].size();
				}
			}
			return 0;
		}
};
const char *ThreeComponent::names[] = { "M0", "3c_T1_a", "3c_T1_b", "3c_T1_c", "3c_T2_a", "3c_T2_b", "3c_T2_c", "3c_f_a", "3c_f_c", "3c_tau_a" };
const double ThreeComponent::lo3Bounds[] = { 0., 0.250, 0.250, 1.500, 0.000, 0.000, 0.150, 0.00, 0.00, 0.025 };
const double ThreeComponent::hi3Bounds[] = { 1.e7, 0.750, 3.500, 7.500, 0.150, 0.250, 1.000, 0.49, 0.75, 1.500 };
const double ThreeComponent::lo7Bounds[] = { 0., 0.500, 1.50, 0.0001, 0.010, 0.0, 0.0, 0., 0., 0. };
const double ThreeComponent::hi7Bounds[] = { 1.e7, 1.000, 3.00, 0.0500, 0.500, 1.0, 1.0, 0., 0., 0. };

//******************************************************************************
#pragma mark Utility functions
//******************************************************************************
template<typename Functor_t>
void write_results(const std::string outPrefix, double **paramsData,
				   double *residualData, NiftiImage &hdr)
{
	std::string outPath;
	hdr.setnt(1);
	hdr.setDatatype(NIFTI_TYPE_FLOAT32);
	for (int p = 0; p < Functor_t::nP; p++)
	{
		outPath = outPrefix + "_" + Functor_t::names[p] + ".nii.gz";
		std::cout << "Writing parameter file: " << outPath << std::endl;
		hdr.open(outPath, 'w');
		hdr.writeVolume(0, paramsData[p]);
		hdr.close();
	}
	
	if (residualData)
	{
		outPath = outPrefix + "_residual.nii.gz";
		std::cout << "Writing residual file: " << outPath << std::endl;
		hdr.open(outPath, 'w');
		hdr.writeVolume(0, residualData);
		hdr.close();
	}
}

#endif
