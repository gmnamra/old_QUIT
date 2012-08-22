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

using namespace Eigen;
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

  const int m_inputs, m_values;

  Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
  Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

  int inputs() const { return m_inputs; }
  int values() const { return m_values; }
  
  virtual int operator()(const VectorXd &params, VectorXd &diffs) { return 0; }
};

class SSFP_1c : public Functor<double>
{
	public:
		const VectorXd &_flipAngles, &_rfPhases;
		const std::vector<VectorXd> &_signals;
		double _TR, _T1, _B1;
		SSFP_1c(const VectorXd &flipAngles, const VectorXd &rfPhases,
		        const std::vector<VectorXd> &signals,
				double TR, double T1, double B1) :
				Functor(3, flipAngles.size() * rfPhases.size()),
				_flipAngles(flipAngles),
				_rfPhases(rfPhases),
				_signals(signals),
				_TR(TR),
				_T1(T1),
				_B1(B1)
				{
				}
	
		int operator()(const VectorXd &params, VectorXd &diffs) const
		{
			Matrix3d A,
					 R_rf = Matrix3d::Zero(),
					 eye = Matrix3d::Identity(),
					 eyema;
			Vector3d M0, Mobs;
			R_rf(0, 0) = 1.;
			M0 << 0., 0., params[0];
			double T2 = params[1],
				   B0 = params[2];
			int index = 0;
			for (int p = 0; p < _rfPhases.size(); p++)
			{
				eigen_assert(diffs.size() == values());
				
				double phase = _rfPhases(p) + (B0 * _TR * 2. * M_PI);
				A << -_TR / T2,     phase,         0.,
					    -phase, -_TR / T2,         0.,
						    0.,        0., -_TR / _T1;
				MatrixExponential<Matrix3d> expA(A);
				expA.compute(A);
				eyema.noalias() = eye - A;
				for (int i = 0; i < _flipAngles.size(); i++)
				{
					double a = _flipAngles[i];
					double ca = cos(_B1 * a), sa = sin(_B1 * a);
					R_rf(1, 1) = R_rf(2, 2) = ca;
					R_rf(1, 2) = sa; R_rf(2, 1) = -sa;
					Mobs = (eye - (A * R_rf)).partialPivLu().solve(eyema) * M0;				
					diffs[index] = Mobs.head(2).norm() - _signals[p][i];
					index++;
				}
			}
			return 0;
		}
};

typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;

class TwoComponent : public Functor<double>
{
	public:
		const VectorXd &_spgrAngles, &_ssfpAngles, &_rfPhases, &_spgrSignals;
		const std::vector<VectorXd> &_ssfpSignals;
		const double _spgrTR, _ssfpTR, _M0, _B0, _B1;
		
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
					 const double M0, const double B0, const double B1) :
					 Functor(TwoComponent::nP, spgrAngles.size() + ssfpAngles.size() * rfPhases.size()),
				     _spgrAngles(spgrAngles), _spgrSignals(spgrSignals),
				     _ssfpAngles(ssfpAngles), _rfPhases(rfPhases),
					 _ssfpSignals(ssfpSignals),
					 _spgrTR(spgrTR), _ssfpTR(ssfpTR), _M0(M0), _B0(B0), _B1(B1)
				     {
						
					 }
	
		int operator()(const VectorXd &params, VectorXd &diffs) const
		{
			double T1_a = params[0],
			       T1_b = params[1],
				   T2_a = params[2],
				   T2_b = params[3],
			       f_a  = params[4],
				   f_b = 1. - f_a,
			       tau_a = params[5],
			       tau_b = f_b * tau_a / f_a,
				   B0 = params[6],
				   k_ab = 1. / tau_a,
				   k_ba = 1. / tau_b;
			
			if (std::isfinite(_B0))
				B0 = _B0;
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
				M0 << _M0 * f_a, _M0 * f_b;
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
					diffs[index] = Mobs.sum() - _spgrSignals[i];
					index++;
				}
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
				M0 << 0., 0., 0., 0., _M0 * f_a, _M0 * f_b;
				
				for (int p = 0; p < _rfPhases.size(); p++)
				{
					double phase = _rfPhases[p] + (B0 * _ssfpTR * 2. * M_PI);
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
						diffs[index] = sqrt(pow(Mobs[0] + Mobs[1], 2.) +
											pow(Mobs[2] + Mobs[3], 2.)) - _ssfpSignals[p][i];
						index++;
					}
				}
			}
			return 0;
		}
};

const char *TwoComponent::names[] = { "T1_a", "T1_b", "T2_a", "T2_b", "f_a", "tau_a", "B0" };
const double TwoComponent::lo3Bounds[] = { 0.1, 0.8, 0.001, 0.01, 0.0, 0.05, 0. };
const double TwoComponent::hi3Bounds[] = { 1.0, 3.0, 0.050, 0.25, 1.0, 2.00, 0. };
const double TwoComponent::lo7Bounds[] = { 0.1, 0.8, 0.001, 0.01, 0.0, 0.05, 0. };
const double TwoComponent::hi7Bounds[] = { 1.0, 3.0, 0.050, 0.25, 1.0, 2.00, 0. };

typedef Matrix<double, 9, 9> Matrix9d;
typedef Matrix<double, 9, 1> Vector9d;

class ThreeComponent : public Functor<double>
{
	public:
		const VectorXd &_spgrAngles, &_ssfpAngles, &_rfPhases, &_spgrSignals;
		const std::vector<VectorXd> &_ssfpSignals;
		const double _spgrTR, _ssfpTR, _M0, _B0, _B1;
		
		static const int nP = 10;
		static const char *names[];
		static const double lo3Bounds[], hi3Bounds[], lo7Bounds[], hi7Bounds[];
		
		static bool f_constraint(const VectorXd &params)
		{
			if ((params[6] + params[7]) > 0.95)
				return false;
			else
				return true;
		}
		
		
		ThreeComponent(const VectorXd &spgrAngles, const VectorXd &spgrSignals,
		               const VectorXd &ssfpAngles, const VectorXd &rfPhases,
					   const std::vector<VectorXd> &ssfpSignals,
					   const double spgrTR, const double ssfpTR, 
					   const double M0, const double B0, const double B1) :
					   Functor(ThreeComponent::nP, spgrAngles.size() + ssfpAngles.size() * rfPhases.size()),
				       _spgrAngles(spgrAngles), _spgrSignals(spgrSignals),
				       _ssfpAngles(ssfpAngles), _rfPhases(rfPhases),
					   _ssfpSignals(ssfpSignals),
					   _spgrTR(spgrTR), _ssfpTR(ssfpTR), _M0(M0), _B0(B0), _B1(B1)
				      {
						
					  }
	
		int operator()(const VectorXd &params, VectorXd &diffs) const
		{
			double T1_a = params[0],
			       T1_b = params[1],
				   T1_c = params[2],
				   T2_a = params[3],
				   T2_b = params[4],
				   T2_c = params[5],
			       f_a  = params[6],
				   f_c  = params[7],
				   f_b  = 1. - f_a - f_c,
			       tau_a = params[8],
			       tau_b = f_b * tau_a / f_a,
				   B0 = params[9],
				   k_ab = 1. / tau_a,
				   k_ba = 1. / tau_b;
			
			if (std::isfinite(_B0))
				B0 = _B0;
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
				M0 << _M0 * f_a, _M0 * f_b, _M0 * f_c;
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
					diffs[index] = Mobs.sum() - _spgrSignals[i];
					index++;
				}
			}
			
			//std::cout << "****************************************************" << std::endl;
			//std::cout << "Now SSFP" << std::endl;
			{
				Matrix9d A, expA, R_rf, eye_mAR;
				const Matrix9d eye9 = Matrix9d::Identity();
				Vector9d M0, Mobs;
				PartialPivLU<Matrix9d> solver;
				A.setZero(); R_rf.setZero();
				R_rf(0, 0) = R_rf(1, 1) = R_rf(2, 2) = 1.;
				M0 << 0., 0., 0., 0., 0., 0., _M0 * f_a, _M0 * f_b, _M0 * f_c;
				
				for (int p = 0; p < _rfPhases.size(); p++)
				{
					double phase = _rfPhases[p] + (B0 * _ssfpTR * 2. * M_PI);
					// Can get away with this because the block structure of the
					// matrix ensures that the zero blocks are always zero after
					// the matrix exponential.
					A(0, 0) = A(3, 3) = -_ssfpTR * (1./T2_a + k_ab);
					A(1, 1) = A(4, 4) = -_ssfpTR * (1./T2_b + k_ba);
					A(2, 2) = A(5, 5) = -_ssfpTR * (1./T2_c);
					A(0, 1) = A(3, 4) = A(6, 7) = _ssfpTR * k_ba;
					A(0, 3) = A(1, 4) = A(2, 5) = phase;
					A(1, 0) = A(4, 3) = A(7, 6) = _ssfpTR * k_ab;
					A(3, 0) = A(4, 1) = A(5, 2) = -phase; 
					A(6, 6) = -_ssfpTR * (1./T1_a + k_ab);
					A(7, 7) = -_ssfpTR * (1./T1_b + k_ba);
					A(8, 8) = -_ssfpTR * (1./T1_c);
					A(0, 2) = A(0, 4) = A(0, 5) =
					A(1, 2) = A(1, 3) = A(1, 5) =
					A(2, 0) = A(2, 1) = A(2, 3) = A(2, 4) =
					A(3, 1) = A(3, 2) = A(3, 5) = 
					A(4, 0) = A(4, 2) = A(4, 5) =
					A(5, 0) = A(5, 1) = A(5, 3) = A(5, 4) =
					A(6, 8) = A(7, 8) = A(8, 6) = A(8, 7) = 0.;
					A.topRightCorner(6, 3).setZero();
					A.bottomLeftCorner(3, 6).setZero();
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
						diffs[index] = sqrt(pow(Mobs[0] + Mobs[1] + Mobs[2], 2.) +
											pow(Mobs[3] + Mobs[4] + Mobs[5], 2.)) - _ssfpSignals[p][i];
						index++;
					}
				}
			}
			return 0;
		}
};

const char *ThreeComponent::names[] = { "T1_a", "T1_b", "T1_c", "T2_a", "T2_b", "T2_c", "f_a", "f_c", "tau_a", "B0" };
const double ThreeComponent::lo3Bounds[] = { 0.250, 0.250, 1.500, 0.000, 0.000, 0.150, 0.00, 0.00, 0.025, 0. };
const double ThreeComponent::hi3Bounds[] = { 0.750, 3.500, 7.500, 0.150, 0.250, 1.000, 0.49, 0.75, 1.500, 0. };
const double ThreeComponent::lo7Bounds[] = { 0.500, 1.50, 0.0001, 0.010, 0.0, 0.0, 0., 0., 0., 0. };
const double ThreeComponent::hi7Bounds[] = { 1.000, 3.00, 0.0500, 0.500, 1.0, 1.0, 0., 0., 0., 0. };

#endif
