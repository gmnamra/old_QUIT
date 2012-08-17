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
				{  }
	
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
		double _spgrTR, _ssfpTR, _M0, _B1;
		TwoComponent(const VectorXd &spgrAngles, const VectorXd &spgrSignals,
		             const VectorXd &ssfpAngles, const VectorXd &rfPhases,
					 const std::vector<VectorXd> &ssfpSignals,
					 const double spgrTR, const double ssfpTR, 
					 const double M0, const double B1) :
					 Functor(7, spgrAngles.size() + ssfpAngles.size() * rfPhases.size()),
				     _spgrAngles(spgrAngles), _spgrSignals(spgrSignals),
				     _ssfpAngles(ssfpAngles), _rfPhases(rfPhases),
					 _ssfpSignals(ssfpSignals),
					 _spgrTR(spgrTR), _ssfpTR(ssfpTR), _M0(M0), _B1(B1)
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
				   B0 = params[6];
			
			eigen_assert(diffs.size() == values());
			
			int index = 0;
			//****************************************************
			// SPGR First
			{
				Matrix2d A;
				const Matrix2d eye2 = Matrix2d::Identity();
				Vector2d M0, Mobs;
				M0 << _M0 * f_a, _M0 * f_b;
				A << -(_spgrTR/T1_a + _spgrTR/tau_a),                  _spgrTR/tau_b,
									   _spgrTR/tau_a, -(_spgrTR/T1_b + _spgrTR/tau_b);
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
				
				double eT2_a = -_ssfpTR * (1./T2_a + 1./tau_a),
					   eT2_b = -_ssfpTR * (1./T2_b + 1./tau_b),
					   eT1_a = -_ssfpTR * (1./T1_a + 1./tau_a),
					   eT1_b = -_ssfpTR * (1./T1_b + 1./tau_b),
					   k_a   = _ssfpTR / tau_a,
					   k_b   = _ssfpTR / tau_b;
				
				for (int p = 0; p < _rfPhases.size(); p++)
				{
					double phase = _rfPhases[p] + (B0 * _ssfpTR * 2. * M_PI);
					// Can get away with this because the block structure of the
					// matrix ensures that the zero blocks are always zero after
					// the matrix exponential.
					A(0, 0) = A(2, 2) = eT2_a;
					A(1, 1) = A(3, 3) = eT2_b;
					A(0, 1) = A(2, 3) = A(4, 5) = k_b;
					A(0, 2) = A(1, 3) = phase;
					A(1, 0) = A(3, 2) = A(5, 4) = k_a;
					A(2, 0) = A(3, 1) = -phase;  
					A(4, 4) = eT1_a;
					A(5, 5) = eT1_b;
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

#endif
