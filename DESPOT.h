/*
 *  DESPOT1.h
 *  MacRI
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright 2011 Tobias Wood. All rights reserved.
 *
 */
#ifndef __DESPOT__
#define __DESPOT__

#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/NonLinearOptimization>

using namespace Eigen;

//******************************************************************************
#pragma mark Convenience stuff
//******************************************************************************
#ifndef TRUE
	#define TRUE 1
#endif
#ifndef FALSE
	#define FALSE 0
#endif
#ifndef M_PI
	#define M_PI 3.14159265358979323846264338327950288
#endif
#ifndef M_SQRT2
	#define M_SQRT2     1.41421356237309504880168872420969808
#endif
extern int MATHS_DEBUG;

double radians(double degrees);
double degrees(double radians);

double clamp(double value, double low, double high);
//******************************************************************************
#pragma mark Parameter Files
//******************************************************************************
int fgetArray(FILE *in, char type, size_t n, void *array);

void linearLeastSquares(double *X, double *Y, int nD,
						double *slope, double *inter, double *res);
double classicDESPOT1(const ArrayXd &flipAngles, const ArrayXd &spgrVals,
				      double TR, double B1, double *M0, double *T1);
double classicDESPOT2(const ArrayXd &flipAngles, const ArrayXd &ssfpVals,
                      double TR, double T1, double B1, double *M0, double *T2);

// From Nonlinear Tests in Eigen 
template<typename _Scalar, int NX=Dynamic, int NY=Dynamic>
struct Functor
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
};

class DESPOT_Functor : public Functor<double>
{
	public:
		int _numValues;
		const ArrayXd &_flipAngles;
		const ArrayXXd &_signals;
		
		DESPOT_Functor(const ArrayXd &flipAngles, const ArrayXXd &signals) :
					   _flipAngles(flipAngles),
					   _signals(signals),
					   _numValues(flipAngles.size())
					   {}
		
		int values() const { return _numValues; }
		int operator()(const VectorXd &params, VectorXd &diffs) const 
		{ return 0; }
		
};


class SSFP_1c : public DESPOT_Functor
{
	public:
		const ArrayXd &_rfPhases;
		double _TR, _T1, _B1;
		SSFP_1c(const ArrayXd &rfPhases, const ArrayXd &flipAngles,
		        const ArrayXXd &signals, double TR, double T1, double B1) :
				DESPOT_Functor(flipAngles, signals),
				_rfPhases(rfPhases),
				_TR(TR),
				_T1(T1),
				_B1(B1)
				{ _numValues *= rfPhases.size(); }
	
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
					
					diffs[index] = sqrt(Mobs(0)*Mobs(0) + Mobs(1)*Mobs(1)) - 
								   _signals(p, i);
					index++;
				}
			}
			return 0;
		}
};

class SPGR_2c : public DESPOT_Functor
{
	public:
		double _TR, _M0, _B1;
		SPGR_2c(const ArrayXd &flipAngles, const ArrayXd &signals,
		        double TR, double M0, double B1) :
				DESPOT_Functor(flipAngles, signals),
				_TR(TR),
				_M0(M0),
				_B1(B1)
				{}
	
		int operator()(const VectorXd &params, VectorXd &diffs) const
		{
			double T1_a = params[0],
			       T1_b = params[1],
			       f_a  = params[4],
				   f_b = 1. - f_a,
			       tau_a = params[5],
			       tau_b = f_b * tau_a / f_a;
			
			Matrix2d A,
					 eye = Matrix2d::Identity(),
					 eyema;
			Vector2d M0, Mobs;
			A << -(_TR/T1_a + _TR/tau_a),               _TR/tau_b,
				               _TR/tau_a, -(_TR/T1_b + _TR/tau_b);
			M0 << _M0 * f_a, _M0 * f_b;
			MatrixExponential<Matrix2d> expA(A);
			expA.compute(A);
			eyema.noalias() = eye - A;
			for (int i = 0; i < _flipAngles.size(); i++)
			{
				double a = _flipAngles[i];
				Mobs = ((eye - A*cos(_B1 * a)).partialPivLu().solve(eyema * sin(_B1 * a))) * M0;
				diffs[i] = Mobs.sum() - _signals(i);
			}
			return 0;
		}
};

class SSFP_2c: public SSFP_1c
{
	public:
		double _M0;
		SSFP_2c(const ArrayXd &rfphases, const ArrayXd &flipAngles,
		        const ArrayXXd &signals, double TR, double M0, double T1, double B1) :
				SSFP_1c(rfphases, flipAngles, signals, TR, T1, B1),
				_M0(M0)
				{}		
		
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
			typedef Matrix<double, 6, 6> Matrix6d;
			typedef Matrix<double, 6, 1> Vector6d;
			Matrix6d A,
			         R_rf = Matrix6d::Zero(),
					 eye = Matrix6d::Identity(),
					 eyema;
			Vector6d M0, Mobs;
			R_rf(0, 0) = R_rf(1, 1) = 1.;
			M0 << 0., 0., 0., 0., _M0 * f_a, _M0 * f_b;
			
			double eT2_a = -_TR * (1./T2_a + 1./tau_a),
				   eT2_b = -_TR * (1./T2_b + 1./tau_b),
				   eT1_a = -_TR * (1./T1_a + 1./tau_a),
				   eT1_b = -_TR * (1./T1_b + 1./tau_b),
				   k_a   = _TR / tau_a,
				   k_b   = _TR / tau_b;
			
			int index = 0;
			for (int p = 0; p < _rfPhases.size(); p++)
			{
				double phase = _rfPhases[p] + (B0 * _TR * 2. * M_PI);
				A << eT2_a,    k_b, phase,    0.,    0.,    0.,
					   k_a,  eT2_b,    0., phase,    0.,    0.,
					-phase,     0., eT2_a,   k_b,    0.,    0.,
						0., -phase,   k_a, eT2_b,    0.,    0.,
						0.,     0.,    0.,    0., eT1_a,   k_b,
						0.,     0.,    0.,    0.,   k_a, eT1_b;
				MatrixExponential<Matrix6d> expM(A);
				expM.compute(A);
				eyema.noalias() = eye - A;
				for (int i = 0; i < _flipAngles.size(); i++)
				{
					double a = _flipAngles[i];
					double ca = cos(_B1 * a), sa = sin(_B1 * a);
					R_rf(2, 2) = R_rf(3, 3) = R_rf(4, 4) = R_rf(5, 5) = ca;
					R_rf(2, 4) = R_rf(3, 5) = sa;
					R_rf(4, 2) = R_rf(5, 3) = -sa;
					
					Mobs = (eye - (A * R_rf)).partialPivLu().solve(eyema) * M0;
					diffs[index] = sqrt(pow(Mobs[0] + Mobs[1], 2.) +
					                    pow(Mobs[2] + Mobs[3], 2.)) - _signals(p, i);
					index++;
				}
			}
			return 0;
		}
};

#endif