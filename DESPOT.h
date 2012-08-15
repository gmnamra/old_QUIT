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

#include "mathsMatrix.h"
#include "mathsOptimisers.h"

#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/NonLinearOptimization>

using namespace Eigen;
int tests();
double SPGR(double flipAngle, double *p, double *c);
void SPGR_Jacobian(double *angles, int nD, double *p, double *c, double *result);
double IRSPGR(double TI, double *p, double *c);
void IRSPGR_Jacobian(double *data, int nD, double *par, double *c, double *result);

double classicDESPOT2(const Map<ArrayXd> &flipAngles, ArrayXd &ssfpVals,
                      double TR, double T1, double B1, double *M0, double *T2);

typedef struct
{
	double TR, M0, B1;
} SPGR_constants;

typedef struct
{
	double TR, M0, T1, B0, B1, rfPhase;
} SSFP_constants;

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

  // you should define that in the subclass :
//  void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const;
};

class a1cSSFP_functor : public Functor<double>
{
	public:
	const int m_values;
	const ArrayXd &_rf_phase, &_angles;
	const ArrayXXd &_signal;
	double _TR, _T1, _B1;
	a1cSSFP_functor(const ArrayXd &phases, const ArrayXd &angles, ArrayXXd &signal,
					double TR, double T1, double B1) :
					_rf_phase(phases),
					_angles(angles),
					_signal(signal),
					_TR(TR),
					_T1(T1),
					_B1(B1),
					m_values(angles.size() * phases.size())
					{}
	
	int values() const {return m_values;}
	int operator()(const VectorXd &params, VectorXd &diffs) const
	{
		Matrix3d A,
		         R_rf = Matrix3d::Zero(),
		         eye = Matrix3d::Identity(),
				 eyema;
		R_rf(0, 0) = 1.;
		Vector3d M0, Mobs;
		M0 << 0., 0., params[0];
		double T2 = params[1],
			   B0 = params[2];
		
		int index = 0;
		for (int p = 0; p < _rf_phase.size(); p++)
		{
			double phase = _rf_phase(p) + (B0 * _TR * 2. * M_PI);
			A << -_TR / T2,    phase,       0.,
				   -phase, -_TR / T2,       0.,
					   0.,       0., -_TR / _T1;
			MatrixExponential<Matrix3d> expA(A);
			expA.compute(A);
			eyema.noalias() = eye - A;
			for (int i = 0; i < _angles.size(); i++)
			{
				double a = _angles[i];
				double ca = cos(_B1 * a), sa = sin(_B1 * a);
				R_rf(1, 1) = R_rf(2, 2) = ca;
				R_rf(1, 2) = sa; R_rf(2, 1) = -sa;
				
				Mobs = (eye - (A * R_rf)).partialPivLu().solve(eyema) * M0;
				
				diffs[index] = sqrt(Mobs(0)*Mobs(0) + Mobs(1)*Mobs(1)) - 
				               _signal(p, i);
				index++;
			}
		}
		return 0;
	}
};

eval_type a2cSPGR;
eval_type a2cSSFP;

double calcDESPOT1(double *flipAngles, double *spgrVals, int n,
				   double TR, double B1, double *M0, double *T1);
void simplexDESPOT2(size_t nPhases, size_t *nD, double *phases, double **angles, double **ssfp,
					double TR, double T1, double B1, double *p);
#endif