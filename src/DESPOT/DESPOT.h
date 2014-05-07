/*
 *  DESPOT.h
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright (c) 2011-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef DESPOT_DESPOT
#define DESPOT_DESPOT

#include <iostream>
#include <exception>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#ifdef AGILENT
#include "Nifti/Nifti.h"
#include "Nifti/ExtensionCodes.h"
#include "procpar.h"
#endif

using namespace std;
using namespace Eigen;

//******************************************************************************
#pragma mark Convenience stuff
//******************************************************************************
#ifdef AGILENT
bool ReadPP(const Nifti &nii, Agilent::ProcPar &pp);
#endif

double clamp(double value, double low, double high);

void linearLeastSquares(const ArrayXd &X, const ArrayXd &Y,
						double &slope, double &inter, double &res);


double classicDESPOT1(const ArrayXd &flipAngles, const ArrayXd &spgrVals,
				      double TR, double B1, double &M0, double &T1);
double classicDESPOT2(const ArrayXd &flipAngles, const ArrayXd &ssfpVals,
                      double TR, double T1, double B1, double &M0, double &T2);

//******************************************************************************
#pragma mark Basic Signal Equations
//******************************************************************************
ArrayXd SPGR(const ArrayXd &flip, const double &TR, const double &B1, const double &M0, const double &T1);
ArrayXd IRSPGR(const ArrayXd &TI, const double &TR, const double &B1,
               const double &flip, const double &eff,
			   const double &M0, const double &T1);

//******************************************************************************
#pragma Magnetisation Evolution Matrices, helper functions etc.
//******************************************************************************
typedef const double cdbl; // To save tedious typing

typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 9, 9> Matrix9d;
typedef Matrix<double, 9, 1> Vector9d;
typedef Matrix<double, 3, Dynamic> MagVector;

const VectorXd SigMag(const MagVector &M_in);
const VectorXcd SigComplex(const MagVector &M_in);
const MagVector SumMC(const MatrixXd &M_in);

const Matrix3d RF(cdbl &alpha, cdbl &beta);
inline const Matrix3d Relax(cdbl &T1, cdbl &T2);
inline const Matrix3d InfinitesimalRF(cdbl &dalpha);
inline const Matrix3d OffResonance(cdbl &inHertz);
inline const Matrix3d Spoiling();
inline const Matrix6d Exchange(cdbl &k_ab, cdbl &k_ba);
const void CalcExchange(cdbl tau_a, cdbl f_a, cdbl f_b, double &k_ab, double &k_ba);

//******************************************************************************
#pragma mark Actual Signal Equations
// Parameters are { T1, T2, f0 }
//******************************************************************************
MagVector One_SPGR(const ArrayXd &flip, cdbl TR, cdbl PD, cdbl T1);
MagVector One_SSFP(const ArrayXd &flip, cdbl TR, cdbl ph, cdbl PD, cdbl T1, cdbl T2, cdbl f0);
MagVector One_SSFP_Finite(const ArrayXd &flip, const bool spoil, cdbl TR, cdbl Trf, cdbl TE, cdbl ph,
                          cdbl PD, cdbl T1, cdbl T2, cdbl f0);
//******************************************************************************
// Parameters are { T1_a, T2_a, T1_b, T2_b, tau_a, f_a, f0 }
//******************************************************************************
MagVector Two_SPGR(const ArrayXd &flip, cdbl TR, cdbl PD, cdbl T1_a, cdbl T1_b, cdbl tau_a, cdbl f_a);
MagVector Two_SSFP(const ArrayXd &flip, cdbl TR, cdbl ph, cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b, cdbl tau_a, cdbl f_a, cdbl f0);
MagVector Two_SSFP_Finite(const ArrayXd &flip, const bool spoil, cdbl TR, cdbl Trf, cdbl TE, cdbl ph,
                          cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b,
						  cdbl tau_a, cdbl f_a, cdbl f0);
//******************************************************************************
// Parameters are { T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c, f0 }
//******************************************************************************
MagVector Three_SPGR(const ArrayXd &flip, cdbl TR, cdbl PD, cdbl T1_a, cdbl T1_b, cdbl T1_c, cdbl tau_a, cdbl f_a, cdbl f_c);
MagVector Three_SSFP(const ArrayXd &flip, cdbl TR, cdbl ph, cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b, cdbl T1_c, cdbl T2_c, cdbl tau_a, cdbl f_a, cdbl f_c, cdbl f0);
MagVector Three_SSFP_Finite(const ArrayXd &flip, const bool spoil, cdbl TR, cdbl Trf, cdbl TE, cdbl ph,
                            cdbl PD, cdbl T1_a, cdbl T2_a, cdbl T1_b, cdbl T2_b, cdbl T1_c, cdbl T2_c,
							cdbl tau_a, cdbl f_a, cdbl f_c, cdbl f0);
#endif
