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
#pragma mark Versioning etc.
//******************************************************************************
#include "version"

const string credit_shared {
"Written by tobias.wood@kcl.ac.uk, based on work by Sean Deoni. \n\
Acknowledgements greatfully received, grant discussions welcome."
};

const string credit_me {
"Written by tobias.wood@kcl.ac.uk. \n\
Acknowledgements greatfully received, grant discussions welcome."
};

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
typedef Matrix<double, 6, 6> Matrix6d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 9, 9> Matrix9d;
typedef Matrix<double, 9, 1> Vector9d;
typedef Matrix<double, 3, Dynamic> MagVector;

const VectorXd SigMag(const MagVector &M_in);
const MagVector SumMC(const MatrixXd &M_in);

const Matrix3d RF(const double &alpha, const double &beta);
inline const Matrix3d Relax(const double &T1, const double &T2);
inline const Matrix3d InfinitesimalRF(const double &dalpha);
inline const Matrix3d OffResonance(const double &inHertz);
inline const Matrix3d Spoiling();
inline const Matrix6d Exchange(const double &k_ab, const double &k_ba);
const void CalcExchange(const double tau_a, const double f_a, const double f_b, double &k_ab, double &k_ba);

//******************************************************************************
#pragma mark Actual Signal Equations
// Parameters are { T1, T2, f0 }
//******************************************************************************
MagVector One_SPGR(const VectorXd &p, const ArrayXd &flip, const double TR);
MagVector One_SSFP(const VectorXd &p, const ArrayXd &flip, const double TR, const double ph);
MagVector One_SSFP_Finite(const VectorXd &p, const ArrayXd &flip, const bool spoil,
                          const double TR, const double Trf, const double TE,
						  const double ph);
//******************************************************************************
// Parameters are { T1_a, T2_a, T1_b, T2_b, tau_a, f_a, f0 }
//******************************************************************************
MagVector Two_SPGR(const VectorXd &p, const ArrayXd &flip, const double TR);
MagVector Two_SSFP(const VectorXd &p, const ArrayXd &flip, const double TR, const double ph);
MagVector Two_SSFP_Finite(const VectorXd &p, const ArrayXd &flip, const bool spoil,
                          const double TR, const double Trf, const double TE,
						  const double ph);
//******************************************************************************
// Parameters are { T1a, T2a, T1b, T2b, T1c, T2c, tau_a, f_a, f_c, f0 }
//******************************************************************************
MagVector Three_SPGR(const VectorXd &p, const ArrayXd &flip, const double TR);
MagVector Three_SSFP(const VectorXd &p, const ArrayXd &flip, const double TR, const double ph);
MagVector Three_SSFP_Finite(const VectorXd &p, const ArrayXd &flip, const bool spoil,
                            const double TR, const double Trf, const double TE,
						    const double ph);
#endif
