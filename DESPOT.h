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

#ifndef __DESPOT__
#define __DESPOT__

#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

//******************************************************************************
#pragma mark Convenience stuff
//******************************************************************************
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
#endif
