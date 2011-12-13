/*
 *  recovery.h
 *  MacRI
 *
 *  Created by Tobias Wood on 07/11/2011.
 *  Copyright 2011 Tobias Wood. All rights reserved.
 *
 */

#ifndef RECOVERY
#define RECOVERY

#include <math.h>

double recovery_Mz(const double TR, const double *par, const double *c);
double recovery_dMzdM0(const double TR, const double *par, const double *c);
double recovery_dMzdalpha(const double TR, const double *par, const double *c);
double recovery_dMzdT1(const double TR, const double *par, const double *c);
double calcRecovery(double *vals, double* TR, int n, double *M0out, double *T1out, double *alpha);

#endif
