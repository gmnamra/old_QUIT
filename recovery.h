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

float recovery_Mz(const float TR, const float *par, const float *c);
float recovery_dMzdM0(const float TR, const float *par, const float *c);
float recovery_dMzdalpha(const float TR, const float *par, const float *c);
float recovery_dMzdT1(const float TR, const float *par, const float *c);
float calcRecovery(float *vals, float* TR, int n, float *M0out, float *T1out, float *alpha);

#endif
