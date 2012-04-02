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

#include "math.h"
#include "mathArray.h"
#include "mathOptimisers.h"
int tests();
float SPGR(float flipAngle, float *p, float *c);
void SPGR_Jacobian(float *angles, int nD, float *p, float *c, float *result);
float IRSPGR(float TI, float *p, float *c);
void IRSPGR_Jacobian(float *data, int nD, float *par, float *c, float *result);
float SSFP(float flipAngle, float *p, float *c);
float n2cSPGR(float alpha, float *p, float *c);
float n2cSSFP(float alpha, float *p, float *c);
void a1cSSFP(float *alpha, float *p, float *c, float *signal, size_t nA);
void a2cSPGR(float *alpha, float *p, float *c, float *signal, size_t nA);
void a2cSSFP(float *alpha, float *p, float *c, float *signal, size_t nA);

float calcHIFI(float *flipAngles, float *spgrVals, int nSPGR, float spgrTR,
				float *TI, float *irVals, int nIR, float irFlipAngle, float irTR,
				float *M0, float *T1, float *B1);
void calcDESPOT1(float *flipAngles, float *spgrVals, int n,
				 float TR, float B1, float *M0, float *T1);
void classicDESPOT2(float *flipAngles, float *ssfpVals, int n,
                    float TR, float T1, float B1, float *p);
void simplexDESPOT2(size_t nPhases, size_t *nD, float *phases, float **angles, float **ssfp,
					float TR, float T1, float B1, float *p);
void contractDESPOT2(size_t nPhases, size_t *nD, float *phases, float **flipAngles, float **ssfp,
					 float TR, float T1, float B1, float *p);
#endif