/*
 *  DESPOT1.h
 *  MacRI
 *
 *  Created by Tobias Wood on 17/10/2011.
 *  Copyright 2011 Tobias Wood. All rights reserved.
 *
 */
#ifndef __DESPOT1__
#define __DESPOT1__

#include "math.h"
#include "mathArray.h"
int tests();
float SPGR(float M0, float B1, float flipAngle, float T1, float TR);
float IRSPGR(float M0, float B1, float flipAngle, float T1, float TI, float TR, float nReadout);

float calcSPGRResiduals(float *spgrVals, float *flipAngles,
					   int n, float TR, float B1,
					   float T1, float M0);
float calcIRSPGRResiduals(float *spgrVals, float *flipAngles, int nSPGR, float spgrTR,
						 float *irVals, float *TI, int nIR, float irFlipAngle, float irTR, float nReadout,
						 float B1);

void calcDESPOT1(float *spgrVals, float *flipAngles, int n, float TR, float B1, float *T1, float *M0);
float calcHIFI(float *spgrVals, float *flipAngles, int nSPGR, float spgrTR,
               float *irVals, float *TI, int nIR, float irFlipAngle, float irTR, float nReadout);


#endif