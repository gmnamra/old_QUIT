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
double SPGR(double flipAngle, double *p, double *c);
double IRSPGR(double TI, double *p, double *c);
double calcHIFI(double *flipAngles, double *spgrVals, int nSPGR, double spgrTR,
				double *TI, double *irVals, int nIR, double irFlipAngle, double irTR, double nReadout,
				double *M0, double *T1, double *B1);
void calcDESPOT1(double *flipAngles, double *spgrVals, int n,
				 double TR, double B1, double *M0, double *T1);
#endif