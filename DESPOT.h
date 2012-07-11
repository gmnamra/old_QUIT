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
#include "mathsArray.h"
#include "mathsOptimisers.h"
int tests();
double SPGR(double flipAngle, double *p, double *c);
void SPGR_Jacobian(double *angles, int nD, double *p, double *c, double *result);
double IRSPGR(double TI, double *p, double *c);
void IRSPGR_Jacobian(double *data, int nD, double *par, double *c, double *result);
void aSSFP(double *flipAngle, double *p, double *c, double *ssfp, size_t nA);
double n2cSPGR(double alpha, double *p, double *c);
double n2cSSFP(double alpha, double *p, double *c);
void a1cSSFP(double *alpha, double *p, double *c, double *signal, size_t nA);
void a1cSSFPB0(double *alpha, double *p, double *c, double *signal, size_t nA);
void a2cSPGR(double *alpha, double *p, double *c, double *signal, size_t nA);
void a2cSSFP(double *alpha, double *p, double *c, double *signal, size_t nA);

double calcHIFI(double *flipAngles, double *spgrVals, int nSPGR, double spgrTR,
				double *TI, double *irVals, int nIR, double irFlipAngle, double irTR,
				double *M0, double *T1, double *B1);
double calcDESPOT1(double *flipAngles, double *spgrVals, int n,
				   double TR, double B1, double *M0, double *T1);
double classicDESPOT2(double *flipAngles, double *ssfpVals, int n,
                      double TR, double T1, double B1, double *p);
void simplexDESPOT2(size_t nPhases, size_t *nD, double *phases, double **angles, double **ssfp,
					double TR, double T1, double B1, double *p);
void contractDESPOT2(size_t nPhases, size_t *nD, double *phases, double **flipAngles, double **ssfp,
					 double TR, double T1, double B1, double *p);
#endif