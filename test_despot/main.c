//
//  main.c
//  test_despot
//
//  Created by Tobias Wood on 17/08/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#include <stdio.h>
#include "DESPOT.h"

int main(int argc, const char * argv[])
{
	SPGR_constants cSPGR;
	SSFP_constants cSSFP;
	
	cSPGR.B1 = 1.;
	cSPGR.M0 = 1.e6;
	cSPGR.TR = 0.0087;
	
	cSSFP.B0 = 0.;
	cSSFP.B1 = 1.;
	cSSFP.M0 = 1.e6;
	cSSFP.rfPhase = M_PI;
	cSSFP.TR = 0.002916;
	
	int nSPGR = 8;
	int nSSFP = 8;
	
	double alphaSPGR[8] = {2, 2.5, 3, 4, 5, 7.5, 10, 15};
	double alphaSSFP[8] = {2, 4, 6, 7, 12, 16, 24, 32};
	
	arrayApply(alphaSPGR, alphaSPGR, radians, nSPGR);
	arrayApply(alphaSSFP, alphaSSFP, radians, nSSFP);
	
	double params[7] = { 0.8, 2., 0.02, 0.03, 0.5, 0.5 };
	double sigSPGR[8], sigSSFP[8];
	
	a2cSPGR(alphaSPGR, params, &cSPGR, sigSPGR, nSPGR);
	a2cSSFP(alphaSSFP, params, &cSSFP, sigSSFP, nSSFP);
	
	ARR_D(sigSPGR, nSPGR);
	ARR_D(sigSSFP, nSSFP);
    return 0;
}

