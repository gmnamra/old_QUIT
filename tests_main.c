/*
 *  tests_main.c
 *  Fitting
 *
 *  Created by Tobias Wood on 20/02/2012.
 *  Copyright 2012 Tobias Wood. All rights reserved.
 *
 */

#include "DESPOT.h"

//******************************************************************************
// Main
//******************************************************************************
int main(int argc, char **argv)
{
	// mcDESPOT Tests
	// 1 - Compare single component model to 2 component with same components
	
	double T1_a = 350, T1_b = 1400, T2_a = 25, T2_b = 120, tau_a = 500, f_a = 0.9;
	double TR = 5., B1 = 1.;
	double alpha[9] = {0., 10., 20., 30., 40., 50., 60., 70., 80.};
	double rf_phase = radians(180.), dO = 0.;
	
	double p1c[3] = { T2_a, dO };
	double c1c[4] = { TR, T1_a, B1, rf_phase };
	
	double p2c[7] = { T1_a, T1_b, T2_a, T2_b, f_a, tau_a, dO };
	double c2c[3] = { TR, B1, rf_phase };
	
	double s1p[3] = { 1., T1_a, B1 };
	double s1c[1] = { TR };
	
	double spgr_2c[9], ssfp_2c[9], spgr_a[9], ssfp_a[9];
	for (size_t i = 0; i < 9; i++)
	{
		alpha[i] = radians(alpha[i]);
		spgr_2c[i]  = n2cSPGR(alpha[i], p2c, c2c);
		ssfp_2c[i] = n2cSSFP(alpha[i], p2c, c2c);
	}
	a2cSPGR(alpha, p2c, c2c, spgr_a, 9);
	a2cSSFP(alpha, p2c, c2c, ssfp_a, 9);
	
	ARR_D(alpha, 9);
	ARR_D(spgr_2c, 9);
	ARR_D(spgr_a, 9);
	ARR_D(ssfp_2c, 9);
	ARR_D(ssfp_a, 9);
	fprintf(stdout, "Difference:\n");
	arraySub(ssfp_a, ssfp_2c, ssfp_a, 9);
	arraySub(spgr_a, spgr_2c, spgr_a, 9);
	ARR_D(ssfp_a, 9);
	ARR_D(spgr_a, 9);
	exit(EXIT_SUCCESS);
}