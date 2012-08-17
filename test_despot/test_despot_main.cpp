//
//  main.c
//  test_despot
//
//  Created by Tobias Wood on 17/08/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#include "DESPOT_Functors.h"

int main(int argc, const char * argv[])
{
	VectorXd alphaSPGR(8); alphaSPGR << 2, 2.5, 3, 4, 5, 7.5, 10, 15;
	VectorXd alphaSSFP(8); alphaSSFP << 2, 4, 6, 7, 12, 16, 24, 32;
	
	alphaSPGR *= M_PI / 180.;
	alphaSSFP *= M_PI / 180.;
	
	VectorXd params(7); params << 0.8, 2., 0.02, 0.03, 0.5, 0.5, 0.;
	VectorXd phases(1); phases << M_PI;
	VectorXd zeros(8); zeros.setZero();
	std::vector<VectorXd> ssfpZeros;
	ssfpZeros.push_back(zeros);
	
	double spgrTR = 0.0087;
	double ssfpTR = 0.002916;
	double M0 = 1.e6;
	double B1 = 1;
	TwoComponent tc(alphaSPGR, zeros, alphaSSFP, phases, ssfpZeros, spgrTR, ssfpTR, M0, B1);
	VectorXd sig(tc.values());
	clock_t loopStart = clock();
	int nLoop = 5000;
	for (int i = 0; i < nLoop; i++)
	{
		tc(params, sig);
	}
	clock_t loopEnd = clock();
	fprintf(stdout, "Combined time per loop was %f s.\n", (loopEnd - loopStart) / ((double)nLoop * CLOCKS_PER_SEC));
	
	
	std::cout << "spgr: " << sig.head(8).transpose() << std::endl;
	std::cout << "ssfp: " << sig.tail(8).transpose() << std::endl;
    return 0;
}

