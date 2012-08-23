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

#include <iostream>
#include <Eigen/Dense>

#include "FSLIO.h"

using namespace Eigen;

//******************************************************************************
#pragma mark Convenience stuff
//******************************************************************************
double clamp(double value, double low, double high);

void linearLeastSquares(double *X, double *Y, int nD,
						double *slope, double *inter, double *res);
double classicDESPOT1(const ArrayXd &flipAngles, const ArrayXd &spgrVals,
				      double TR, double B1, double *M0, double *T1);
double classicDESPOT2(const ArrayXd &flipAngles, const ArrayXd &ssfpVals,
                      double TR, double T1, double B1, double *M0, double *T2);

template<typename Functor_t>
void write_results(const std::string outPrefix, double **paramsData, double *residualData, FSLIO *hdr)
{
	std::string outPath;
	short nx, ny, nz, nvol;
	FslGetDim(hdr, &nx, &ny, &nz, &nvol);
	FslSetDim(hdr, nx, ny, nz, 1);
	FslSetDimensionality(hdr, 3);
	FslSetDataType(hdr, NIFTI_TYPE_FLOAT32);	
	for (int p = 0; p < Functor_t::nP; p++)
	{
		outPath = outPrefix + "_" + Functor_t::names[p] + ".nii.gz";
		std::cout << "Writing parameter file: " << outPath << std::endl;
		FSLIO *outFile = FslOpen(outPath.c_str(), "wb");
		FslCloneHeader(outFile, hdr);
		FslWriteHeader(outFile);
		FslWriteVolumeFromDouble(outFile, paramsData[p], 0);
		FslClose(outFile);
	}
	
	outPath = outPrefix + "_residual.nii.gz";
	std::cout << "Writing residual file: " << outPath << std::endl;
	FSLIO *outFile = FslOpen(outPath.c_str(), "wb");
	FslCloneHeader(outFile, hdr);
	FslWriteHeader(outFile);
	FslWriteVolumeFromDouble(outFile, residualData, 0);
	FslClose(outFile);
}

#endif
