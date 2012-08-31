//
//  main.c
//  test_despot
//
//  Created by Tobias Wood on 17/08/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#include "DESPOT_Functors.h"
#include "fslio.h"
#include "procpar.h"

int main(int argc, const char * argv[])
{
	// Generate test images
	
	// Set up all the parameters
	double spgrTR = 0.013;
	double ssfpTR = 0.0105;
	double M0 = 1.;
	double B1 = 1;
	
	VectorXd alphaSPGR(8); alphaSPGR << 2, 3, 4, 5, 6, 7, 10, 14;
	VectorXd alphaSSFP(8); alphaSSFP << 11, 16, 19, 23, 27, 35, 50, 70;
	VectorXd phases(2); phases << 0., 180.;
	alphaSPGR *= M_PI / 180.;
	alphaSSFP *= M_PI / 180.;
	phases    *= M_PI / 180.;
	
	std::cout << "RF Phases: " << phases.transpose() << std::endl;
	std::cout << "SPGR flip: " << alphaSPGR.transpose() << std::endl;
	std::cout << "SSFP flip: " << alphaSSFP.transpose() << std::endl;
	
	const char *param_names[] =
	 { "T1_a", "T1_b", "T1_c", "T2_a", "T2_b", "T2_c", "f_a", "f_c", "tau_a", "B0" };
	VectorXd params(10); params << 0.465, 0.965, 3.5, 0.012, 0.09, 0.25, 0., 0., 0.125, 0.;
	VectorXd zeros(8); zeros.setZero();
	std::vector<VectorXd> ssfpZeros;
	for (int i = 0; i < phases.size(); i++)
		ssfpZeros.push_back(zeros);
	
	// Now create our images
	int nx = 5, ny = 6, nz = 2;
	int totalVoxels = nx * ny * nz;
	double *dataSPGR = (double *)malloc(alphaSPGR.size() * totalVoxels * sizeof(double));
	double *dataSSFP[phases.size()];
	for (int p = 0; p < phases.size(); p++)
		dataSSFP[p] = (double *)malloc(alphaSSFP.size() * totalVoxels * sizeof(double));
	double *dataParams[params.size()];
	for (int p = 0; p < params.size(); p++)
		dataParams[p] = (double *)malloc(totalVoxels * sizeof(double));
	double *dataM0 = (double *)malloc(totalVoxels * sizeof(double));
	double *dataB1 = (double *)malloc(totalVoxels * sizeof(double));
	
	ThreeComponent tc(alphaSPGR, zeros, alphaSSFP, phases, ssfpZeros, spgrTR, ssfpTR, M0, INFINITY, B1, false);
	VectorXd sig(tc.values());
	
	std::cout << "M0 " << M0 << " B1 " << B1 << std::endl;
	
	clock_t loopStart = clock();
	for (short z = 0; z < nz; z++)
	{	for (short y = 0; y < ny; y++)
		{	
			params[7] = ((double)y / (ny - 1)) * 0.6;
			for (short x = 0; x < nx; x++)
			{
				params[6] = ((double)x / (nx - 1)) * 0.2 + 0.05;
				//std::cout << "P: " << params.transpose() << std::endl;				
				tc(params, sig);
				//std::cout << "SPGR: " << sig.head(alphaSPGR.size()).transpose() << std::endl;
				int index = 0;
				// Results are in one long vector
				for (int i = 0; i < alphaSPGR.size(); i++)
					dataSPGR[i*totalVoxels + (z*ny + y)*nx + x] = sig[index++];
				for (int p = 0; p < phases.size(); p++)
				{
					//std::cout  << "SSFP Phase " << phases[p] << ": " << sig.segment(index, alphaSSFP.size()).transpose() << std::endl;
					for (int i = 0; i < alphaSSFP.size(); i++)
					{	dataSSFP[p][i * totalVoxels + (z*ny + y)*nx + x] = sig[index++];	}
				}
				for (int p = 0; p < params.size(); p++)
					dataParams[p][(z*ny + y)*nx + x] = params[p];
				dataM0[(z*ny + y)*nx + x] = M0;
				dataB1[(z*ny + y)*nx + x] = B1;
			}
		}
	}

	clock_t loopEnd = clock();
	std::cout << "Time to generate data was " << ((loopEnd - loopStart) / ((float)totalVoxels * CLOCKS_PER_SEC)) << " seconds per voxel." << std::endl;

	float volSize = 40.;
	// Reset to degrees
	alphaSPGR *= 180. / M_PI;
	alphaSSFP *= 180. / M_PI;
	phases *= 180. / M_PI;
	FSLIO *outFile;
	char outName[1024];
	FILE *procparfile;
	par_t *pars[3];
	double spgrPhase = 0.;
	pars[0] = createPar("tr", PAR_REAL, 1, &spgrTR);
	pars[1] = createPar("flip1", PAR_REAL, alphaSPGR.size(), alphaSPGR.data());
	pars[2] = createPar("rfphase", PAR_REAL, 1, &spgrPhase);
	snprintf(outName, 1024, "SPGR_Noiseless.nii.gz");
	outFile = FslOpen(outName, "wb");
	FslSetDim(outFile, nx, ny, nz, alphaSPGR.size());
	FslSetVoxDim(outFile, volSize / nx, volSize / ny, volSize / nz, spgrTR);
	FslWriteHeader(outFile);
	FslWriteAllVolumesFromDouble(outFile, dataSPGR);
	FslClose(outFile);
	snprintf(outName, 1024, "SPGR_Noiseless.procpar");
	procparfile = fopen(outName, "w");
	for (int i = 0; i < 3; i++)
		fprintPar(procparfile, pars[i]);
	fclose(procparfile);
	for (int p = 0; p < phases.size(); p++)
	{	
		double phaseDeg = phases[p];
		snprintf(outName, 1024, "SSFP-%.00f_Noiseless.nii.gz", phaseDeg);
		outFile = FslOpen(outName, "wb");
		FslSetDim(outFile, nx, ny, nz, alphaSSFP.size());
		FslSetVoxDim(outFile, volSize / nx, volSize / ny, volSize / nz, spgrTR);
		FslWriteHeader(outFile);
		FslWriteAllVolumesFromDouble(outFile, dataSSFP[p]);
		FslClose(outFile);
		snprintf(outName, 1024, "SSFP-%.00f_Noiseless.procpar", phaseDeg);
		setPar(pars[0], 1, &ssfpTR);
		setPar(pars[1], alphaSSFP.size(), alphaSSFP.data());
		setPar(pars[2], 1, &phaseDeg);
		procparfile = fopen(outName, "w");
		for (int i = 0; i < 3; i++)
			fprintPar(procparfile, pars[i]);
		fclose(procparfile);
	}
	
	for (int p = 0; p < params.size(); p++)
	{	
		snprintf(outName, 1024, "%s_Noiseless.nii.gz", param_names[p]);
		outFile = FslOpen(outName, "wb");
		FslSetDim(outFile, nx, ny, nz, 1);
		FslSetVoxDim(outFile, volSize / nx, volSize / ny, volSize / nz, 0.);
		FslWriteHeader(outFile);
		FslWriteAllVolumesFromDouble(outFile, dataParams[p]);
		FslClose(outFile);
	}
	snprintf(outName, 1024, "M0_Noiseless.nii.gz");
	outFile = FslOpen(outName, "wb");
	FslSetDim(outFile, nx, ny, nz, 1);
	FslSetVoxDim(outFile, volSize / nx, volSize / ny, volSize / nz, 0.);
	FslWriteHeader(outFile);
	FslWriteAllVolumesFromDouble(outFile, dataM0);
	FslClose(outFile);

	snprintf(outName, 1024, "B1_Noiseless.nii.gz");
	outFile = FslOpen(outName, "wb");
	FslSetDim(outFile, nx, ny, nz, 1);
	FslSetVoxDim(outFile, volSize / nx, volSize / ny, volSize / nz, 0.);
	FslWriteHeader(outFile);
	FslWriteAllVolumesFromDouble(outFile, dataB1);
	FslClose(outFile);
		
    return 0;
}

