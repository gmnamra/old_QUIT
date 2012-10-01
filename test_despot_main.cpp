//
//  main.c
//  test_despot
//
//  Created by Tobias Wood on 17/08/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#include <getopt.h>
#include <stdio.h>
#include "DESPOT_Functors.h"
#include "NiftiImage.h"
#include "procpar.h"

static int components = 1;
static std::string outPrefix;
static struct option long_options[] =
{
	{"components", required_argument, 0, 'c'},
	{0, 0, 0, 0}
};

int main(int argc, char **argv)
{
	int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "c:", long_options, &indexptr)) != -1)
	{
		switch (c)
		{
			case 'c':
				components = atoi(optarg);
				outPrefix = std::string("test_") + optarg + "c";
				if ((components < 1) || (components > 3))
				{
					std::cout << "Valid models are 1 to 3 components." << std::endl;
					exit(EXIT_FAILURE);
				}
				break;
			default:
				std::cout << "Unknown command line switch " << c << std::endl;
		}
	}

	// Generate test images
	
	// Set up all the parameters
	double spgrTR = 0.0087;
	double ssfpTR = 0.002916;
	double B0 = 0.;
	double B1 = 1;
	
	VectorXd alphaSPGR(6); alphaSPGR << 2, 4, 8, 16, 24, 32;
	VectorXd alphaSSFP(6); alphaSSFP << 2, 4, 8, 16, 32, 64;
	VectorXd phases(3); phases << 180., 135., 225.;
	alphaSPGR *= M_PI / 180.;
	alphaSSFP *= M_PI / 180.;
	phases    *= M_PI / 180.;
	VectorXd zeros(6); zeros.setZero();
	std::vector<VectorXd> ssfpZeros;
	for (int i = 0; i < phases.size(); i++)
		ssfpZeros.push_back(zeros);
	
	std::cout << "RF Phases: " << phases.transpose() << std::endl;
	std::cout << "SPGR flip: " << alphaSPGR.transpose() << std::endl;
	std::cout << "SSFP flip: " << alphaSSFP.transpose() << std::endl;
	
	VectorXd params;
	Functor<double> *model = NULL;
	switch (components)
	{
		case 1:
			params.resize(OneComponent::nP, 1);
			params << 1.e5, 1., 0.015;
			break;
		case 2:
			params.resize(TwoComponent::nP, 1);
			params << 0.7, 2., 0.008, 0.02, 0., 0.125, 0.;
			model = new TwoComponent(alphaSPGR, zeros, alphaSSFP, phases, ssfpZeros, spgrTR, ssfpTR, 1., 0., B1);
			break;
		case 3:
			params.resize(ThreeComponent::nP, 1);
			params << 0.7, 2., 8., 0.008, 0.02, 0.25, 0., 0., 0.125, 0.;
			model = new ThreeComponent(alphaSPGR, zeros, alphaSSFP, phases, ssfpZeros, spgrTR, ssfpTR, 1., INFINITY, B1);
			break;
	}
	
	// Now create our images
	int nx = 3, ny = 3, nz = 3;
	int totalVoxels = nx * ny * nz;
	double *dataSPGR = (double *)malloc(alphaSPGR.size() * totalVoxels * sizeof(double));
	double *dataSSFP[phases.size()];
	for (int p = 0; p < phases.size(); p++)
		dataSSFP[p] = (double *)malloc(alphaSSFP.size() * totalVoxels * sizeof(double));
	double *dataParams[params.size()];
	for (int p = 0; p < params.size(); p++)
		dataParams[p] = (double *)malloc(totalVoxels * sizeof(double));
	double *dataB0 = (double *)malloc(totalVoxels * sizeof(double));
	double *dataB1 = (double *)malloc(totalVoxels * sizeof(double));
	
	clock_t loopStart = clock();
	for (short z = 0; z < nz; z++)
	{	for (short y = 0; y < ny; y++)
		{	for (short x = 0; x < nx; x++)
			{
				switch (components)
				{
					case 1:
						B0 = ((double)x / (nx - 1)) * 80. - 40.;
						B1 = ((double)y / (ny - 1)) * 0.2 + 0.9;
						model = new OneComponent(alphaSPGR, zeros, alphaSSFP, phases, ssfpZeros, spgrTR, ssfpTR, B0, B1);
						break;
					case 3:
						params[6] = ((double)x / (nx - 1)) * 0.2 + 0.05;
						params[7] = ((double)y / (ny - 1)) * 0.6;
						break;
				}
				VectorXd sig(model->values());
				model->operator()(params, sig);
				int index = 0;
				// Results are in one long vector
				for (int i = 0; i < alphaSPGR.size(); i++)
					dataSPGR[i*totalVoxels + (z*ny + y)*nx + x] = sig[index++];
				for (int p = 0; p < phases.size(); p++)
				{
					for (int i = 0; i < alphaSSFP.size(); i++)
					{	dataSSFP[p][i * totalVoxels + (z*ny + y)*nx + x] = sig[index++];	}
				}
				for (int p = 0; p < params.size(); p++)
					dataParams[p][(z*ny + y)*nx + x] = params[p];
				dataB0[(z*ny + y)*nx + x] = B0;
				dataB1[(z*ny + y)*nx + x] = B1;
				delete model;
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
	NiftiImage outFile(nx, ny, nz, alphaSPGR.size(),
					   volSize / nx, volSize / ny, volSize / nz, spgrTR,
					   NIFTI_TYPE_FLOAT32);
	std::string outName;
	FILE *procparfile;
	par_t *pars[3];
	double spgrPhase = 0.;
	pars[0] = createPar("tr", PAR_REAL, 1, &spgrTR);
	pars[1] = createPar("flip1", PAR_REAL, alphaSPGR.size(), alphaSPGR.data());
	pars[2] = createPar("rfphase", PAR_REAL, 1, &spgrPhase);
	outName = "SPGR_Noiseless.nii.gz";
	outFile.open(outName, NIFTI_WRITE);
	outFile.writeAllVolumes(dataSPGR);
	outFile.close();
	outName = "SPGR_Noiseless.procpar";
	procparfile = fopen(outName.c_str(), "w");
	for (int i = 0; i < 3; i++)
		fprintPar(procparfile, pars[i]);
	fclose(procparfile);
	for (int p = 0; p < phases.size(); p++)
	{	
		double phaseDeg = phases[p];
		outName = (std::stringstream("SSFP-") << outName << phaseDeg << "_Noiseless.nii.gz").str();
		outFile.setnt(alphaSSFP.size());
		outFile.open(outName, NIFTI_WRITE);
		outFile.writeAllVolumes(dataSSFP[p]);
		outFile.close();
		outName = (std::stringstream("SSFP-") << outName << phaseDeg << "_Noiseless.procpar").str();
		setPar(pars[0], 1, &ssfpTR);
		setPar(pars[1], alphaSSFP.size(), alphaSSFP.data());
		setPar(pars[2], 1, &phaseDeg);
		procparfile = fopen(outName.c_str(), "w");
		for (int i = 0; i < 3; i++)
			fprintPar(procparfile, pars[i]);
		fclose(procparfile);
	}
	
	outName = "B0_Noiseless.nii.gz";
	outFile.setnt(1);
	outFile.open(outName, NIFTI_WRITE);
	outFile.writeAllVolumes(dataB0);
	outFile.close();
	outName = "B1_Noiseless.nii.gz";
	outFile.open(outName, NIFTI_WRITE);
	outFile.writeAllVolumes(dataB1);
	outFile.close();
	switch (components)
	{
		case 1: write_results<OneComponent>(outPrefix, dataParams, NULL, outFile); break;
		case 2: write_results<TwoComponent>(outPrefix, dataParams, NULL, outFile); break;
		case 3: write_results<ThreeComponent>(outPrefix, dataParams, NULL, outFile); break;
	}
    return 0;
}

