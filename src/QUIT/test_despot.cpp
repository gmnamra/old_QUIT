//
//  main.c
//  test_despot
//
//  Created by Tobias Wood on 17/08/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#include <getopt.h>

#include <vector>
#include <string>

#include "Nifti/Nifti.h"
#include "DESPOT_Functors.h"

static int components = 1, tesla = 7;
static std::string outPrefix;
static bool testSpeed = false, testFiles = false;
static struct option long_options[] =
{
	{"components", required_argument, 0, 'c'},
	{"speed", no_argument, 0, 's'},
	{"files", no_argument, 0, 'f'},
	{"tesla", required_argument, 0, 't'},
	{0, 0, 0, 0}
};

using namespace std;

int main(int argc, char **argv)
{
	/*int indexptr = 0, c;
	while ((c = getopt_long(argc, argv, "c:sf", long_options, &indexptr)) != -1) {
		switch (c) {
			case 'c':
				components = atoi(optarg);
				outPrefix = std::string("test_") + optarg + "c";
				if ((components < 1) || (components > 3)) {
					std::cout << "Valid models are 1 to 3 components." << std::endl;
					return EXIT_FAILURE;
				}
				break;
			case 's': testSpeed = true; break;
			case 'f': testFiles = true; break;
			case 't': tesla = atoi(optarg);
			default:
				break;
		}
	}*/
		
	// Set up all the parameters	
	ArrayXd alphaSPGR(8); alphaSPGR << 3, 4, 5, 6, 7, 9, 13, 18;
	ArrayXd alphaSSFP(8); alphaSSFP << 12.4, 16.5, 21.6, 27.8, 34.0, 41.2, 52.5, 70.0;
	ArrayXd phases(2); phases << M_PI, 0;
	alphaSPGR *= M_PI / 180.;
	alphaSSFP *= M_PI / 180.;
	VectorXd p(11); p << 0.127495, 0.0174758,   1.82927,  0.101974,    3.8992,   3.32707,   0.27656,  0.151555,  0.848445,   22.7945,   94.2626;
	
	SPGRFinite spgr(Pools::Three, alphaSPGR, 0.0081, 0.0005, 0.0037);
	SSFPFinite ssfp(Pools::Three, alphaSSFP, 0.003912, 0.0005, phases);
	auto spgr_s = spgr.signal(p.head(10), 1.0, p(10));
	auto ssfp_s = ssfp.signal(p.head(10), 1.0, p(10));
	cout << spgr_s.transpose() << endl;
	cout << (spgr_s.array() / spgr_s.mean()).transpose() << endl;
	cout << ssfp_s.transpose() << endl;
	cout << (ssfp_s.array() / ssfp_s.mean()).transpose() << endl;
	//angles.push_back(alphaSSFP); signals.push_back(sSSFP); consts.emplace_back(false, ssfpTR, Trf, 0., M_PI);
	
	/*long loops = 1;
	if (testSpeed) {
		cout << "Alpha SPGR: " << alphaSPGR.transpose() << endl;
		cout << "Alpha SSFP: " << alphaSSFP.transpose() << endl;
		cout << "SPGR TR: " << spgrTR << " SSFP TR: " << ssfpTR << " Trf: " << Trf << endl;
		for (int c = 1; c < 4; c++) {
			ArrayXd p(mcDESPOT::nP(c) + mcDESPOT::nB0(mcDESPOT::B0_Single, signals.size()));
			p.head(mcDESPOT::nP(c)) = (mcDESPOT::defaultLo(c, tesla) + mcDESPOT::defaultHi(c, tesla)) / 2.;
			p.tail(mcDESPOT::nB0(mcDESPOT::B0_Single, signals.size())).setConstant(0.);
			mcDESPOT mcd(c, data, mcDESPOT::B0_Single, false, false);
			VectorXd signal;
			clock_t start = clock();
			for (int l = 0; l < loops; l++)
				signal = mcd.theory(p);
			clock_t end = clock();
			cout << c << "-component instant model average time " << ((end - start) / ((float)loops * CLOCKS_PER_SEC)) * 1000 << " ms" << endl;
			cout << "Params: " << p.transpose() << endl;
			cout << "Signal: " << signal.transpose() << endl;
		}
		for (int c = 1; c < 4; c++) {
			ArrayXd p(mcDESPOT::nP(c) + mcDESPOT::nB0(mcDESPOT::B0_Single, signals.size()));
			p.head(mcDESPOT::nP(c)) = (mcDESPOT::defaultLo(c, tesla) + mcDESPOT::defaultHi(c, tesla)) / 2.;
			p.tail(mcDESPOT::nB0(mcDESPOT::B0_Single, signals.size())).setConstant(0.);
			mcFinite mcd(c, data, mcDESPOT::B0_Single, false, false);
			VectorXd signal;
			clock_t start = clock();
			for (int l = 0; l < loops; l++)
				signal = mcd.theory(p);
			clock_t end = clock();
			cout << c << "-component finite model average time " << ((end - start) / ((float)loops * CLOCKS_PER_SEC)) * 1000 << " ms" << endl;
			cout << "Params: " << p.transpose() << endl;
			cout << "Signal: " << signal.transpose() << endl;
		}
	}
	
	/*
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
	Nifti outFile(nx, ny, nz, alphaSPGR.size(),
					   volSize / nx, volSize / ny, volSize / nz, spgrTR,
					   NIFTI_TYPE_FLOAT32);
	std::string outName;
	FILE *procparfile;
	par_t *pars[3];
	double spgrPhase = 0.;
	pars[0] = createPar("tr", PAR_REAL, 1, &spgrTR);
	pars[1] = createPar("flip1", PAR_REAL, alphaSPGR.size(), alphaSPGR.data());
	pars[2] = createPar("rfphase", PAR_REAL, 1, &spgrPhase);
	outName = "SPGR_Noiseless" + OutExt();
	outFile.open(outName, WRITE);
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
		outName = (std::stringstream("SSFP-") << outName << phaseDeg << "_Noiseless" + OutExt()).str();
		outFile.setnt(alphaSSFP.size());
		outFile.open(outName, WRITE);
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
	
	outName = "B0_Noiseless" + OutExt();
	outFile.setnt(1);
	outFile.open(outName, WRITE);
	outFile.writeAllVolumes(dataB0);
	outFile.close();
	outName = "B1_Noiseless" + OutExt();
	outFile.open(outName, WRITE);
	outFile.writeAllVolumes(dataB1);
	outFile.close();
	switch (components)
	{
		case 1: write_results<OneComponent>(outPrefix, dataParams, NULL, outFile); break;
		case 2: write_results<TwoComponent>(outPrefix, dataParams, NULL, outFile); break;
		case 3: write_results<ThreeComponent>(outPrefix, dataParams, NULL, outFile); break;
	}
	*/
    return 0;
}

