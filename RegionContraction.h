/*
 *  RegionContraction.h
 *
 *  Created by Tobias Wood on 17/08/2012.
 *  Copyright (c) 2012-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef DESPOT_RegionContraction_h
#define DESPOT_RegionContraction_h

#include <vector>
#include <random>
#include <iostream>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

typedef Array<bool, Dynamic, Dynamic> ArrayXb;

template<typename Derived>
vector<size_t> index_partial_sort(const ArrayBase<Derived> &x, size_t N)
{
<<<<<<< HEAD
	eigen_assert(x.size() >= N);
    vector<size_t> allIndices(x.size()), indices(N);
    for(int i=0;i<allIndices.size();i++) {
		allIndices[i] = i;
    }
	partial_sort(allIndices.begin(), allIndices.begin() + N, allIndices.end(),
	             [&x](size_t i1, size_t i2) { return x[i1] < x[i2]; });
	for (size_t i = 0; i < N; i++) {
		indices[i] = allIndices[i];
	}
    return indices;
}

template <typename Functor_t, typename Derived>
ArrayXd regionContraction(ArrayBase<Derived> &params, Functor_t &f,
                          const ArrayBase<Derived> &loStart, const ArrayBase<Derived> &hiStart,
					      const int nS = 5000, const int nR = 50, const int maxContractions = 10,
						  const double thresh = 0.05, const double expand = 0., const int seed = 0)
{
	static bool finiteWarning = false;
	static bool constraintWarning = false;
	int nP = static_cast<int>(params.size());
	ArrayXXd samples(nP, nS);
	ArrayXXd retained(nP, nR);
	ArrayXd sampleRes(nS);
	vector<size_t> indices(nR);
	ArrayXd loBounds = loStart;
	ArrayXd hiBounds = hiStart;
	ArrayXd regionSize = (hiBounds - loBounds);
	ArrayXd diffs(f.values()); diffs.setZero();
	size_t c;
	
	mt19937 twist(seed);
	uniform_real_distribution<double> uniform(0., 1.);
	
	for (c = 0; c < maxContractions; c++) {
		for (int s = 0; s < nS; s++) {
			ArrayXd tempSample(nP);
			size_t nTries = 0;
			do {
				for (int p = 0; p < nP; p++)
					tempSample(p) = uniform(twist);
				tempSample.array() *= regionSize.array();
				tempSample += loBounds;
				nTries++;
				if (nTries > 100) {
					if (!constraintWarning) {
						constraintWarning = true;
						cout << "Warning: Cannot fulfill sample constraints after " << to_string(nTries) << " attempts, giving up." << endl;
						cout << "Last attempt was: " << tempSample.transpose() << endl;
						cout << "This warning will only be printed once." << endl;
					}
					params.setZero();
					return diffs;
				}
			} while (!f.constraint(tempSample));
			f(tempSample, diffs);
			sampleRes(s) = diffs.square().sum();
			if (!isfinite(diffs.square().sum())) {
				if (!finiteWarning) {
					finiteWarning = true;
					cout << "Warning: Non-finite residual found!" << endl
							  << "Result may be meaningless. This warning will only be printed once." << endl;
					cout << "Parameters were " << tempSample.transpose() << endl;
					cout << "Signal " << f.signals().transpose() << endl;
					cout << "Theory " << f.theory(tempSample).transpose() << endl;
				}
				params = retained.col(0);
				diffs.setConstant(numeric_limits<double>::infinity());
				return diffs;
			}
			samples.col(s) = tempSample;
		}
		indices = index_partial_sort(sampleRes, nR);
		for (int i = 0; i < nR; i++)
			retained.col(i) = samples.col(indices[i]);
		
		// Find the min and max for each parameter in the top nR samples
		loBounds = retained.rowwise().minCoeff();
		hiBounds = retained.rowwise().maxCoeff();
		regionSize = (hiBounds - loBounds);	
		// Terminate if ALL the distances between bounds are under the threshold
		if (((regionSize.array() / hiBounds.array()).abs() < thresh).all())
			break;
		
		// Expand the boundaries back out in case we just missed a minima,
		// but don't go past initial boundaries
		loBounds -= (regionSize * expand);
		hiBounds += (regionSize * expand);
		for (int p = 0; p < nP; p++)
		{
			if (loBounds[p] < loStart[p])
				loBounds[p] = loStart[p];
			if (hiBounds[p] > hiStart[p])
				hiBounds[p] = hiStart[p];
		}
		regionSize = hiBounds - loBounds;
	}
	// Return the best evaluated solution so far
	params = retained.col(0);
	// Calculate the residual in %
	f(params, diffs);
	diffs /= f.signals();
	return diffs;
}

#endif
