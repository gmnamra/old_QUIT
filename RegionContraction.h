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
#include <atomic>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

typedef Array<bool, Dynamic, Dynamic> ArrayXb;

vector<size_t> index_partial_sort(const Ref<ArrayXd> &x, size_t N)
{
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

template <typename Functor_t>
ArrayXd regionContraction(Ref<ArrayXd> params, Functor_t &f,
                          const Ref<ArrayXXd> &startBounds, const Ref<ArrayXd> &weights,
					      const int nS = 5000, const int nR = 50, const int maxContractions = 10,
						  const double thresh = 0.05, const double expand = 0., const int seed = 0)
{
	eigen_assert(params.size() == startBounds.rows());
	eigen_assert(startBounds.cols() == 2);
	
	static atomic<bool> finiteWarning(false);
	static atomic<bool> constraintWarning(false);
	int nP = static_cast<int>(params.size());
	ArrayXXd samples(nP, nS);
	ArrayXXd retained(nP, nR);
	ArrayXd sampleRes(nS);
	vector<size_t> indices(nR);
	ArrayXXd bounds = startBounds;
	ArrayXd regionSize = (bounds.col(1) - bounds.col(0));
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
				tempSample += bounds.col(0);
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
			sampleRes(s) = (diffs * weights).square().sum();
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
		bounds.col(0) = retained.rowwise().minCoeff();
		bounds.col(1) = retained.rowwise().maxCoeff();
		regionSize = (bounds.col(1) - bounds.col(0));
		// Terminate if ALL the distances between bounds are under the threshold
		if (((regionSize.array() / bounds.col(1)).abs() < thresh).all())
			break;
		
		// Expand the boundaries back out in case we just missed a minima,
		// but don't go past initial boundaries
		bounds.col(0) = (bounds.col(0) - regionSize * expand).max(startBounds.col(0));
		bounds.col(1) = (bounds.col(1) + regionSize * expand).min(startBounds.col(1));
		regionSize = bounds.col(1) - bounds.col(0);
	}
	// Return the best evaluated solution so far
	params = retained.col(0);
	// Calculate the residuals
	f(params, diffs);
	//diffs /= f.signals();
	return diffs;
}

#endif
