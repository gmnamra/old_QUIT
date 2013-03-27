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

#include <random>

typedef Array<bool, Dynamic, Dynamic> ArrayXb;
typedef std::pair<int, double> argsort_pair;

bool argsort_comp(const argsort_pair& left, const argsort_pair& right) {
    return left.second < right.second;
}

template<typename Derived>
VectorXi arg_partial_sort(const MatrixBase<Derived> &x, int middle)
{
    VectorXi indices(middle);
    std::vector<argsort_pair> data(x.size());
    for(int i=0;i<x.size();i++)
	{
        data[i].first = i;
        data[i].second = x(i);
    }
    std::partial_sort(data.begin(), data.begin() + middle, data.end(), argsort_comp);
    for(int i=0;i<middle;i++) {
        indices(i) = data[i].first;
    }    
    return indices;
}

template <typename Functor_t, typename Derived>
ArrayXd regionContraction(ArrayBase<Derived> &params, Functor_t &f,
                          const ArrayBase<Derived> &loStart, const ArrayBase<Derived> &hiStart,
					      const int nS = 5000, const int nR = 50, const int maxContractions = 10,
						  const double thresh = 0.05, const double expand = 0., const int seed = 0)
{
	static bool haveWarned = false;
	int nP = static_cast<int>(params.size());
	MatrixXd samples(nP, nS);
	MatrixXd retained(nP, nR);
	VectorXd sampleRes(nS);
	VectorXi indices(nR);
	VectorXd loBounds = loStart;
	VectorXd hiBounds = hiStart;
	VectorXd regionSize = (hiBounds - loBounds);
	ArrayXd diffs(f.values());
	size_t c;
	
	std::mt19937 twist(seed);
	std::uniform_real_distribution<double> uniform(0., 1.);
	
	for (c = 0; c < maxContractions; c++)
	{
		for (int s = 0; s < nS; s++)
		{
			VectorXd tempSample(nP);
			do
			{
				for (int p = 0; p < nP; p++)
					tempSample(p) = uniform(twist);
				tempSample.array() *= regionSize.array();
				tempSample += loBounds;
			} while (!f.constraint(tempSample));
			f(tempSample, diffs);
			sampleRes(s) = diffs.square().sum();
			if (!std::isfinite(diffs.square().sum())) {
				if (!haveWarned) {
					haveWarned = true;
					std::cout << "Warning: Non-finite residual found!" << std::endl
							  << "Result may be meaningless. This warning will only be printed once." << std::endl;
					std::cout << "Parameters were " << tempSample.transpose() << std::endl;
					std::cout << "Signal " << f.signals().transpose() << std::endl;
					std::cout << "Theory " << f.theory(tempSample).transpose() << std::endl;
				}
				params = retained.col(0);
				diffs.setConstant(std::numeric_limits<double>::infinity());
				return diffs;
			}
			samples.col(s) = tempSample;
		}
		indices = arg_partial_sort(sampleRes, nR);
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
