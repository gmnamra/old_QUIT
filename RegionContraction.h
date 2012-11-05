//
//  RegionContraction.h
//  DESPOT
//
//  Created by Tobias Wood on 17/08/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

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

template <typename Functor_t>
double regionContraction(VectorXd &params, Functor_t &f,
                         const VectorXd &loStart, const VectorXd &hiStart,
					     const VectorXi &loConstrained, const VectorXi &hiConstrained,
					     const int nS, const int nR, const int maxContractions,
						 const double thresh, const double expand, const int seed = 0)
{
	int nP = static_cast<int>(params.size());
	MatrixXd samples(nP, nS);
	MatrixXd retained(nP, nR);
	VectorXd sampleRes(nS);
	VectorXi indices(nR);
	VectorXd loBounds = loStart;
	VectorXd hiBounds = hiStart;
	VectorXd regionSize = (hiBounds - loBounds);
	VectorXd diffs(f.values());
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
			samples.col(s) = tempSample;
			f(samples.col(s), diffs);
			sampleRes(s) = diffs.norm();
			if (!std::isfinite(diffs.norm()))
			{
				std::cout << "Infinite Diff" << std::endl;
				std::cout << "Sample = " << samples.col(s).transpose() << std::endl;
				std::cout << "Lo Bnds= " << loBounds.transpose() << std::endl;
				std::cout << "Hi Bnds= " << hiBounds.transpose() << std::endl;
				std::cout << "Signals= " << f.signals().transpose() << std::endl;
				std::cout << "Diffs  = " << diffs.transpose() << std::endl;
				abort();
			}
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
		// but don't go past initial boundaries if constrained
		loBounds -= (regionSize * expand);
		hiBounds += (regionSize * expand);
		for (int p = 0; p < nP; p++)
		{
			if (loConstrained[p] && (loBounds[p] < loStart[p]))
				loBounds[p] = loStart[p];
			if (hiConstrained[p] && (hiBounds[p] > hiStart[p]))
				hiBounds[p] = hiStart[p];
		}
		regionSize = hiBounds - loBounds;
	}
	// Return the best evaluated solution so far
	params = retained.col(0);
	//std::cout << "Final: " << params.transpose() << " Res: " << sampleRes[indices[0]] << std::endl << std::endl;
	return sampleRes[indices[0]];
}

#endif
