//
//  RegionContraction.h
//  DESPOT
//
//  Created by Tobias Wood on 17/08/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#ifndef DESPOT_RegionContraction_h
#define DESPOT_RegionContraction_h

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
						 const double thresh, const double expand)
{
	size_t nP = params.size();
	
	MatrixXd samples(nP, nS);
	MatrixXd retained(nP, nR);
	VectorXd sampleRes(nS);
	VectorXi indices(nR);
	VectorXd loBounds = loStart;
	VectorXd hiBounds = hiStart;
	VectorXd regionSize = (hiBounds - loBounds);
	size_t c;
	srand(time(NULL));
	for (c = 0; c < maxContractions; c++)
	{
		samples.setRandom();
		// Get in the range 0 to 1
		samples.array() += 1;
		samples.array() /= 2.;
		//std::cout << "Contraction : " << c << std::endl;		
		for (int s = 0; s < nS; s++)
		{
			samples.col(s).array() *= regionSize.array();
			samples.col(s) += loBounds;
			VectorXd diffs(f.values());
			f(samples.col(s), diffs);
			sampleRes(s) = diffs.squaredNorm();
		}
		indices = arg_partial_sort(sampleRes, nR);
		//std::cout << std::endl;
		//std::cout << "First 10 residuals: " << sampleRes.head(10).transpose() << std::endl;
		//std::cout << "Retained residuals: ";
		for (int i = 0; i < nR; i++)
		{
			retained.col(i) = samples.col(indices[i]);
			//std::cout << indices[i] << " " << sampleRes[indices[i]] << " ";
		}
		//std::cout << std::endl;*/
			
		// Find the min and max for each parameter in the top nR samples
		loBounds = retained.rowwise().minCoeff();
		hiBounds = retained.rowwise().maxCoeff();
		regionSize = (hiBounds - loBounds);	
		// Terminate if ALL the distances between bounds are under the threshold
		if ((((regionSize.array() / hiBounds.array()).abs() - thresh) < 0).all())
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
		//std::cout << "Lo:   " << loBounds.transpose() << std::endl;
		//std::cout << "Hi:   " << hiBounds.transpose() << std::endl;
		//std::cout << "Size: " << regionSize.transpose() << std::endl;
	}
	// Return the best evaluated solution so far
	params = retained.col(0);
	//std::cout << "End P: " << params.transpose() << " Res: " << sampleRes[indices[0]] << std::endl;
	return sampleRes[indices[0]];
}

#endif
