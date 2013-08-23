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
class RegionContraction {
	public:
		enum class Status {
			NotStarted = -1,
			DidNotConverge, Converged, ErrorConstraints, ErrorInfiniteResidual
		};
	
	private:
		Functor_t &m_f;
		ArrayXXd m_startBounds, m_currentBounds;
		ArrayXd m_weights, m_residuals;
		size_t m_nS, m_nR, m_maxContractions, m_contractions;
		double m_thresh, m_expand;
		Status m_status;
		bool m_debug;
	
	public:
	
		RegionContraction(Functor_t &f, const Ref<ArrayXXd> &startBounds, const Ref<ArrayXd> &weights,
						  const int nS = 5000, const int nR = 50, const int maxContractions = 10,
						  const double thresh = 0.05, const double expand = 0., const bool debug = false) :
				m_f(f), m_startBounds(startBounds), m_currentBounds(startBounds),
				m_nS(nS), m_nR(nR), m_maxContractions(maxContractions),
				m_thresh(thresh), m_expand(expand), m_residuals(f.values()), m_contractions(0),
				m_status(Status::NotStarted), m_weights(weights), m_debug(debug)
		{
			eigen_assert(f.inputs() == startBounds.rows());
			eigen_assert(startBounds.cols() == 2);
			eigen_assert(weights.rows() == f.values());
		}
		
		const ArrayXXd &startBounds() const { return m_startBounds; }
		void setBounds(const Ref<ArrayXXd> &b) {
			eigen_assert(m_f.inputs() == b.rows());
			eigen_assert(b.cols() == 2);
			m_startBounds = b;
		}
		const ArrayXd &weights() const { return m_weights; }
		void setWeights(const Ref<ArrayXd> &w) {
			eigen_assert(w.rows() == m_f.values());
			m_weights = w;
		}
		const ArrayXd &residuals() const { return m_residuals; }
		const size_t contractions() const { return m_contractions; }
		const Status status() const { return m_status; }
		const ArrayXXd &currentBounds() const { return m_currentBounds; }
		const ArrayXd width() const { return m_currentBounds.col(1) - m_currentBounds.col(0); }
		const ArrayXd midPoint() const { return (m_currentBounds.rowwise().sum() / 2.); }
		
		void optimise(Ref<ArrayXd> params, const int seed = 0) {
			static atomic<bool> finiteWarning(false);
			static atomic<bool> constraintWarning(false);
			eigen_assert(m_f.inputs() == params.size());
			int nP = static_cast<int>(params.size());
			ArrayXXd samples(m_f.inputs(), m_nS);
			ArrayXXd retained(m_f.inputs(), m_nR);
			ArrayXXd residuals(m_f.values(), m_nS);
			vector<size_t> indices(m_nR);
			m_currentBounds = m_startBounds;
			m_residuals.setZero();
			
			if (m_debug) {
				cout << "Starting bounds: " << endl << m_startBounds.transpose() << endl;
				cout << "Mid-point: " << midPoint().transpose() << endl;
				cout << "Width:     " << width().transpose() << endl;
			}
			
			mt19937 twist(seed);
			uniform_real_distribution<double> uniform(0., 1.);
			
			m_status = Status::DidNotConverge;
			for (m_contractions = 0; m_contractions < m_maxContractions; m_contractions++) {
				for (int s = 0; s < m_nS; s++) {
					ArrayXd tempSample(nP);
					size_t nTries = 0;
					do {
						for (int p = 0; p < nP; p++)
							tempSample(p) = uniform(twist);
						tempSample *= width();
						tempSample += m_currentBounds.col(0);
						nTries++;
						if (nTries > 100) {
							if (!constraintWarning) {
								constraintWarning = true;
								cout << "Warning: Cannot fulfill sample constraints after " << to_string(nTries) << " attempts, giving up." << endl;
								cout << "Last attempt was: " << tempSample.transpose() << endl;
								cout << "This warning will only be printed once." << endl;
							}
							params.setZero();
							m_status = Status::ErrorConstraints;
							return;
						}
					} while (!m_f.constraint(tempSample));
					m_f(tempSample, residuals.col(s));
					if (!isfinite(m_residuals.square().sum())) {
						if (!finiteWarning) {
							finiteWarning = true;
							cout << "Warning: Non-finite residual found!" << endl
								 << "Result may be meaningless. This warning will only be printed once." << endl;
							cout << "Parameters were " << tempSample.transpose() << endl;
							cout << "Signal " << m_f.signals().transpose() << endl;
							cout << "Theory " << m_f.theory(tempSample).transpose() << endl;
						}
						params = retained.col(0);
						m_residuals.setConstant(numeric_limits<double>::infinity());
						m_status = Status::ErrorInfiniteResidual;
						return;
					}
					residuals.col(s) *= m_weights;
					samples.col(s) = tempSample;
				}
				ArrayXd toSort = residuals.square().colwise().sum();
				indices = index_partial_sort(toSort, m_nR);
				for (int i = 0; i < m_nR; i++) {
					retained.col(i) = samples.col(indices[i]);
					if (m_debug) {
						cout << "Sample " << indices[i] << " RES " << residuals.col(indices[i]).square().sum() << "\t PARAMS " << retained.col(i).transpose();
						cout << "\t " << residuals.col(indices[i]).transpose() << endl;
					}
				}
				// Find the min and max for each parameter in the top nR samples
				if (m_debug) {
					cout << "RES Min: " << toSort.minCoeff() << " Max: " << toSort.maxCoeff() << endl;
					cout << "Before search MID: " << midPoint().transpose() << " WIDTH: " << width().transpose() << endl;
				}
				m_currentBounds.col(0) = retained.rowwise().minCoeff();
				m_currentBounds.col(1) = retained.rowwise().maxCoeff();
				// Terminate if ALL the distances between bounds are under the threshold
				if (((width() / midPoint()).abs() < m_thresh).all()) {
					m_status = Status::Converged;
					break;
				}
				
				// Expand the boundaries back out in case we just missed a minima,
				// but don't go past initial boundaries
				if (m_debug)
					cout << "After search MID: " << midPoint().transpose() << " WIDTH: " << width().transpose() << endl;
				ArrayXd tempW = width(); // Because altering .col(0) will change width
				m_currentBounds.col(0) = (m_currentBounds.col(0) - tempW * m_expand).max(m_startBounds.col(0));
				m_currentBounds.col(1) = (m_currentBounds.col(1) + tempW * m_expand).min(m_startBounds.col(1));
				if (m_debug)
					cout << "After expand MID: " << midPoint().transpose() << " WIDTH: " << width().transpose() << endl;
			}
			// Return the best evaluated solution so far
			params = retained.col(0);
			// Calculate the residuals
			m_f(params, m_residuals);
			//diffs /= f.signals();
			m_contractions++; // Adjust for 0-based looping
			if (m_debug) {
				cout << "Finished, contractions = " << m_contractions << endl;
				cout << "Mid-point: " << midPoint().transpose() << endl;
				cout << "Width:     " << width().transpose() << endl;
			}
		}
};
#endif
