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
#include <mutex>

#include <Eigen/Dense>

#include "QUIT/Util.h"

using namespace std;
using namespace Eigen;

typedef Array<bool, Dynamic, 1> ArrayXb;

vector<size_t> index_partial_sort(const Ref<ArrayXd> &x, ArrayXd::Index N);
vector<size_t> index_partial_sort(const Ref<ArrayXd> &x, ArrayXd::Index N)
{
	eigen_assert(x.size() >= N);
    vector<size_t> allIndices(x.size()), indices(N);
    for(size_t i = 0; i < allIndices.size(); i++) {
		allIndices[i] = i;
    }
	partial_sort(allIndices.begin(), allIndices.begin() + N, allIndices.end(),
	             [&x](size_t i1, size_t i2) { return x[i1] < x[i2]; });
	for (ArrayXd::Index i = 0; i < N; i++) {
		indices[i] = allIndices[i];
	}
    return indices;
}

enum class RCStatus {
	NotStarted = -1,
	Converged, IterationLimit, ErrorInvalid, ErrorResidual
};

ostream& operator<<(ostream &os, const RCStatus &s) {
	switch (s) {
		case RCStatus::NotStarted: os << "Not Started"; break;
		case RCStatus::Converged: os << "Converged"; break;
		case RCStatus::IterationLimit: os << "Reached iteration limit"; break;
		case RCStatus::ErrorInvalid: os << "Could not generate valid sample"; break;
		case RCStatus::ErrorResidual: os << "Infinite residual found"; break;
	}
	return os;
}

template <typename Functor_t>
class RegionContraction {
	private:
		Functor_t &m_f;
		mt19937_64 m_rng;
		ArrayXXd m_startBounds, m_currentBounds;
		ArrayXd m_weights, m_residuals, m_threshes;
		size_t m_nS, m_nR, m_maxContractions, m_contractions;
		double m_expand, m_SoS;
		RCStatus m_status;
		bool m_gaussian, m_debug;
	
	public:
		RegionContraction(Functor_t &f,
						  const Ref<ArrayXXd> &startBounds, const ArrayXd &weights, const ArrayXd &thresh,
						  const int nS = 5000, const int nR = 50, const int maxContractions = 10,
						  const double expand = 0., const bool gauss = false, const bool debug = false, const int seed = -1) :
				m_f(f), m_startBounds(startBounds), m_currentBounds(startBounds),
				m_nS(nS), m_nR(nR), m_maxContractions(maxContractions),
				m_threshes(thresh), m_expand(expand), m_residuals(f.values()), m_contractions(0),
				m_status(RCStatus::NotStarted), m_weights(weights), m_gaussian(gauss), m_debug(debug)
		{
			eigen_assert(f.inputs() == startBounds.rows());
			eigen_assert(startBounds.cols() == 2);
			eigen_assert(weights.rows() == f.values());
			eigen_assert(thresh.rows() == f.inputs());
			eigen_assert((thresh >= 0.).all() && (thresh <= 1.).all());

			if (seed < 0) {
				m_rng = mt19937_64(QUIT::RandomSeed());
			} else {
				m_rng = mt19937_64(seed);
			}
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
		const ArrayXb &thresholds() const { return m_threshes; }
		void setThresholds(const Ref<ArrayXd> &t) {
			eigen_assert(t.rows() == m_f.inputs());
			eigen_assert((t >= 0.).all() && (t <= 1.).all());
			m_threshes = t;
		}
		const ArrayXd  &residuals() const { return m_residuals; }
		const size_t   contractions() const { return m_contractions; }
		RCStatus       status() const { return m_status; }
		const ArrayXXd &currentBounds() const { return m_currentBounds; }
		const double   SoS() const { return m_SoS; }
		const ArrayXd  startWidth() const { return m_startBounds.col(1) - m_startBounds.col(0); }
		const ArrayXd  width() const { return m_currentBounds.col(1) - m_currentBounds.col(0); }
		const ArrayXd  midPoint() const { return (m_currentBounds.rowwise().sum() / 2.); }
		
		void optimise(Ref<ArrayXd> params) {
			static atomic<bool> finiteWarning(false);
			static atomic<bool> constraintWarning(false);
			static atomic<bool> boundsWarning(false);
			mutex warn_mtx;

			eigen_assert(m_f.inputs() == params.size());
			int nP = static_cast<int>(params.size());
			ArrayXXd samples(m_f.inputs(), m_nS);
			ArrayXXd retained(m_f.inputs(), m_nR);
			ArrayXXd residuals(m_f.values(), m_nS);
			ArrayXXd retainedRes(m_f.values(), m_nR);
			ArrayXd gauss_mu(m_f.inputs()), gauss_sigma(m_f.inputs());
			vector<size_t> indices(m_nR);
			m_currentBounds = m_startBounds;
			m_residuals.setZero();

			if ((m_startBounds != m_startBounds).any() || (m_startBounds >= numeric_limits<double>::infinity()).any()) {
				warn_mtx.lock();
				if (!boundsWarning) {
					boundsWarning = true;
					cerr << "Warning: Starting boundaries do not make sense." << endl;
					cerr << "Bounds were: " << m_startBounds.transpose() << endl;
					cerr << "This warning will only be printed once." << endl;
				}
				warn_mtx.unlock();
				params.setZero();
				m_status = RCStatus::ErrorInvalid;

				return;
			}

			if (m_debug) {
				cout << endl;
				cout << "Start Boundaries: " << endl << m_startBounds.transpose() << endl;
				cout << "Weights:        " << m_weights.transpose() << endl;
			}
			
			uniform_real_distribution<double> uniform(0., 1.);
			m_status = RCStatus::IterationLimit;
			for (m_contractions = 0; m_contractions < m_maxContractions; m_contractions++) {
				size_t startSample = 0;
				if (m_contractions > 0) {
					// Keep the retained samples to prevent boundary contracting too fast
					for (size_t s = 0; s < m_nR; s++) {
						samples.col(s) = retained.col(s);
						residuals.col(s) = retainedRes.col(s);
					}
					startSample = m_nR;
				}
				for (size_t s = startSample; s < m_nS; s++) {
					ArrayXd tempSample(nP);
					size_t nTries = 0;
					do {
						if (!m_gaussian || (m_contractions == 0)) {
							for (int p = 0; p < nP; p++) {
								tempSample(p) = uniform(m_rng);
							}
							tempSample *= width();
							tempSample += m_currentBounds.col(0);
						} else {
							for (int p = 0; p < nP; p++) {
								normal_distribution<double> gauss(gauss_mu(p), gauss_sigma(p));
								tempSample(p) = gauss(m_rng);
								if (tempSample(p) < m_currentBounds(p, 0))
									tempSample(p) = m_currentBounds(p, 0);
								if (tempSample(p) > m_currentBounds(p, 1))
									tempSample(p) = m_currentBounds(p, 1);
							}
						}
						nTries++;
						if (nTries > 100) {
							warn_mtx.lock();
							if (!constraintWarning) {
								constraintWarning = true;
								cerr << "Warning: Cannot fulfill sample constraints after " << to_string(nTries) << " attempts, giving up." << endl;
								cerr << "Last attempt was: " << tempSample.transpose() << endl;
								cerr << "This warning will only be printed once." << endl;
							}
							warn_mtx.unlock();
							params.setZero();
							m_status = RCStatus::ErrorInvalid;
							return;
						}
					} while (!m_f.constraint(tempSample));
					
					m_f(tempSample, residuals.col(s));
					if (!isfinite(residuals.col(s).square().sum())) {
						warn_mtx.lock();
						if (!finiteWarning) {
							finiteWarning = true;
							cout << "Warning: Non-finite residual found!" << endl
								 << "Result may be meaningless. This warning will only be printed once." << endl;
							cout << "Parameters were " << tempSample.transpose() << endl;
						}
						warn_mtx.unlock();
						params = retained.col(0);
						m_residuals.setConstant(numeric_limits<double>::infinity());
						m_status = RCStatus::ErrorResidual;
						return;
					}
					samples.col(s) = tempSample;
				}
				ArrayXd toSort = (residuals.colwise() * m_weights).square().colwise().sum();
				indices = index_partial_sort(toSort, m_nR);
				ArrayXd previousBest = retained.col(0);
				for (size_t i = 0; i < m_nR; i++) {
					retained.col(i) = samples.col(indices[i]);
					retainedRes.col(i) = residuals.col(indices[i]);
				}
				// Find the min and max for each parameter in the top nR samples
				m_currentBounds.col(0) = retained.rowwise().minCoeff();
				m_currentBounds.col(1) = retained.rowwise().maxCoeff();
				if (m_gaussian) {
					gauss_mu = retained.rowwise().mean();
					gauss_sigma = ((retained.colwise() - gauss_mu).square().rowwise().sum() / m_f.inputs()).sqrt();
				}
				// Terminate if all the desired parameters have converged
				if (m_debug) {
					cout << "Best sample:    " << retained.col(0).transpose() << endl;
					cout << "Best residual:  " << retainedRes.col(0).transpose() << endl;
					cout << "Best SoS:       " << retainedRes.col(0).square().sum() << endl;
					cout << "Current bounds: " << endl << m_currentBounds.transpose() << endl;
					cout << "Current width:  " << width().transpose() << endl;
					cout << "Current thresh: " << (m_threshes * startWidth()).transpose() << endl;
					cout << "Width < Thresh: " << (width() <= (m_threshes * startWidth())).transpose() << endl;
					cout << "Converged:      " << (width() <= (m_threshes * startWidth())).all() << endl;
					if (m_gaussian) {
						cout << "Gaussian mu:    " << gauss_mu.transpose() << endl;
						cout << "Gaussian sigma: "<<  gauss_sigma.transpose() << endl;
					}
				}
				if (((width() <= (m_threshes * startWidth())).all()) ||
				    (previousBest == retained.col(0)).all()) {
					m_status = RCStatus::Converged;
					m_contractions++; // Just to give an accurate contraction count.
					break;
				}
				
				// Expand the boundaries back out in case we just missed a minima,
				// but don't go past initial boundaries
				ArrayXd tempW = width(); // Because altering .col(0) will change width
				m_currentBounds.col(0) = (m_currentBounds.col(0) - tempW * m_expand).max(m_startBounds.col(0));
				m_currentBounds.col(1) = (m_currentBounds.col(1) + tempW * m_expand).min(m_startBounds.col(1));
			}
			// Return the best evaluated solution so far
			params = retained.col(0);
			m_residuals = retainedRes.col(0);
			m_SoS = m_residuals.square().sum();
			if (m_debug) {
				cout << "Finished, contractions = " << m_contractions << endl;
			}
		}
};

#endif
