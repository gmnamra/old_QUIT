/*
 *  DESPOT_Functors.h
 *
 *  Created by Tobias Wood on 16/08/2012.
 *  Copyright (c) 2012-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef DESPOT_Functors_h
#define DESPOT_Functors_h

#include <vector>
#include <array>
#include <map>
#include <iostream>
#include <exception>
#include <string>
#include <memory>
#include <Eigen/Dense>
#include <unsupported/Eigen/LevenbergMarquardt>
#include "unsupported/Eigen/NumericalDiff"

#include "Nifti/Nifti.h"
#include "DESPOT.h"
#include "Model.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
#pragma mark Basic DESPOT Functor
//******************************************************************************
class DESPOTFunctor : public DenseFunctor<double> {
	protected:
		SequenceBase &m_sequence;
		Pools m_p;
		size_t m_nV;
		ArrayXcd m_data;
		
	public:
		bool m_complex, m_debug;
		double m_B1;
		
		const long inputs() const { return PoolInfo::nParameters(m_p); }
		const long values() const { return m_nV; }
		
		DESPOTFunctor(SequenceBase &s, const Pools p, const ArrayXcd &d, const double B1, const bool fitComplex, const bool debug = false);
		
		const bool constraint(const VectorXd &params) const;
		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const;
};

//******************************************************************************
// T2 Only Functor
//******************************************************************************
class D2Functor : public DESPOTFunctor {
	public:
		double m_T1;
		D2Functor(const double T1, SequenceBase &s, const Pools p, const ArrayXcd &d, const double B1, const bool fitComplex, const bool debug = false);
		const long inputs() const { return 3; }
		const bool constraint(const VectorXd &params) const;
		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const;
};

// HIFI Functor - includes optimising B1
class HIFIFunctor : public DenseFunctor<double> {
	protected:
		const SequenceBase &m_sequence;
		const ArrayXd m_data;
		const bool m_debug;

	public:
		HIFIFunctor(SequenceBase &s, const ArrayXd &d, const bool debug = false);
		const bool constraint(const VectorXd &params) const;
		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const;
};

#endif
