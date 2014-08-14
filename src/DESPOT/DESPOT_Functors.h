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

#include "Nifti/Nifti.h"
#include "DESPOT.h"
#include "Model.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
#pragma mark Optimisation Functor Base Class
//******************************************************************************
// From Nonlinear Tests in Eigen 
template<typename _Scalar, int NX=Dynamic, int NY=Dynamic>
class Functor {
	public:
		typedef _Scalar Scalar;
		enum {
			InputsAtCompileTime = NX,
			ValuesAtCompileTime = NY
		};
		typedef Matrix<Scalar,InputsAtCompileTime,1> InputType;
		typedef Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
		typedef Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;
		
		const long m_inputs, m_values;
		
		Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
		Functor(long inputs, long values) : m_inputs(inputs), m_values(values) {}
		Functor(long values) : m_inputs(InputsAtCompileTime), m_values(values) {}
		
		virtual ~Functor() {};
		
		virtual const long inputs() const { return m_inputs; }
		virtual const long values() const { return m_values; }
		
		virtual int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const = 0;
};

//******************************************************************************
#pragma mark Basic DESPOT Functor
//******************************************************************************
class DESPOTFunctor : public Functor<double> {
	protected:
		Sequences &m_sequences;
		size_t m_nV;
		ArrayXcd m_data;
		
	public:
		bool m_complex, m_debug;
		double m_B1;
		
		const long inputs() const override { return m_sequences.nParameters(); }
		const long values() const override { return m_nV; }
		
		DESPOTFunctor(Sequences &cs, const ArrayXcd &d, const double B1, const bool fitComplex, const bool debug = false);
		
		const bool constraint(const VectorXd &params) const;
		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const override;
};

#endif
