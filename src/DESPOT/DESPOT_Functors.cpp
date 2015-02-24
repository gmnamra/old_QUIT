//
//  DESPOT_Functors.cpp
//  DESPOT
//
//  Created by Tobias Wood on 25/09/2013.
//
//

#include "DESPOT_Functors.h"
//******************************************************************************
#pragma mark DESPOTFunctor
//******************************************************************************
DESPOTFunctor::DESPOTFunctor(SequenceBase &cs, const Pools np, const ArrayXcd &data, const double B1, const bool fitComplex, const bool debug) :
	m_sequence(cs), m_p(np), m_data(data), m_B1(B1), m_complex(fitComplex), m_debug(debug)
{
	m_nV = m_sequence.size();
	assert(static_cast<size_t>(m_data.rows()) == m_nV);
}

int DESPOTFunctor::operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
	eigen_assert(diffs.size() == values());
	ArrayXcd s = m_sequence.signal(m_p, params, m_B1);
	if (m_complex) {
		diffs = (s - m_data).abs();
	} else {
		diffs = s.abs() - m_data.abs();
	}
	if (m_debug) {
		cout << endl << __PRETTY_FUNCTION__ << endl;
		cout << "p:     " << params.transpose() << endl;
		cout << "s:     " << s.transpose() << endl;
		cout << "data:  " << m_data.transpose() << endl;
		cout << "diffs: " << diffs.transpose() << endl;
	}
	return 0;
}

const bool DESPOTFunctor::constraint(const VectorXd &params) const {
	return PoolInfo::ValidParameters(m_p, params);
}
//******************************************************************************
// D2 Functor
//******************************************************************************
D2Functor::D2Functor(const double T1, SequenceBase &cs, const Pools np, const ArrayXcd &data, const double B1, const bool fitComplex, const bool debug) :
	DESPOTFunctor(cs,np,data,B1,fitComplex,debug),
	m_T1(T1)
{}

int D2Functor::operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
	Array4d fullparams;
	fullparams << params(0), m_T1, params(1), params(2);
	return DESPOTFunctor::operator()(fullparams, diffs);
}

const bool D2Functor::constraint(const VectorXd &params) const {
	Array4d fullparams;
	fullparams << params(0), m_T1, params(1), params(2);
	return PoolInfo::ValidParameters(m_p, fullparams);
}

// HIFI Functor
HIFIFunctor::HIFIFunctor(SequenceBase &cs, const ArrayXd &data, const bool debug) :
	DenseFunctor<double>(3, cs.size()),
	m_sequence(cs), m_data(data), m_debug(debug) {
	assert(static_cast<size_t>(m_data.rows()) == values());
}

int HIFIFunctor::operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
	eigen_assert(diffs.size() == values());
	ArrayXcd s = m_sequence.signal(Pools::One, params.head(2), params(2));
	diffs = s.abs() - m_data;
	if (m_debug) {
		cout << endl << __PRETTY_FUNCTION__ << endl;
		cout << "p:     " << params.transpose() << endl;
		cout << "s:     " << s.transpose() << endl;
		cout << "data:  " << m_data.transpose() << endl;
		cout << "diffs: " << diffs.transpose() << endl;
	}
	return 0;
}
