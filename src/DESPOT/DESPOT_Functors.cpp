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
DESPOTFunctor::DESPOTFunctor(Model &m, const ArrayXcd &data, const double B1, const bool fitComplex, const bool debug) :
	m_model(m), m_data(data), m_B1(B1), m_complex(fitComplex), m_debug(debug)
{
	m_nV = m_model.size();
	assert(static_cast<size_t>(m_data.rows()) == m_nV);
}

int DESPOTFunctor::operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
	eigen_assert(diffs.size() == values());
	ArrayXcd s = m_model.signal(params, m_B1);
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
	return m_model.validParameters(params);
}
