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
DESPOTFunctor::DESPOTFunctor(shared_ptr<Model> m, const ArrayXd &data, const double B1, const bool debug) :
	m_model(m), m_data(data), m_B1(B1), m_debug(debug)
{
	m_nV = m_model->size();
	assert(m_data.rows() == m_nV);
}

int DESPOTFunctor::operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) const {
	eigen_assert(diffs.size() == values());
	ArrayXd s = m_model->signal(params, m_B1);
	diffs = s - m_data;
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
	return m_model->validParameters(params);
}
