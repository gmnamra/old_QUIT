/*
 *  Sequence.cpp
 *
 *  Created by Tobias Wood on 14/11/2012.
 *  Copyright (c) 2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Sequence.h"

const string to_string(const Scale &p) {
	static const string sn{"None"}, snm{"Normalised to Mean"};
	switch (p) {
		case Scale::None: return sn;
		case Scale::NormToMean: return snm;
	}
}

//******************************************************************************
#pragma mark Sequence Functors
//******************************************************************************
Sequence::Sequence() : SequenceBase() {}
Sequence::Sequence(const ArrayXd &flip, const double TR) :
	SequenceBase(), m_flip(flip), m_TR(TR)
{}

ostream& operator<<(ostream& os, const SequenceBase& s) {
	s.write(os);
	return os;
}

ArrayXd Sequence::B1flip(const double B1) const {
	return B1 * m_flip;
}

// Parameters are PD, T1, T2, f0
//                PD, T1_a, T2_a, T1_b, T2_b, tau_a, f_a, f0
//				  PD, T1_a, T2_a, T1_b, T2_b, T1_c, T2_c, tau_a, f_a, f_c, f0
SPGRSimple::SPGRSimple(const ArrayXd &flip, const double TR) :
	Sequence(flip, TR)
{}
SPGRSimple::SPGRSimple(const bool prompt, const Agilent::ProcPar &pp) :
	Sequence()
{
	if (pp) {
		m_flip = pp.realValues("flip1") * M_PI / 180.;
		m_TR = pp.realValue("tr");
	} else {
		size_t nFlip;
		if (prompt) cout << "Enter number of flip-angles: " << flush;
		QUIT::Read<size_t>::FromLine(cin, nFlip);
		ArrayXd inAngles(nFlip);
		if (prompt) cout << "Enter " << inAngles.size() << " flip-angles (degrees): " << flush;
		QUIT::ReadEigenFromLine(cin, inAngles);
		m_flip = inAngles * M_PI / 180.;
		if (prompt) cout << "Enter TR (seconds): " << flush;
		QUIT::Read<double>::FromLine(cin, m_TR);
	}
}

void SPGRSimple::write(ostream &os) const {
	os << "SPGR Simple" << endl;
	os << "TR: " << m_TR << endl;
	os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

ArrayXcd SPGRSimple::signal(shared_ptr<Model> m, const VectorXd &p, const double B1) const {
	return m->SPGR(p, B1flip(B1), m_TR);
}

SPGRFinite::SPGRFinite(const ArrayXd &flip, const double TR, const double Trf, const double TE) :
	SPGRSimple(flip, TR), m_Trf(Trf), m_TE(TE)
{}
SPGRFinite::SPGRFinite(const bool prompt, const Agilent::ProcPar &pp) :
	SPGRSimple(prompt, pp)
{
	if (pp) {
		m_Trf = pp.realValue("p1") / 1.e6; // p1 is in microseconds
		m_TE = pp.realValue("te");
	} else {
		if (prompt) cout << "Enter RF Pulse Length (seconds): " << flush;
		QUIT::Read<double>::FromLine(cin, m_Trf);
		if (prompt) cout << "Enter TE (seconds): " << flush;
		QUIT::Read<double>::FromLine(cin, m_TE);
	}
}

void SPGRFinite::write(ostream &os) const {
	os << "SPGR Finite" << endl;
	os << "TR: " << m_TR << "\tTrf: " << m_Trf << "\tTE: " << m_TE << endl;
	os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

ArrayXcd SPGRFinite::signal(shared_ptr<Model> m, const VectorXd &p, const double B1) const {
	return m->SPGRFinite(p, B1flip(B1), m_TR, m_Trf, m_TE);
}

MPRAGE::MPRAGE(const ArrayXd &TI, const double TD, const double TR, const int N, const double flip) :
	Sequence(), m_TI(TI), m_TD(TD), m_N(N) {
	m_TR = TR;
	m_flip.resize(1); m_flip[0] = flip;
}

MPRAGE::MPRAGE(const bool prompt, const Agilent::ProcPar &pp) : Sequence() {
	if (pp) {
		throw(runtime_error("MPRAGE Procpar reader not implemented."));
	} else {
		if (prompt) cout << "Enter read-out flip-angle (degrees): " << flush;
		ArrayXd inFlip(1);
		QUIT::ReadEigenFromLine(cin, inFlip);
		m_flip = inFlip * M_PI / 180.;
		if (prompt) cout << "Enter read-out TR (seconds): " << flush;
		QUIT::Read<double>::FromLine(cin, m_TR);
		if (prompt) cout << "Enter segment size: " << flush;
		QUIT::Read<int>::FromLine(cin, m_N);
		size_t nTI;
		if (prompt) cout << "Enter number of inversion times: " << flush;
		QUIT::Read<size_t>::FromLine(cin, nTI);
		m_TI.resize(nTI);
		if (prompt) cout << "Enter " << m_TI.size() << " inversion times (seconds): " << flush;
		QUIT::ReadEigenFromLine(cin, m_TI);
		if (prompt) cout << "Enter delay time (seconds): " << flush;
		QUIT::Read<double>::FromLine(cin, m_TD);
	}
}

IRSPGR::IRSPGR(const bool prompt, const Agilent::ProcPar &pp) : MPRAGE() {
	if (pp)
		throw(runtime_error("IRSPGR Procpar reader not implemented."));
	if (prompt) cout << "Enter read-out flip-angle (degrees): " << flush;
	ArrayXd inFlip(1);
	QUIT::ReadEigenFromLine(cin, inFlip);
	m_flip = inFlip * M_PI / 180.;
	if (prompt) cout << "Enter read-out TR (seconds): " << flush;
	QUIT::Read<double>::FromLine(cin, m_TR);

	int NPE2;
	if (prompt) cout << "Enter original number of slices (PE2):";
	QUIT::Read<int>::FromLine(cin, NPE2);
	m_N = (NPE2 / 2) + 2;

	int nTI;
	if (prompt) cout << "Enter number of TIs: " << flush;
	QUIT::Read<int>::FromLine(cin, nTI);
	m_TI.resize(nTI);
	if (prompt) cout << "Enter " << m_TI.size() << " TIs (seconds): " << flush;
	QUIT::ReadEigenFromLine(cin, m_TI);
}

ArrayXcd MPRAGE::signal(shared_ptr<Model> m, const VectorXd &par, const double B1) const {
	return m->MPRAGE(par, m_flip[0] * B1, m_TR, m_N, m_TI, m_TD);
}

void MPRAGE::write(ostream &os) const {
	os << name() << endl;
	os << "TR: " << m_TR << "\tN: " << m_N << "\tAlpha: " << m_flip[0] * 180 / M_PI << "\tTD: " << m_TD << endl;
	os << "TI: " << m_TI.transpose() << endl;
	os << "TS: " << (m_TI + m_N*m_TR + m_TD).transpose() << endl;
}

SSFPSimple::SSFPSimple(const ArrayXd &flip, const double TR, const ArrayXd &phases) :
	Sequence(flip, TR), m_phases(phases)
{}
SSFPSimple::SSFPSimple(const bool prompt, const Agilent::ProcPar &pp) :
	Sequence()
{
	if (pp) {
		m_phases = pp.realValues("rfphase") * M_PI / 180.;
		m_flip = pp.realValues("flip1") * M_PI / 180.;
		m_TR = pp.realValue("tr");
	} else {
		size_t nFlip;
		if (prompt) cout << "Enter number of flip-angles: " << flush;
		QUIT::Read<size_t>::FromLine(cin, nFlip);
		ArrayXd inAngles(nFlip);
		if (prompt) cout << "Enter " << inAngles.size() << " flip-angles (degrees): " << flush;
		QUIT::ReadEigenFromLine(cin, inAngles);
		m_flip = inAngles * M_PI / 180.;
		size_t nPhases;
		if (prompt) cout << "Enter number of phase-cycles: " << flush;
		QUIT::Read<size_t>::FromLine(cin, nPhases);
		ArrayXd inPhases(nPhases);
		if (prompt) cout << "Enter " << inPhases.size() << " phase-cycles (degrees): " << flush;
		QUIT::ReadEigenFromLine(cin, inPhases);
		m_phases = inPhases * M_PI / 180.;
		if (prompt) cout << "Enter TR (seconds): " << flush;
		QUIT::Read<double>::FromLine(cin, m_TR);
	}
}

void SSFPSimple::write(ostream &os) const {
	os << "SSFP Simple" << endl;
	os << "TR: " << m_TR << "\tPhases: " << (m_phases * 180. / M_PI).transpose() << endl;
	os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

size_t SSFPSimple::phases() const { return m_phases.rows(); }

ArrayXcd SSFPSimple::signal(shared_ptr<Model> m, const VectorXd &p, const double B1) const {
	ArrayXcd s(size());
	ArrayXcd::Index start = 0;
	for (ArrayXcd::Index i = 0; i < m_phases.rows(); i++) {
		s.segment(start, m_flip.rows()) = m->SSFP(p, B1flip(B1), m_TR, m_phases(i));
		start += m_flip.rows();
	}
	return s;
}

SSFPFinite::SSFPFinite(const ArrayXd &flip, const double TR, const double Trf, const ArrayXd &phases) :
	SSFPSimple(flip, TR, phases), m_Trf(Trf)
{}
SSFPFinite::SSFPFinite(const bool prompt, const Agilent::ProcPar &pp) :
	SSFPSimple(prompt, pp)
{
	if (pp) {
		m_Trf = pp.realValue("p1") / 1.e6; // p1 is in microseconds
	} else {
		if (prompt) cout << "Enter RF Pulse Length (seconds): " << flush;
		QUIT::Read<double>::FromLine(cin, m_Trf);
	}
}

void SSFPFinite::write(ostream &os) const {
	os << "SSFP Finite" << endl;
	os << "TR: " << m_TR << "\tTrf: " << m_Trf << "\tPhases: " << (m_phases * 180. / M_PI).transpose() << endl;
	os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

ArrayXcd SSFPFinite::signal(shared_ptr<Model> m, const VectorXd &p, const double B1) const {
	ArrayXcd s(size());
	ArrayXcd::Index start = 0;
	for (ArrayXcd::Index i = 0; i < m_phases.rows(); i++) {
		s.segment(start, m_flip.rows()) = m->SSFPFinite(p, B1flip(B1), m_TR, m_Trf, m_phases(i));
		start += m_flip.rows();
	}
	return s;
}

SSFPEllipse::SSFPEllipse(const bool prompt, const Agilent::ProcPar &pp) :
	Sequence()
{
	if (pp) {
		m_flip = pp.realValues("flip1") * M_PI / 180.;
		m_TR = pp.realValue("tr");
	} else {
		size_t nFlip;
		if (prompt) cout << "Enter number of flip-angles: " << flush;
		QUIT::Read<size_t>::FromLine(cin, nFlip);
		ArrayXd inAngles(nFlip);
		if (prompt) cout << "Enter " << inAngles.size() << " flip-angles (degrees): " << flush;
		QUIT::ReadEigenFromLine(cin, inAngles);
		m_flip = inAngles * M_PI / 180.;
		if (prompt) cout << "Enter TR (seconds): " << flush;
		QUIT::Read<double>::FromLine(cin, m_TR);
	}
}

void SSFPEllipse::write(ostream &os) const {
	os << "SSFP Ellipse" << endl;
	os << "TR: " << m_TR << endl;
	os << "Angles: " << (m_flip * 180. / M_PI).transpose() << endl;
}

ArrayXcd SSFPEllipse::signal(shared_ptr<Model> m, const VectorXd &p, const double B1) const {
	return m->SSFPEllipse(p, B1*m_flip, m_TR);
}

//******************************************************************************
#pragma mark Sequences Class
//******************************************************************************
SequenceGroup::SequenceGroup(const Scale s) :
	SequenceBase(), m_scaling(s)
{}

void SequenceGroup::write(ostream &os) const {
	os << "Combined Sequence Count: " << m_sequences.size() << "\tCombined size: " << size() << endl;
	for (auto& s : m_sequences)
		os << *s;
}

size_t SequenceGroup::count() const {
	return m_sequences.size();
}

shared_ptr<Sequence> SequenceGroup::sequence(const size_t i) const {
	return m_sequences.at(i);
}

vector<shared_ptr<Sequence>> &SequenceGroup::sequences() {
	return m_sequences;
}

size_t SequenceGroup::size() const {
	size_t sz = 0;
	for (auto& sig : m_sequences)
		sz += sig->size();
	return sz;
}

ArrayXcd SequenceGroup::signal(shared_ptr<Model> m, const VectorXd &p, const double B1) const {
	ArrayXcd result(size());
	size_t start = 0;
	for (auto &sig : m_sequences) {
		ArrayXcd thisResult = sig->signal(m, p, B1);
		switch (m_scaling) {
			case Scale::None :       break;
			case Scale::NormToMean : thisResult /= thisResult.abs().mean();
		}
		result.segment(start, sig->size()) = thisResult;
		start += sig->size();
	}
	return result;
}

double SequenceGroup::minTR() const {
	double minTR = numeric_limits<double>::max();
	for (auto &s : m_sequences) {
		if (s->m_TR < minTR)
			minTR = s->m_TR;
	}
	return minTR;
}

ArrayXcd SequenceGroup::loadSignals(vector<QUIT::MultiArray<complex<float>, 4>> &sigs,
                                const size_t i, const size_t j, const size_t k,
                                const bool flip) const {
	ArrayXcd signal(size());
	size_t start = 0;
	for (size_t s = 0; s < m_sequences.size(); s++) {
		ArrayXcd thisSig = sigs.at(s).slice<1>({i,j,k,0},{0,0,0,-1}).asArray().cast<complex<double>>();
		if (flip) {
			ArrayXXcd flipped = Map<ArrayXXcd>(thisSig.data(), m_sequences.at(s)->phases(), m_sequences.at(s)->angles()).transpose();
			thisSig = Map<ArrayXcd>(flipped.data(), thisSig.rows(), 1);
		}
		if (m_scaling == Scale::NormToMean)
			thisSig /= thisSig.abs().mean();
		signal.segment(start, thisSig.rows()) = thisSig;
		start += thisSig.rows();
	}
	return signal;
}

void SequenceGroup::addSequence(const shared_ptr<Sequence> &seq) {
	m_sequences.push_back(seq);
}
