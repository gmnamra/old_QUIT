//
//  DESPOT_Functors.cpp
//  DESPOT
//
//  Created by Tobias Wood on 25/09/2013.
//
//

#include "DESPOT_Functors.h"

//******************************************************************************
#pragma mark Signal Functors
//******************************************************************************
Signal::Signal(const Components nC, const ArrayXd &flip, const double TR, const double weight) :
	m_flip(flip), m_TR(TR), m_weight(weight), m_nC(nC)
{}

SPGR::SPGR(const Components nC, const ArrayXd &flip, const double TR, const double weight) :
	Signal(nC, flip, TR, weight)
{}
ArrayXd SPGR::signal(const VectorXd &p, const double B1, const double f0) const {
	switch (m_nC) {
		case (Components::One) : return SigMag(One_SPGR(p, m_flip, m_TR, B1));
		case (Components::Two) : return SigMag(Two_SPGR(p, m_flip, m_TR, B1));
		case (Components::Three) : return SigMag(Three_SPGR(p, m_flip, m_TR, B1));
	}
}

SPGREcho::SPGREcho(const Components nC, const ArrayXd &flip, const double TR, const double TE, const double weight) :
	Signal(nC, flip, TR, weight), m_TE(TE)
{}
ArrayXd SPGREcho::signal(const VectorXd &p, const double B1, const double f0) const {
	switch (m_nC) {
		case (Components::One) : return SigMag(One_SPGR_Echo(p, m_flip, m_TR, m_TE, B1));
		case (Components::Two) : return SigMag(Two_SPGR_Echo(p, m_flip, m_TR, m_TE, B1));
		case (Components::Three) : return SigMag(Three_SPGR_Echo(p, m_flip, m_TR, m_TE, B1));
	}
}

SPGRFinite::SPGRFinite(const Components nC, const ArrayXd &flip, const double TR, const double Trf, const double TE, const double weight) :
	Signal(nC, flip, TR, weight), m_Trf(Trf), m_TE(TE)
{}
ArrayXd SPGRFinite::signal(const VectorXd &p, const double B1, double f0) const {
	switch (m_nC) {
		case (Components::One) : return SigMag(One_SSFP_Finite(p, m_flip, true, m_TR, m_Trf, m_TE, 0, B1, f0));
		case (Components::Two) : return SigMag(Two_SSFP_Finite(p, m_flip, true, m_TR, m_Trf, m_TE, 0, B1, f0));
		case (Components::Three) : return SigMag(Three_SSFP_Finite(p, m_flip, true, m_TR, m_Trf, m_TE, 0, B1, f0));
	}
}


SSFP::SSFP(const Components nC, const ArrayXd &flip, const double TR, const ArrayXd &phases, const double weight) :
	Signal(nC, flip, TR, weight), m_phases(phases)
{}
size_t SSFP::size() const {
	return m_flip.rows() * m_phases.rows();
}
ArrayXd SSFP::signal(const VectorXd &p, const double B1, double f0) const {
	ArrayXd s(size());
	ArrayXd::Index start = 0;
	for (ArrayXd::Index i = 0; i < m_phases.rows(); i++) {
		switch (m_nC) {
			case (Components::One) : s.segment(start, m_flip.rows()) = SigMag(One_SSFP(p, m_flip, m_TR, m_phases(i), B1, f0)); break;
			case (Components::Two) : s.segment(start, m_flip.rows()) = SigMag(Two_SSFP(p, m_flip, m_TR, m_phases(i), B1, f0)); break;
			case (Components::Three) : s.segment(start, m_flip.rows()) = SigMag(Three_SSFP(p, m_flip, m_TR, m_phases(i), B1, f0)); break;
		}
		start += m_flip.rows();
	}
	return s;
}

SSFPEcho::SSFPEcho(const Components nC, const ArrayXd &flip, const double TR, const ArrayXd &phases, const double weight) :
	SSFP(nC, flip, TR, phases, weight)
{}
ArrayXd SSFPEcho::signal(const VectorXd &p, const double B1, const double f0) const {
	ArrayXd s(size());
	ArrayXd::Index start = 0;
	for (ArrayXd::Index i = 0; i < m_phases.rows(); i++) {
		switch (m_nC) {
			case (Components::One) : s.segment(start, m_flip.rows()) = SigMag(One_SSFP_Echo(p, m_flip, m_TR, m_phases(i), B1, f0)); break;
			case (Components::Two) : s.segment(start, m_flip.rows()) = SigMag(Two_SSFP_Echo(p, m_flip, m_TR, m_phases(i), B1, f0)); break;
			case (Components::Three) : s.segment(start, m_flip.rows()) = SigMag(Three_SSFP_Echo(p, m_flip, m_TR, m_phases(i), B1, f0)); break;
		}
		start += m_flip.rows();
	}
	return s;
}	

SSFPFinite::SSFPFinite(const Components nC, const ArrayXd &flip, const double TR, const double Trf, const ArrayXd &phases, const double weight) :
	SSFP(nC, flip, TR, phases, weight), m_Trf(Trf)
{}
ArrayXd SSFPFinite::signal(const VectorXd &p, const double B1, double f0) const {
	ArrayXd s(size());
	ArrayXd::Index start = 0;
	for (ArrayXd::Index i = 0; i < m_phases.rows(); i++) {
		switch (m_nC) {
			case (Components::One) : s.segment(start, m_flip.rows()) = SigMag(One_SSFP_Finite(p, m_flip, false, m_TR, m_Trf, 0., m_phases(i), B1, f0)); break;
			case (Components::Two) : s.segment(start, m_flip.rows()) = SigMag(Two_SSFP_Finite(p, m_flip, false, m_TR, m_Trf, 0., m_phases(i), B1, f0)); break;
			case (Components::Three) : s.segment(start, m_flip.rows()) = SigMag(Three_SSFP_Finite(p, m_flip, false, m_TR, m_Trf, 0., m_phases(i), B1, f0)); break;
		}
		start += m_flip.rows();
	}
	return s;
}

//******************************************************************************
#pragma mark Parsing Functions
//******************************************************************************
shared_ptr<Signal> parseSPGR(const Components nC, const Model mdl, const size_t nFlip,
									const bool prompt, const bool use_weights) {
	double inTR = 0., inTrf = 0., inTE = 0., inWeight = 1.;
	ArrayXd inAngles(nFlip);
	if (prompt) cout << "Enter TR (seconds): " << flush; cin >> inTR;
	if (mdl == Model::Finite) {
		if (prompt) cout << "Enter RF Pulse Length (seconds): " << flush; cin >> inTrf;
	}
	if (mdl == Model::Echo || mdl == Model::Finite) {
		if (prompt) cout << "Enter TE (seconds): " << flush; cin >> inTE;
	}
	if (prompt) cout << "Enter " << inAngles.size() << " Flip-angles (degrees): " << flush;
	for (int i = 0; i < inAngles.size(); i++) cin >> inAngles[i];
	string temp; getline(cin, temp); // Just to eat the newline
	if (use_weights) {
		if (prompt) cout << "Enter weighting: " << flush; cin >> inWeight;
		string temp; getline(cin, temp); // Just to eat the newline
	}
	shared_ptr<Signal> f;
	switch (mdl) {
		case Model::Simple: f = make_shared<SPGR>(nC, inAngles * M_PI / 180., inTR, inWeight); break;
		case Model::Echo:   f = make_shared<SPGREcho>(nC, inAngles * M_PI / 180, inTR, inTE, inWeight); break;
		case Model::Finite: f = make_shared<SPGRFinite>(nC, inAngles * M_PI / 180, inTR, inTrf, inTE, inWeight); break;
	}
	return f;
}

shared_ptr<Signal> parseSSFP(const Components nC, const Model mdl, const size_t nVols,
                                    const bool prompt, const bool use_weights) {
	double inTR = 0., inTrf = 0., inWeight = 1.;
	ArrayXd inPhases, inAngles;
	size_t nPhases = 0;
	if (prompt) cout << "Enter number of phase-cycling patterns: " << flush; cin >> nPhases;
	inPhases.resize(nPhases);
	inAngles.resize(nVols / nPhases);
	if (prompt) cout << "Enter " << nPhases << " phase-cycles (degrees): " << flush;
	for (size_t i = 0; i < nPhases; i++) cin >> inPhases(i);
	if (prompt) cout << "Enter TR (seconds): " << flush; cin >> inTR;
	if (mdl == Model::Finite) {
		if (prompt) cout << "Enter RF Pulse Length (seconds): " << flush; cin >> inTrf;
	}
	if (prompt) cout << "Enter " << inAngles.size() << " Flip-angles (degrees): " << flush;
	for (ArrayXd::Index i = 0; i < inAngles.size(); i++) cin >> inAngles(i);
	string temp; getline(cin, temp); // Just to eat the newline
	if (use_weights) {
		if (prompt) cout << "Enter weighting: " << flush; cin >> inWeight;
		string temp; getline(cin, temp); // Just to eat the newline
	}
	shared_ptr<Signal> f;
	switch (mdl) {
		case Model::Simple: f = make_shared<SSFP>(nC, inAngles * M_PI / 180., inTR, inPhases * M_PI / 180., inWeight); break;
		case Model::Echo:   f = make_shared<SSFPEcho>(nC, inAngles * M_PI / 180., inTR, inPhases * M_PI / 180., inWeight); break;
		case Model::Finite: f = make_shared<SSFPFinite>(nC, inAngles * M_PI / 180., inTR, inTrf, inPhases * M_PI / 180., inWeight); break;
	}
	return f;
}

#ifdef AGILENT
shared_ptr<Signal> procparseSPGR(const Agilent::ProcPar &pp, const Components nC, const Model mdl,
										const bool prompt, const bool use_weights) {
	double inTR = 0., inTrf = 0., inTE = 0., inWeight = 1.;
	ArrayXd inAngles = pp.realValues("flip1");
	inTR = pp.realValue("tr");
	inAngles = pp.realValues("flip1");
	inTE = pp.realValue("te");         // Just read these anyway to save a switch
	inTrf = pp.realValue("p1") / 1.e6; // p1 is in microseconds
	if (use_weights) {
		if (prompt) cout << "Enter weighting: " << flush; cin >> inWeight;
		string temp; getline(cin, temp); // Just to eat the newline
	}
	shared_ptr<Signal> f;
	switch (mdl) {
		case Model::Simple: f = make_shared<SPGR>(nC, inAngles * M_PI / 180., inTR, inWeight); break;
		case Model::Echo:   f = make_shared<SPGREcho>(nC, inAngles * M_PI / 180, inTR, inTE, inWeight); break;
		case Model::Finite: f = make_shared<SPGRFinite>(nC, inAngles * M_PI / 180, inTR, inTrf, inTE, inWeight); break;
	}
	return f;
}

shared_ptr<Signal> procparseSSFP(const Agilent::ProcPar &pp, const Components nC, const Model mdl,
										const bool prompt, const bool use_weights) {
	double inTR = 0., inTrf = 0., inWeight = 1.;
	ArrayXd inPhases, inAngles;
	inPhases = pp.realValues("rfphase");
	inTR = pp.realValue("tr");
	inAngles = pp.realValues("flip1");
	inTrf = pp.realValue("p1") / 1.e6; // p1 is in microseconds
	if (use_weights) {
		if (prompt) cout << "Enter weighting: " << flush; cin >> inWeight;
		string temp; getline(cin, temp); // Just to eat the newline
	}
	shared_ptr<Signal> f;
	switch (mdl) {
		case Model::Simple: f = make_shared<SSFP>(nC, inAngles * M_PI / 180., inTR, inPhases * M_PI / 180., inWeight); break;
		case Model::Echo:   f = make_shared<SSFPEcho>(nC, inAngles * M_PI / 180., inTR, inPhases * M_PI / 180., inWeight); break;
		case Model::Finite: f = make_shared<SSFPFinite>(nC, inAngles * M_PI / 180., inTR, inTrf, inPhases * M_PI / 180., inWeight); break;
	}
	return f;
}
#endif

//******************************************************************************
#pragma mark DESPOTFunctor
//******************************************************************************
DESPOTFunctor::DESPOTFunctor(vector<shared_ptr<Signal>> &signals_in,
			                 const FieldStrength tesla, const OffRes offRes,
			                 const Scaling s, const bool debug) :
	m_signals(signals_in),
	m_fieldStrength(tesla), m_offRes(offRes),
	m_scaling(s), m_debug(debug),
	m_f0(0), m_B1(1)
{
	m_actual.reserve(m_signals.size());
	m_theory.reserve(m_signals.size());
	m_nV = 0;
	for (size_t i = 0; i < m_signals.size(); i++) {
		m_actual.emplace_back(m_signals.at(i)->size());
		m_theory.emplace_back(m_signals.at(i)->size());
		m_nV += m_signals.at(i)->size();
	}
}

const size_t DESPOTFunctor::nOffRes() const {
	switch (m_offRes) {
		case OffRes::Map: return 1;
		case OffRes::MapLoose: return 1;
		case OffRes::Fit: return 1;
		case OffRes::FitSym: return 1;
	}
}

const ArrayXXd DESPOTFunctor::offResBounds() const {
	double minTR = numeric_limits<double>::max();
	for (auto &s : m_signals) {
		if (s->m_TR < minTR)
			minTR = s->m_TR;
	}
	ArrayXXd b(1, 2);
	switch (m_offRes) {
		case OffRes::Map:             b << m_f0, m_f0; break;
		case OffRes::MapLoose:        b << m_f0 * 0.9, m_f0 * 1.1; break;
		case OffRes::Fit:          b << -0.5 / minTR, 0.5 / minTR; break;
		case OffRes::FitSym: b << 0., 0.5 / minTR; break;
	}
	return b;
}

const size_t DESPOTFunctor::nPD() const {
	switch (m_scaling) {
		case Scaling::Global:     return 1;
		case Scaling::NormToMean: return 0;
	}
}

const ArrayXXd DESPOTFunctor::PDBounds() const {
	ArrayXXd b(nPD(), 2);
	for (size_t i = 0; i < nPD(); i++) {
		b(i, 0) = 1.e4;
		b(i, 1) = 5.e6;
	}
	return b;
}

const ArrayXd DESPOTFunctor::actual() const {
	ArrayXd v(values());
	int index = 0;
	if (m_debug) cout << __PRETTY_FUNCTION__ << endl;
	for (size_t i = 0; i < m_actual.size(); i++) {
		v.segment(index, m_actual.at(i).size()) = m_actual.at(i);
		index += m_actual.at(i).size();
		if (m_debug) cout << v.transpose() << endl;
	}
	return v;
}

int DESPOTFunctor::operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) {
	eigen_assert(diffs.size() == values());
	ArrayXd t = theory(params);
	ArrayXd s = actual();
	diffs = t - s;
	if (m_debug) {
		cout << endl << __PRETTY_FUNCTION__ << endl;
		cout << "p:      " << params.transpose() << endl;
		cout << "t:      " << t.transpose() << endl;
		cout << "s:      " << s.transpose() << endl;
		cout << "Diffs:  " << diffs.transpose() << endl;
		cout << "Sum:    " << diffs.square().sum() << endl;
	}
	return 0;
}

ArrayXd DESPOTFunctor::weights() const {
	ArrayXd w(values());
	ArrayXd::Index index = 0;
	for (auto &s : m_signals) {
		w.segment(index, s->size()).setConstant(s->m_weight);
		index += s->size();
	}
	return w;
}