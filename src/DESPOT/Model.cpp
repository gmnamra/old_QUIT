/*
 *  Sequences.cpp
 *
 *  Created by Tobias Wood on 14/11/2012.
 *  Copyright (c) 2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "Model.h"

const string to_string(const Components& c) {
	switch (c) {
		case Components::One: return "1";
		case Components::Two: return "2";
		case Components::Three: return "3";
	}
};

const string to_string(const FieldStrength& f) {
	static const string f3{"3T"}, f7{"7T"}, fu{"User"};
	switch (f) {
		case FieldStrength::Three: return f3;
		case FieldStrength::Seven: return f7;
		case FieldStrength::User: return fu;
	}
}

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
Sequence::Sequence() {}
Sequence::Sequence(const ArrayXd &flip, const double TR) : m_flip(flip), m_TR(TR) {}

ostream& operator<<(ostream& os, const Sequence& s) {
	s.write(os);
	return os;
}

ArrayXd Sequence::B1flip(const double B1) const {
	return B1 * m_flip;
}

// Parameters are PD, T1, T2, f0
//                PD, T1_a, T2_a, T1_b, T2_b, tau_a, f_a, f0
//				  PD, T1_a, T2_a, T1_b, T2_b, T1_c, T2_c, tau_a, f_a, f_c, f0
SPGRSimple::SPGRSimple(const ArrayXd &flip, const double TR) : Sequence(flip, TR) {}
SPGRSimple::SPGRSimple(const bool prompt, const Agilent::ProcPar &pp) {
	if (pp) {
		m_flip = pp.realValues("flip1") * M_PI / 180.;
		m_TR = pp.realValue("tr");
	} else {
		size_t nFlip;
		if (prompt) cout << "Enter number of flip-angles: " << flush;
		QUIT::Read<size_t>::FromLine(cin, nFlip);
		ArrayXd inAngles(nFlip);
		if (prompt) cout << "Enter " << inAngles.size() << " flip-angles (degrees): " << flush;
		QUIT::Read<ArrayXd>::FromLine(cin, inAngles);
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

ArrayXcd SPGRSimple::signal(const Components nC, const VectorXd &p, const double B1) const {
	switch (nC) {
		case (Components::One) : return SigComplex(One_SPGR(B1flip(B1), m_TR, p[0], p[1]));
		case (Components::Two) : return SigComplex(Two_SPGR(B1flip(B1), m_TR, p[0], p[1], p[3], p[5], p[6]));
		case (Components::Three) : return SigComplex(Three_SPGR(B1flip(B1), m_TR, p[0], p[1], p[3], p[5], p[7], p[8], p[9]));
	}
}

SPGRFinite::SPGRFinite(const ArrayXd &flip, const double TR, const double Trf, const double TE) : SPGRSimple(flip, TR), m_Trf(Trf), m_TE(TE) {}
SPGRFinite::SPGRFinite(const bool prompt, const Agilent::ProcPar &pp) : SPGRSimple(prompt, pp) {
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

ArrayXcd SPGRFinite::signal(const Components nC, const VectorXd &p, const double B1) const {
	switch (nC) {
		case (Components::One) : return SigComplex(One_SSFP_Finite(B1flip(B1), true, m_TR, m_Trf, m_TE, 0, p[0], p[1], p[2], p[3]));
		case (Components::Two) : return SigComplex(Two_SSFP_Finite(B1flip(B1), true, m_TR, m_Trf, m_TE, 0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]));
		case (Components::Three) : return SigComplex(Three_SSFP_Finite(B1flip(B1), true, m_TR, m_Trf, m_TE, 0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10]));
	}
}


SSFPSimple::SSFPSimple(const ArrayXd &flip, const double TR, const ArrayXd &phases) : Sequence(flip, TR), m_phases(phases) {}
SSFPSimple::SSFPSimple(const bool prompt, const Agilent::ProcPar &pp) {
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
		QUIT::Read<ArrayXd>::FromLine(cin, inAngles);
		m_flip = inAngles * M_PI / 180.;
		size_t nPhases;
		if (prompt) cout << "Enter number of phase-cycles: " << flush;
		QUIT::Read<size_t>::FromLine(cin, nPhases);
		ArrayXd inPhases(nPhases);
		if (prompt) cout << "Enter " << inPhases.size() << " phase-cycles (degrees): " << flush;
		QUIT::Read<ArrayXd>::FromLine(cin, inPhases);
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

size_t SSFPSimple::size() const { return m_flip.rows() * m_phases.rows(); }
ArrayXcd SSFPSimple::signal(const Components nC, const VectorXd &p, const double B1) const {
	ArrayXcd s(size());
	ArrayXcd::Index start = 0;
	for (ArrayXcd::Index i = 0; i < m_phases.rows(); i++) {
		switch (nC) {
			case (Components::One) : s.segment(start, m_flip.rows()) = SigComplex(One_SSFP(B1flip(B1), m_TR, m_phases(i), p[0], p[1], p[2], p[3])); break;
			case (Components::Two) : s.segment(start, m_flip.rows()) = SigComplex(Two_SSFP(B1flip(B1), m_TR, m_phases(i), p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7])); break;
			case (Components::Three) : s.segment(start, m_flip.rows()) = SigComplex(Three_SSFP(B1flip(B1), m_TR, m_phases(i), p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10])); break;
		}
		start += m_flip.rows();
	}
	return s;
}

SSFPFinite::SSFPFinite(const ArrayXd &flip, const double TR, const double Trf, const ArrayXd &phases) : SSFPSimple(flip, TR, phases), m_Trf(Trf) {}
SSFPFinite::SSFPFinite(const bool prompt, const Agilent::ProcPar &pp) : SSFPSimple(prompt, pp) {
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

ArrayXcd SSFPFinite::signal(const Components nC, const VectorXd &p, const double B1) const {
	ArrayXcd s(size());
	ArrayXcd::Index start = 0;
	for (ArrayXcd::Index i = 0; i < m_phases.rows(); i++) {
		switch (nC) {
			case (Components::One) : s.segment(start, m_flip.rows()) = SigComplex(One_SSFP_Finite(B1flip(B1), false, m_TR, m_Trf, 0., m_phases(i), p[0], p[1], p[2], p[3])); break;
			case (Components::Two) : s.segment(start, m_flip.rows()) = SigComplex(Two_SSFP_Finite(B1flip(B1), false, m_TR, m_Trf, 0., m_phases(i), p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7])); break;
			case (Components::Three) : s.segment(start, m_flip.rows()) = SigComplex(Three_SSFP_Finite(B1flip(B1), false, m_TR, m_Trf, 0., m_phases(i), p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10])); break;
		}
		start += m_flip.rows();
	}
	return s;
}

SSFPEllipse::SSFPEllipse(const bool prompt, const Agilent::ProcPar &pp) {
	if (pp) {
		m_flip = pp.realValues("flip1") * M_PI / 180.;
		m_TR = pp.realValue("tr");
	} else {
		size_t nFlip;
		if (prompt) cout << "Enter number of flip-angles: " << flush;
		QUIT::Read<size_t>::FromLine(cin, nFlip);
		ArrayXd inAngles(nFlip);
		if (prompt) cout << "Enter " << inAngles.size() << " flip-angles (degrees): " << flush;
		QUIT::Read<ArrayXd>::FromLine(cin, inAngles);
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

ArrayXcd SSFPEllipse::signal(const Components nC, const VectorXd &p, const double B1) const {
	switch (nC) {
		case (Components::One) : return SigComplex(One_SSFP_Ellipse(B1*m_flip, m_TR, p[0], p[1], p[2], p[3])); break;
		case (Components::Two) : throw(logic_error("Two component SSFP Ellipse not implemented.")); break;
		case (Components::Three) : throw(logic_error("Three component SSFP Ellipse not implemented.")); break;
	}
}

//******************************************************************************
#pragma mark Sequences Class
//******************************************************************************
Sequences::Sequences(const Components c, const Scale s) : m_nC(c), m_scaling(s) {}

ostream& operator<<(ostream &os, const Sequences& m) {
	os << "Sequences Parameters: " << m.nParameters() << endl;
	os << "Names:\t"; for (auto& n : m.names()) os << n << "\t"; os << endl;
	os << "Signals: " << m.m_sequences.size() << "\tCombined size: " << m.combinedSize() << endl;
	for (auto& sig : m.m_sequences)
		os << *sig;
	return os;
}

size_t Sequences::count() const {
	return m_sequences.size();
}

shared_ptr<Sequence> Sequences::sequence(const size_t i) const {
	return m_sequences.at(i);
}

const ArrayXcd Sequences::signal(const size_t i, const VectorXd &p, const double B1) const {
	return m_sequences.at(i)->signal(m_nC, p, B1);
}

size_t Sequences::combinedSize() const {
	size_t sz = 0;
	for (auto& sig : m_sequences)
		sz += sig->size();
	return sz;
}

const ArrayXcd Sequences::combinedSignal(const VectorXd &p, const double B1) const {
	ArrayXcd result(combinedSize());
	size_t start = 0;
	for (auto &sig : m_sequences) {
		ArrayXcd thisResult = sig->signal(m_nC, p, B1);
		switch (m_scaling) {
			case Scale::None :       break;
			case Scale::NormToMean : thisResult /= thisResult.abs().mean();
		}
		result.segment(start, sig->size()) = thisResult;
		start += sig->size();
	}
	return result;
}

const size_t Sequences::nParameters() const {
	switch (m_nC) {
		case Components::One: return 4;
		case Components::Two: return 8;
		case Components::Three: return 11;
	}
}

const vector<string> &Sequences::names() const {
	static vector<string> n1 {"PD", "T1", "T2", "f0"},
	                      n2 {"PD", "T1_a", "T2_a", "T1_b", "T2_b", "tau_a", "f_a", "f0"},
				          n3 {"PD", "T1_a", "T2_a", "T1_b", "T2_b", "T1_c", "T2_c", "tau_a", "f_a", "f_c", "f0"};
	
	switch (m_nC) {
		case Components::One: return n1;
		case Components::Two: return n2;
		case Components::Three: return n3;
	}
}

const ArrayXXd Sequences::bounds(const FieldStrength f) const {
	size_t nP = nParameters();
	ArrayXXd b(nP, 2);
	switch (f) {
		case FieldStrength::Three:
			switch (m_nC) {
				case Components::One:   b.block(0, 0, nP - 1, 2) << 1.0, 1.0, 0.1, 4.5, 0.010, 2.500; break;
				case Components::Two:   b.block(0, 0, nP - 1, 2) << 1.0, 1.0, 0.1, 0.5, 0.001, 0.030, 0.700, 4.500, 0.050, 0.200, 0.025, 0.600, 0.00, 1.0; break;
				case Components::Three: b.block(0, 0, nP - 1, 2) << 1.0, 1.0, 0.1, 0.5, 0.001, 0.030, 0.700, 2.000, 0.050, 0.200, 3.000, 4.500, 1.50, 2.50, 0.025, 0.600, 0.0, 1.0, 0, 1.0; break;
			} break;
		case FieldStrength::Seven:
			switch (m_nC) {
				case Components::One:   b.block(0, 0, nP - 1, 2) << 1.0, 1.0, 0.1, 4.5, 0.010, 2.500; break;
				case Components::Two:   b.block(0, 0, nP - 1, 2) << 1.0, 1.0, 0.1, 0.5, 0.001, 0.025, 1.5, 4.5, 0.04, 0.20, 0.025, 0.600, 0.0, 1.0; break;
				case Components::Three: b.block(0, 0, nP - 1, 2) << 1.0, 1.0, 0.1, 0.5, 0.001, 0.025, 1.5, 2.5, 0.04, 0.20, 3.000, 4.500, 1.5, 2.5, 0.025, 0.600, 0.0, 1.0, 0.0, 1.0; break;
			} break;
		case FieldStrength::User:
			switch (m_nC) {
				case Components::One:   b.block(0, 0, nP - 1, 2).setZero(); break;
				case Components::Two:   b.block(0, 0, nP - 1, 2).setZero(); break;
				case Components::Three: b.block(0, 0, nP - 1, 2).setZero(); break;
			} break;
	}
	
	double minTR = numeric_limits<double>::max();
	for (auto &s : m_sequences) {
		if (s->m_TR < minTR)
			minTR = s->m_TR;
	}
	b.block(nP - 1, 0, 1, 2) << -0.5 / minTR, 0.5 / minTR;
	
	return b;
}

const bool Sequences::validParameters(const VectorXd &params) const {
	// Negative T1/T2 makes no sense
	if ((params[1] <= 0.) || (params[2] <= 0.))
		return false;
	
	switch (m_nC) {
		case Components::One : return true;
		case Components::Two :
			// Check that T1_a, T2_a < T1_b, T2_b and that f_a makes sense
			if ((params[1] < params[3]) &&
				(params[2] < params[4]) &&
				(params[6] <= 1.0))
				return true;
			else
				return false;
		case Components::Three :
			// Check that T1/2_a < T1/2_b < T1/2_c and that f_a + f_c makes sense
			if ((params[1] < params[3]) &&
				(params[2] < params[4]) &&
				(params[3] < params[5]) &&
				(params[4] < params[6]) &&
				((params[8] + params[9]) <= 1.0))
				return true;
			else
				return false;
	}
}

ArrayXcd Sequences::loadSignals(vector<QUIT::MultiArray<complex<float>, 4>> &sigs, const size_t i, const size_t j, const size_t k) const {
	ArrayXcd signal(combinedSize());
	size_t start = 0;
	for (size_t s = 0; s < m_sequences.size(); s++) {
		ArrayXcd thisSig = sigs.at(s).slice<1>({i,j,k,0},{0,0,0,-1}).asArray().cast<complex<double>>();
		if (m_scaling == Scale::NormToMean)
			thisSig /= thisSig.abs().mean();
		signal.segment(start, thisSig.rows()) = thisSig;
		start += thisSig.rows();
	}
	return signal;
}


void Sequences::addSequence(const SequenceType &st, const bool prompt, const Agilent::ProcPar &pp) {
	switch (st) {
		case SequenceType::SPGR:         m_sequences.push_back(make_shared<SPGRSimple>(prompt, pp)); break;
		case SequenceType::SPGR_Finite:  m_sequences.push_back(make_shared<SPGRFinite>(prompt, pp)); break;
		case SequenceType::SSFP:         m_sequences.push_back(make_shared<SSFPSimple>(prompt, pp)); break;
		case SequenceType::SSFP_Finite:  m_sequences.push_back(make_shared<SSFPFinite>(prompt, pp)); break;
		case SequenceType::SSFP_Ellipse: m_sequences.push_back(make_shared<SSFPEllipse>(prompt, pp)); break;
	}
}
