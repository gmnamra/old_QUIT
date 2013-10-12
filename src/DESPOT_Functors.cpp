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
SignalFunctor::SignalFunctor(const ArrayXd &flip, const double TR, const double weight, const Components nC) :
	m_flip(flip), m_TR(TR), m_weight(weight), m_nC(nC) {}

SPGR_Functor::SPGR_Functor(const ArrayXd &flip, const double TR, const double weight, const Components nC) :
	SignalFunctor(flip, TR, weight, nC) {}
ArrayXd SPGR_Functor::signal(const VectorXd &p, const double B1, double f0) const {
	switch (m_nC) {
		case (Components::One) : return SigMag(One_SPGR(p, m_flip, m_TR, B1));
		case (Components::Two) : return SigMag(Two_SPGR(p, m_flip, m_TR, B1));
		case (Components::Three) : return SigMag(Three_SPGR(p, m_flip, m_TR, B1));
	}
}

SPGR_Finite_Functor::SPGR_Finite_Functor(const ArrayXd &flip, const double TR, double Trf, double TE, const double weight, const Components nC) :
	SignalFunctor(flip, TR, weight, nC), m_Trf(Trf), m_TE(TE) {}
ArrayXd SPGR_Finite_Functor::signal(const VectorXd &p, const double B1, double f0) const {
	switch (m_nC) {
		case (Components::One) : return SigMag(One_SSFP_Finite(p, m_flip, true, m_TR, m_Trf, m_TE, 0, B1, f0));
		case (Components::Two) : return SigMag(Two_SSFP_Finite(p, m_flip, true, m_TR, m_Trf, m_TE, 0, B1, f0));
		case (Components::Three) : return SigMag(Three_SSFP_Finite(p, m_flip, true, m_TR, m_Trf, m_TE, 0, B1, f0));
	}
}


SSFP_Functor::SSFP_Functor(const ArrayXd &flip, const double TR, const ArrayXd &phases, const double weight, const Components nC) :
	SignalFunctor(flip, TR, weight, nC), m_phases(phases) {}
size_t SSFP_Functor::size() const {
	return m_flip.rows() * m_phases.rows();
}
ArrayXd SSFP_Functor::signal(const VectorXd &p, const double B1, double f0) const {
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

SSFP_Finite_Functor::SSFP_Finite_Functor(const ArrayXd &flip, const double TR, const double Trf, const ArrayXd &phases, const double weight, const Components nC) :
	SSFP_Functor(flip, TR, phases, weight, nC), m_Trf(Trf) {}
ArrayXd SSFP_Finite_Functor::signal(const VectorXd &p, const double B1, double f0) const {
	ArrayXd s(size());
	ArrayXd::Index start = 0;
	for (ArrayXd::Index i = 0; i < m_phases.rows(); i++) {
		switch (m_nC) {
			case (Components::One) : s.segment(start, m_flip.rows()) = SigMag(One_SSFP_Finite(p, m_flip, false, m_TR, m_Trf, 0, m_phases(i), B1, f0)); break;
			case (Components::Two) : s.segment(start, m_flip.rows()) = SigMag(Two_SSFP_Finite(p, m_flip, false, m_TR, m_Trf, 0, m_phases(i), B1, f0)); break;
			case (Components::Three) : s.segment(start, m_flip.rows()) = SigMag(Three_SSFP_Finite(p, m_flip, false, m_TR, m_Trf, 0, m_phases(i), B1, f0)); break;
		}
		start += m_flip.rows();
	}
	return s;
}

shared_ptr<SignalFunctor> parseSPGR(const Nifti &img, const bool prompt, const Components nC,
                                    const bool use_finite, const bool use_weights) {
	double inTR = 0., inTrf = 0., inTE = 0., inWeight = 1.;
	ArrayXd inAngles(img.dim(4));
	#ifdef AGILENT
	Agilent::ProcPar pp;
	if (ReadPP(img, pp)) {
		inTR = pp.realValue("tr");
		inAngles = pp.realValues("flip1");
		if (use_finite) {
			inTE = pp.realValue("te");
			inTrf = pp.realValue("p1") / 1.e6; // p1 is in microseconds
		}
	} else
	#endif
	{
		if (prompt) cout << "Enter TR (seconds): " << flush; cin >> inTR;
		if (use_finite) {
			if (prompt) cout << "Enter RF Pulse Length (seconds): " << flush; cin >> inTrf;
			if (prompt) cout << "Enter TE (seconds): " << flush; cin >> inTE;
		}
		if (prompt) cout << "Enter " << inAngles.size() << " Flip-angles (degrees): " << flush;
		for (int i = 0; i < inAngles.size(); i++) cin >> inAngles[i];
		string temp; getline(cin, temp); // Just to eat the newline
	}
	if (use_weights) {
		if (prompt) cout << "Enter weighting: " << flush; cin >> inWeight;
		string temp; getline(cin, temp); // Just to eat the newline
	}
	shared_ptr<SignalFunctor> f;
	if (use_finite)
		f = make_shared<SPGR_Finite_Functor>(inAngles * M_PI / 180, inTR, inTrf, inTE, inWeight, nC);
	else
		f = make_shared<SPGR_Functor>(inAngles * M_PI / 180., inTR, inWeight, nC);
	return f;
}

shared_ptr<SignalFunctor> parseSSFP(const Nifti &img, const bool prompt, const Components nC,
                                    const bool use_finite, const bool use_weights) {
	double inTR = 0., inTrf = 0., inWeight = 1.;
	ArrayXd inPhases, inAngles;
	#ifdef AGILENT
	Agilent::ProcPar pp;
	if (ReadPP(img, pp)) {
		inPhases = pp.realValues("rfphase");
		inTR = pp.realValue("tr");
		inAngles = pp.realValues("flip1");
		if (use_finite) {
			inTrf = pp.realValue("p1") / 1.e6; // p1 is in microseconds
		}
	} else
	#endif
	{
		size_t nPhases = 0;
		if (prompt) cout << "Enter number of phase-cycling patterns: " << flush; cin >> nPhases;
		inPhases.resize(nPhases);
		inAngles.resize(img.dim(4) / nPhases);
		if (prompt) cout << "Enter " << nPhases << " phase-cycles (degrees): " << flush;
		for (size_t i = 0; i < nPhases; i++) cin >> inPhases(i);
		if (prompt) cout << "Enter TR (seconds): " << flush; cin >> inTR;
		if (use_finite) {
			if (prompt) cout << "Enter RF Pulse Length (seconds): " << flush; cin >> inTrf;
		}
		if (prompt) cout << "Enter " << inAngles.size() << " Flip-angles (degrees): " << flush;
		for (ArrayXd::Index i = 0; i < inAngles.size(); i++) cin >> inAngles(i);
		string temp; getline(cin, temp); // Just to eat the newline
	}
	if (use_weights) {
		if (prompt) cout << "Enter weighting: " << flush; cin >> inWeight;
		string temp; getline(cin, temp); // Just to eat the newline
	}
	shared_ptr<SignalFunctor> f;
	if (use_finite)
		f = make_shared<SSFP_Finite_Functor>(inAngles * M_PI / 180., inTR, inTrf, inPhases * M_PI / 180., inWeight, nC);
	else
		f = make_shared<SSFP_Functor>(inAngles * M_PI / 180., inTR, inPhases * M_PI / 180., inWeight, nC);
	return f;
}