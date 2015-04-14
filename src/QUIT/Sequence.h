/*
 *  Sequence.h
 *
 *  Created by Tobias Wood on 14/11/2012.
 *  Copyright (c) 2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef DESPOT_SEQUENCE
#define DESPOT_SEQUENCE

#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <Eigen/Dense>

#include "QUIT/QUIT.h"
#include "SignalEquations.h"
#include "Models.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
// Some convenience Enums
//******************************************************************************
enum class Scale { None, NormToMean };
const string to_string(const Scale &p);

//******************************************************************************
// SteadyState Functors
//******************************************************************************
class SequenceBase {
	public:
		virtual ArrayXcd signal(const shared_ptr<Model> m, const VectorXd &p) const = 0;
		virtual size_t size() const = 0;
		virtual void write(ostream &os) const = 0;
		virtual string name() const = 0;
};
ostream& operator<<(ostream& os, const SequenceBase& s);

class MultiEcho : public SequenceBase {
	public:
		double m_ESP;
		ArrayXd m_TE;

		MultiEcho();
		MultiEcho(const ArrayXd &te);
		MultiEcho(const bool prompt, const Agilent::ProcPar &pp = Agilent::ProcPar());

		size_t size() const override { return m_TE.rows(); };
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream &os) const override;
		string name() const override { return "MultiEcho"; };
};

class SteadyState : public SequenceBase {
	public:
		double m_TR;
		ArrayXd m_flip;

		SteadyState();
		SteadyState(const ArrayXd &flip, const double TR);

		virtual size_t size() const override { return angles() * phases(); };
		virtual size_t angles() const { return m_flip.rows(); }
		virtual size_t phases() const { return 1; };
};

class SPGRSimple : public SteadyState {
	public:
		SPGRSimple(const ArrayXd &flip, const double TR);
		SPGRSimple(const bool prompt, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream &os) const override;
		string name() const override { return "SPGR"; };
};
class SPGRFinite : public SPGRSimple {
	public:
		double m_Trf, m_TE;
		SPGRFinite(const ArrayXd &flip, const double TR, const double Trf, const double TE);
		SPGRFinite(const bool prompt, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream &os) const override;
		string name() const override { return "SPGR_Finite"; };
};
class MPRAGE : public SteadyState {
	public:
		ArrayXd m_TI;
		double m_TD;
		int m_N;
		MPRAGE() : SteadyState() {};
		MPRAGE(const ArrayXd &TI, const double TD, const double TR, const int N, const double flip);
		MPRAGE(const bool prompt, const Agilent::ProcPar &pp = Agilent::ProcPar());
		size_t size() const override { return m_TI.size(); };
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream &os) const override;
		string name() const override { return "MPRAGE"; };
};

// Special class for GE IRSPGR, for backwards compatibility
class IRSPGR : public MPRAGE {
	public:
		IRSPGR(const bool prompt, const Agilent::ProcPar &pp = Agilent::ProcPar());
		string name() const override { return "IRSPGR"; };
};

class SSFPSimple : public SteadyState {
	public:
		ArrayXd m_phases;
		SSFPSimple(const ArrayXd &flip, const double TR, const ArrayXd &phases);
		SSFPSimple(const bool prompt, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		size_t phases() const override;
		void write(ostream& os) const override;
		string name() const override { return "SSFP"; } ;
};
class SSFPFinite : public SSFPSimple {
	public:
		double m_Trf;
		SSFPFinite(const ArrayXd &flip, const double TR, const double Trf, const ArrayXd &phases);
		SSFPFinite(const bool prompt, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream& os) const override;
		string name() const override { return "SSFP_Finite"; } ;
};
class SSFPEllipse : public SteadyState {
	public:
		SSFPEllipse(const bool prompt, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
		void write(ostream& os) const override;
		string name() const override { return "SSFP_Ellipse"; };
};

enum class OffRes { Fit, FitSym, Map }; // Put this here so mcdespot and despot2fm can access it

class SequenceGroup : public SequenceBase {
private:
	Scale m_scaling;
	vector<shared_ptr<SteadyState>> m_sequences;

public:
	SequenceGroup(const Scale s);
	void write(ostream &os) const override;
	string name() const override { return "Sequences"; } ;

	size_t count() const;
	shared_ptr<SteadyState> sequence(const size_t i) const;
	vector<shared_ptr<SteadyState>> &sequences();

	size_t size() const override;
	ArrayXcd signal(shared_ptr<Model> m, const VectorXd &par) const override;
	
	double minTR() const;
	ArrayXcd loadSignals(vector<QUIT::MultiArray<complex<float>, 4>> &sigs, const size_t i, const size_t j, const size_t k, bool needsFlip = false) const;
	
	void addSequence(const shared_ptr<SteadyState> &seq);
};

#endif
