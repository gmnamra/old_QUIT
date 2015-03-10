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

using namespace std;
using namespace Eigen;

//******************************************************************************
// Some convenience Enums
//******************************************************************************
enum class Pools { One, Two, Three };
const string to_string(const Pools& c);

enum class Scale { None, NormToMean };
const string to_string(const Scale &p);

enum class FieldStrength { Three, Seven, User };
const string to_string(const FieldStrength& f);

class PoolInfo {
public:
	static const size_t nParameters(const Pools p);
	static const bool ValidParameters(const Pools p, const VectorXd &params);
	static const vector<string> &Names(const Pools p);
	static const ArrayXXd Bounds(const Pools p, const FieldStrength f, const double TR);
};

//******************************************************************************
// Sequence Functors
//******************************************************************************
class SequenceBase {
	public:
		virtual ArrayXcd signal(const Pools c, const VectorXd &p, const double B1 = 1.) const = 0;
		virtual size_t size() const = 0;
		virtual void write(ostream &os) const = 0;
		virtual string name() const = 0;
};
ostream& operator<<(ostream& os, const SequenceBase& s);

class Sequence : public SequenceBase {
	public:
		double m_TR;
		ArrayXd m_flip;

		Sequence();
		Sequence(const ArrayXd &flip, const double TR);
		ArrayXd B1flip(const double B1) const;

		virtual size_t size() const override { return angles() * phases(); };
		virtual size_t angles() const { return m_flip.rows(); }
		virtual size_t phases() const { return 1; };
};


class SPGRSimple : public Sequence {
	public:
		SPGRSimple(const ArrayXd &flip, const double TR);
		SPGRSimple(const bool prompt = false, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(const Pools p, const VectorXd &par, const double B1 = 1.) const override;
		void write(ostream &os) const override;
		string name() const override { return "SPGR"; };
};
class SPGRFinite : public SPGRSimple {
	public:
		double m_Trf, m_TE;
		SPGRFinite(const ArrayXd &flip, const double TR, const double Trf, const double TE);
		SPGRFinite(const bool prompt = false, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(const Pools p, const VectorXd &par, const double B1 = 1.) const override;
		void write(ostream &os) const override;
		string name() const override { return "SPGR_Finite"; };
};
class MPRAGE : public Sequence {
	public:
		ArrayXd m_TI;
		double m_TD;
		int m_N;
		MPRAGE(const ArrayXd &TI, const double TD, const double TR, const int N, const double flip);
		MPRAGE(const bool prompt = false, const Agilent::ProcPar &pp = Agilent::ProcPar());
		size_t size() const override { return m_TI.size(); };
		ArrayXcd signal(const Pools p, const VectorXd &par, const double B1 = 1.) const override;
		void write(ostream &os) const override;
		string name() const override { return "MPRAGE"; };
};
class SSFPSimple : public Sequence {
	public:
		ArrayXd m_phases;
		SSFPSimple(const ArrayXd &flip, const double TR, const ArrayXd &phases);
		SSFPSimple(const bool prompt = false, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(const Pools p, const VectorXd &par, const double B1 = 1.) const override;
		size_t phases() const override;
		void write(ostream& os) const override;
		string name() const override { return "SSFP"; } ;
};
class SSFPFinite : public SSFPSimple {
	public:
		double m_Trf;
		SSFPFinite(const ArrayXd &flip, const double TR, const double Trf, const ArrayXd &phases);
		SSFPFinite(const bool prompt = false, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(const Pools p, const VectorXd &par, const double B1 = 1.) const override;
		void write(ostream& os) const override;
		string name() const override { return "SSFP_Finite"; } ;
};
class SSFPEllipse : public Sequence {
	public:
		SSFPEllipse(const bool prompt = false, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(const Pools p, const VectorXd &par, const double B1 = 1.) const override;
		void write(ostream& os) const override;
		string name() const override { return "SSFP_Ellipse"; };
};

enum class OffRes { Fit, FitSym, Map }; // Put this here so mcdespot and despot2fm can access it

class SequenceGroup : public SequenceBase {
private:
	Pools m_nC;
	Scale m_scaling;
	vector<shared_ptr<Sequence>> m_sequences;

public:
	SequenceGroup(const Scale s);
	void write(ostream &os) const override;
	string name() const override { return "Sequences"; } ;

	size_t count() const;
	shared_ptr<Sequence> sequence(const size_t i) const;
	vector<shared_ptr<Sequence>> &sequences();

	size_t size() const override;
	ArrayXcd signal(const Pools p, const VectorXd &par, const double B1) const override;
	
	double minTR() const;
	ArrayXcd loadSignals(vector<QUIT::MultiArray<complex<float>, 4>> &sigs, const size_t i, const size_t j, const size_t k, bool needsFlip = false) const;
	
	//void addSequence(const SequenceType &st, const bool prompt = false, const Agilent::ProcPar &pp = Agilent::ProcPar());
	void addSequence(const shared_ptr<Sequence> &seq);
};

#endif
