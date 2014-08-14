/*
 *  Model.h
 *
 *  Created by Tobias Wood on 14/11/2012.
 *  Copyright (c) 2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef DESPOT_Model
#define DESPOT_Model

#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <Eigen/Dense>

#include "QUIT/QUIT.h"
#include "DESPOT.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
// Some convenience Enums
//******************************************************************************
enum class Components { One, Two, Three };
const string to_string(const Components& c);

enum class Scale { None, NormToMean };
const string to_string(const Scale &p);

enum class FieldStrength { Three, Seven, User };
const string to_string(const FieldStrength& f);
//******************************************************************************
// Sequence Functors
//******************************************************************************
class Sequence {
	public:
		double m_TR;
		ArrayXd m_flip;
		Sequence();
		Sequence(const ArrayXd &flip, const double TR);
		ArrayXd B1flip(const double B1) const;
		virtual ArrayXcd signal(const Components nC, const VectorXd &p, const double B1 = 1.) const = 0;
		virtual size_t size() const { return m_flip.rows(); }
		virtual void write(ostream &os) const = 0;
		virtual string name() const = 0;
};
ostream& operator<<(ostream& os, const Sequence& s);

class SPGRSimple : public Sequence {
	public:
		SPGRSimple(const ArrayXd &flip, const double TR);
		SPGRSimple(const bool prompt = false, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(const Components nC, const VectorXd &p, const double B1 = 1.) const override;
		void write(ostream& os) const override;
		string name() const override { return "SPGR"; } ;
};
class SPGRFinite : public SPGRSimple {
	public:
		double m_Trf, m_TE;
		SPGRFinite(const ArrayXd &flip, const double TR, const double Trf, const double TE);
		SPGRFinite(const bool prompt = false, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(const Components nC, const VectorXd &p, const double B1 = 1.) const override;
		void write(ostream& os) const override;
		string name() const override { return "SPGR_Finite"; } ;
};
class SSFPSimple : public Sequence {
	public:
		ArrayXd m_phases;
		SSFPSimple(const ArrayXd &flip, const double TR, const ArrayXd &phases);
		SSFPSimple(const bool prompt = false, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(const Components nC, const VectorXd &p, const double B1 = 1.) const override;
		size_t size() const override;
		void write(ostream& os) const override;
		string name() const override { return "SSFP"; } ;
};
class SSFPFinite : public SSFPSimple {
	public:
		double m_Trf;
		SSFPFinite(const ArrayXd &flip, const double TR, const double Trf, const ArrayXd &phases);
		SSFPFinite(const bool prompt = false, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(const Components nC, const VectorXd &p, const double B1 = 1.) const override;
		void write(ostream& os) const override;
		string name() const override { return "SSFP_Finite"; } ;
};
class SSFPEllipse : public Sequence {
	public:
		SSFPEllipse(const bool prompt = false, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(const Components nC, const VectorXd &p, const double B1 = 1.) const override;
		void write(ostream& os) const override;
		string name() const override { return "SSFP_Ellipse"; };
};

enum class SequenceType { SPGR, SPGR_Finite, SSFP, SSFP_Finite, SSFP_Ellipse };
enum class OffRes { Fit, FitSym, Map }; // Put this here so mcdespot and despot2fm can access it

class Sequences {
private:
	Components m_nC;
	Scale m_scaling;
	vector<shared_ptr<Sequence>> m_sequences;

public:
	Sequences(const Components c, const Scale s);
	friend ostream& operator<<(ostream& os, const Sequences& m);

	size_t count() const;
	shared_ptr<Sequence> sequence(const size_t i) const;
	const ArrayXcd signal(const size_t i, const VectorXd &p, const double B1) const;

	size_t combinedSize() const;
	const ArrayXcd combinedSignal(const VectorXd &p, const double B1) const;

	const size_t nParameters() const;
	const bool validParameters(const VectorXd &params) const;
	const vector<string> &names() const;
	const ArrayXXd bounds(const FieldStrength f) const;
	
	ArrayXcd loadSignals(vector<QUIT::MultiArray<complex<float>, 4>> &sigs, const size_t i, const size_t j, const size_t k) const;
	
	void addSequence(const SequenceType &st, const bool prompt = false, const Agilent::ProcPar &pp = Agilent::ProcPar());
};

#endif
