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
#pragma mark Signal Functors
//******************************************************************************
class Signal {
	public:
		enum class Components { One, Two, Three };
		static const string to_string(const Components& c);
		
		double m_TR;
		ArrayXd m_flip;
		Signal();
		Signal(const ArrayXd &flip, const double TR);
		ArrayXd B1flip(const double B1) const;
		virtual ArrayXcd signal(const Components nC, const VectorXd &p, const double B1 = 1.) const = 0;
		virtual size_t size() const { return m_flip.rows(); }
		virtual void write(ostream &os) const = 0;
		virtual string name() const = 0;
};
ostream& operator<<(ostream& os, const Signal& s);

class SPGRSimple : public Signal {
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
class SSFPSimple : public Signal {
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
class SSFPEllipse : public Signal {
	public:
		SSFPEllipse(const bool prompt = false, const Agilent::ProcPar &pp = Agilent::ProcPar());
		ArrayXcd signal(const Components nC, const VectorXd &p, const double B1 = 1.) const override;
		void write(ostream& os) const override;
		string name() const override { return "SSFP_Ellipse"; };
};

enum class SignalType { SPGR, SPGR_Finite, SSFP, SSFP_Finite, SSFP_Ellipse };
enum class OffRes { Fit, FitSym, Map }; // Put this here so mcdespot and despot2fm can access it

class Model {
public:
	enum class Scaling { None, NormToMean };
	enum class FieldStrength { Three, Seven, User };
	
	static const string to_string(const FieldStrength& f);
	static const string to_string(const Scaling &p);

	Signal::Components m_nC;
	Scaling m_scaling;
	vector<shared_ptr<Signal>> m_signals;
	
	Model(const Signal::Components c, const Scaling s);
	friend ostream& operator<<(ostream& os, const Model& m);

	const ArrayXcd signal(const VectorXd &p, const double B1) const;
	size_t nSignals() const;
	size_t size() const;
	
	const size_t nParameters() const;
	const bool validParameters(const VectorXd &params) const;
	const vector<string> &names() const;
	const ArrayXXd bounds(const FieldStrength f) const;
	
	ArrayXcd loadSignals(vector<MultiArray<complex<float>, 4>> &sigs, const size_t i, const size_t j, const size_t k) const;
	
	void addSignal(const SignalType &st, const bool prompt = false, const Agilent::ProcPar &pp = Agilent::ProcPar());
};

#endif
