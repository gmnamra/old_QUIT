/*
 *  DESPOT_Functors.h
 *
 *  Created by Tobias Wood on 16/08/2012.
 *  Copyright (c) 2012-2013 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef DESPOT_Functors_h
#define DESPOT_Functors_h

#include <vector>
#include <array>
#include <map>
#include <iostream>
#include <exception>
#include <string>
#include <memory>
#include <Eigen/Dense>

#include "DESPOT.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
#pragma mark Signal Functors
//******************************************************************************
enum class Components {
	One, Two, Three
};

enum class Model {
	Simple, Echo, Finite
};

class SignalFunctor {
	public:
		Components m_nC;
		double m_TR, m_weight;
		ArrayXd m_flip;
		SignalFunctor(const Components nC, const ArrayXd &flip, const double TR, const double weight = 1.);
		virtual ArrayXd signal(const VectorXd &p, const double B1 = 1., const double f0 = 0.) const = 0;
		virtual size_t size() const { return m_flip.rows(); }
};

class SPGR_Functor : public SignalFunctor {
	public:
		SPGR_Functor(const Components nC, const ArrayXd &flip, const double TR, const double weight = 1.);
		ArrayXd signal(const VectorXd &p, const double B1 = 1., const double f0 = 0.) const override;
};
class SPGR_Echo_Functor : public SignalFunctor {
	public:
		double m_TE;
		SPGR_Echo_Functor(const Components nC, const ArrayXd &flip, const double TR, const double TE, const double weight = 1.);
		ArrayXd signal(const VectorXd &p, const double B1 = 1., const double f0 = 0.) const override;
};
class SPGR_Finite_Functor : public SignalFunctor {
	public:
		double m_Trf, m_TE;
		SPGR_Finite_Functor(const Components nC, const ArrayXd &flip, const double TR, const double Trf, const double TE, const double weight = 1.);
		ArrayXd signal(const VectorXd &p, const double B1 = 1., const double f0 = 0.) const override;
};
class SSFP_Functor : public SignalFunctor {
	public:
		ArrayXd m_phases;
		SSFP_Functor(const Components nC, const ArrayXd &flip, const double TR, const ArrayXd &phases, const double weight = 1.);
		ArrayXd signal(const VectorXd &p, const double B1 = 1., const double f0 = 0.) const override;
		size_t size() const override;
};
class SSFP_Echo_Functor : public SSFP_Functor {
	public:
		SSFP_Echo_Functor(const Components nC, const ArrayXd &flip, const double TR, const ArrayXd &phases, const double weight = 1.);
		ArrayXd signal(const VectorXd &p, const double B1 = 1., const double f0 = 0.) const override;
};
class SSFP_Finite_Functor : public SSFP_Functor {
	public:
		double m_Trf, m_TE;
		SSFP_Finite_Functor(const Components nC, const ArrayXd &flip, const double TR, const double Trf, const ArrayXd &phases, const double weight = 1.);
		ArrayXd signal(const VectorXd &p, const double B1 = 1., const double f0 = 0.) const override;
};

//******************************************************************************
#pragma mark Parsing Functions
//******************************************************************************
shared_ptr<SignalFunctor> parseSPGR(const Components nC, const Model mdl, const size_t nFlip,
									const bool prompt, const bool use_weights);
shared_ptr<SignalFunctor> parseSSFP(const Components nC, const Model mdl, const size_t nVols,
                                    const bool prompt, const bool use_weights);
#ifdef AGILENT
shared_ptr<SignalFunctor> procparseSPGR(const Agilent::ProcPar &pp, const Components nC, const Model mdl,
										const bool prompt, const bool use_weights);
shared_ptr<SignalFunctor> procparseSSFP(const Agilent::ProcPar &pp, const Components nC, const Model mdl,
										const bool prompt, const bool use_weights);
#endif
//******************************************************************************
#pragma mark Optimisation Functor Base Class
//******************************************************************************
// From Nonlinear Tests in Eigen 
template<typename _Scalar, int NX=Dynamic, int NY=Dynamic>
class Functor {
	public:
		typedef _Scalar Scalar;
		enum {
			InputsAtCompileTime = NX,
			ValuesAtCompileTime = NY
		};
		typedef Matrix<Scalar,InputsAtCompileTime,1> InputType;
		typedef Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
		typedef Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;
		
		const long m_inputs, m_values;
		
		Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
		Functor(long inputs, long values) : m_inputs(inputs), m_values(values) {}
		Functor(long values) : m_inputs(InputsAtCompileTime), m_values(values) {}
		
		virtual ~Functor() {};
		
		virtual const long inputs() const { return m_inputs; }
		virtual const long values() const { return m_values; }
		
		virtual int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) = 0;
		
		virtual const ArrayXd theory(const Ref<VectorXd> &params) = 0;
		virtual const ArrayXd actual() const = 0;
};

//******************************************************************************
#pragma mark Basic DESPOT Functor
//******************************************************************************
class DESPOTFunctor : public Functor<double> {
	public:
		enum class FieldStrength {
			Three, Seven, Unknown
		};
		static const string to_string(const FieldStrength& f) {
			switch (f) {
				case FieldStrength::Three: return "3";
				case FieldStrength::Seven: return "7";
				case FieldStrength::Unknown: return "User";
			}
		}
		enum class Scaling {
			Global, NormToMean
		};
		static const string to_string(const Scaling &p) {
			switch (p) {
				case Scaling::Global: return "global";
				case Scaling::NormToMean: return "normalised to per signal mean";
			}
		}
		enum class OffRes {
			Map = 0, MapLoose, Fit, FitSym
		};
	
	protected:
		const FieldStrength m_fieldStrength;
		const OffRes m_offRes;
		const Scaling m_scaling;
		size_t m_nV;
		vector<shared_ptr<SignalFunctor>> m_signals;
		vector<ArrayXd> m_actual, m_theory;
		vector<string> m_names; // Subclasses responsible for initialising this
	
		const ArrayXXd offResBounds() const;
		const ArrayXXd PDBounds() const;
		
		virtual const ArrayXXd defaultBounds() = 0;
		
	public:
		bool m_debug;
		double m_f0, m_B1;
		const size_t nOffRes() const;
		const size_t nPD() const;
		
		virtual const size_t nP() const = 0;
		virtual const bool constraint(const VectorXd &params) const = 0;
		const long inputs() const override { return this->nP() + nOffRes() + nPD(); }
		const long values() const override { return m_nV; }
		
		DESPOTFunctor(vector<shared_ptr<SignalFunctor>> &signals_in,
					  const FieldStrength tesla, const OffRes offRes,
					  const Scaling s, const bool debug = false);
		
		const vector<string> &names() { return m_names; }
		shared_ptr<SignalFunctor> &signal(const size_t i) { return m_signals.at(i); }
		ArrayXd &actual(const size_t s) { return m_actual.at(s); }
		const ArrayXd actual() const;
		int operator()(const Ref<VectorXd> &params, Ref<ArrayXd> diffs) override;		
		ArrayXd weights() const;
};

//******************************************************************************
#pragma mark mcDESPOT Functor
//******************************************************************************
class mcDESPOT : public DESPOTFunctor {
	public:
		static const string to_string(const Components& c) {
			switch (c) {
				case Components::One: return "1";
				case Components::Two: return "2";
				case Components::Three: return "3";
			}
		};
	protected:
		const Components m_components;
		const bool m_finite;
	
	public:
		const size_t nP() const override {
			size_t n = 0;
			switch (m_components) {
				case Components::One: n = 2; break;
				case Components::Two: n = 6; break;
				case Components::Three: n = 9; break;
			}
			if (m_finite)
				n += 1;
			return n;
		}
		
		const ArrayXXd defaultBounds() override {
			ArrayXXd b(inputs(), 2);
			switch (m_fieldStrength) {
				case FieldStrength::Three:
					switch (m_components) {
						case Components::One:   b.block(0, 0, 2, 2) << 0.1, 4.0, 0.01, 1.5; break;
						case Components::Two:   b.block(0, 0, 6, 2) << 0.200, 0.350, 0.005, 0.015, 0.700, 2.000, 0.050, 0.120, 0.050, 0.200, 0.0, 0.5; break;
						case Components::Three: b.block(0, 0, 9, 2) << 0.200, 0.350, 0.005, 0.015, 0.700, 2.000, 0.050, 0.120, 3.500, 7.000, 3.000, 7.000, 0.050, 0.200, 0, 0.5, 0, 1.0; break;
					} break;
				case FieldStrength::Seven:
					switch (m_components) {
						case Components::One:   b.block(0, 0, 2, 2) << 0.1, 4.0, 0.01, 2.0; break;
						case Components::Two:   b.block(0, 0, 6, 2) << 0.1, 0.5, 0.001, 0.025, 1.0, 4.0, 0.04, 0.08, 0.01, 0.25, 0.001, 1.0; break;
						case Components::Three: b.block(0, 0, 9, 2) << 0.1, 0.5, 0.001, 0.025, 1.0, 2.5, 0.04, 0.08, 3., 4.5, 0.5, 2.0, 0.05, 0.200, 0.0, 0.5, 0.0, 1.0; break;
					} break;
				case FieldStrength::Unknown:
					switch (m_components) {
						case Components::One:   b.block(0, 0, 2, 2).setZero(); break;
						case Components::Two:   b.block(0, 0, 6, 2).setZero(); break;
						case Components::Three: b.block(0, 0, 9, 2).setZero(); break;
					} break;
			}
			if (m_finite) {
				b.block(nP() - 1, 0, 1, 2) << 0., 25.;
			}
			b.block(nP(), 0, nOffRes(), 2) = offResBounds();
			b.block(nP() + nOffRes(), 0, nPD(), 2) = PDBounds();
			return b;
		}
		
		const ArrayXd defaultThresholds() {
			ArrayXd m(inputs());
			switch (m_components) {
				case Components::One: m.head(nP()) << 0.05, 0.05; break;
				case Components::Two: m.head(nP()) << 0.5, 0.5, 0.5, 0.5, 0.5, 0.05; break;
				case Components::Three: m.head(nP()) << 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.05, 0.05; break;
			}
			if (m_finite)
				m.segment(nP() - 1, 1).setConstant(0.2);
			m.segment(nP(), nOffRes()).setConstant(0.1);
			m.tail(nPD()).setConstant(0.1);
			return m;
		}
		
		mcDESPOT(const Components &c, vector<shared_ptr<SignalFunctor>> &data,
				 const FieldStrength &tesla, const OffRes &offRes, const Scaling &s = Scaling::NormToMean,
				 const bool &isFinite = false, const bool &debug = false) :
			DESPOTFunctor(data, tesla, offRes, s, debug),
			m_components(c), m_finite(isFinite)
		{
			m_names.reserve(inputs());
			switch (c) {
				case Components::One: m_names = {"T1", "T2"}; break;
				case Components::Two: m_names = {"T1_a", "T2_a", "T1_b", "T2_b", "tau_a", "f_a"}; break;
				case Components::Three: m_names = {"T1_a", "T2_a", "T1_b", "T2_b", "T1_c", "T2_c", "tau_a", "f_a", "f_c"}; break;
			}
			if (m_finite)
				m_names.emplace_back("deltaf0_");
			for (size_t i = 0; i < nOffRes(); i++)
				m_names.emplace_back("f0_" + std::to_string(i));
			for (size_t i = 0; i < nPD(); i++)
				m_names.emplace_back("PD_" + std::to_string(i));
		}
		
		const bool constraint(const VectorXd &params) const override {
			// Negative PD or T1/T2 makes no sense
			if ((params[0] <= 0.) || (params[1] <= 0.))
				return false;
			
			if (m_components == Components::One) {
				return true;
			} else if (m_components == Components::Two) {
				// Check that T1_a, T2_a < T1_b, T2_b and that f_a makes sense
				if ((params[0] < params[2]) &&
					(params[1] < params[3]) &&
					(params[5] <= 1.0))
					return true;
				else
					return false;
			} else if (m_components == Components::Three) {
				// Check that T1/2_a < T1/2_b < T1/2_c and that f_a + f_c makes sense 
				if ((params[0] < params[2]) &&
					(params[1] < params[3]) &&
					(params[2] < params[4]) &&
					(params[3] < params[5]) &&
					((params[7] + params[8]) <= 1.0))
					return true;
				else
					return false;
			} else
				return true;
		}
		
		const ArrayXd theory(const Ref<VectorXd> &params) override {
			ArrayXd t(values());
			int index = 0;
			if (m_debug) cout << endl << __PRETTY_FUNCTION__ << endl << "Params: " << params.transpose() << endl;
			for (size_t i = 0; i < m_signals.size(); i++) {
				ArrayXd sig = m_signals.at(i)->signal(params.head(nP()), m_B1, params[nP()]);
				if (m_debug)
					cout << "Raw signal " << i << ": " << sig.transpose() << endl;
				switch (m_scaling) {
					case (Scaling::NormToMean) : sig /= sig.mean(); break;
					case (Scaling::Global)     : sig = sig * params(nP() + nOffRes()); break;
				}
				if (m_debug)
					cout << "Scaled signal " << i << ": " << sig.transpose() << endl;
				m_theory.at(i) = sig;
				t.segment(index, sig.rows()) = sig;
				index += sig.rows();
			}
			return t;
		}
};

#endif
