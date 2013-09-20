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
#include <Eigen/Dense>

#include "DESPOT.h"

using namespace std;
using namespace Eigen;

//******************************************************************************
#pragma mark Functor Base Class
//******************************************************************************
// From Nonlinear Tests in Eigen 
template<typename _Scalar, int NX=Dynamic, int NY=Dynamic>
class Functor
{
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
		
		virtual int operator()(const VectorXd &params, Ref<ArrayXd> diffs) = 0;
		
		virtual const ArrayXd theory(const VectorXd &params) = 0;
		virtual const ArrayXd signals() const = 0;
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
			Global, PerSignal, MeanPerSignal, MeanPerType
		};
		static const string to_string(const Scaling &p) {
			switch (p) {
				case Scaling::Global: return "global";
				case Scaling::PerSignal: return "per signal";
				case Scaling::MeanPerSignal: return "normalised to per signal mean";
				case Scaling::MeanPerType: return "normalised to mean for each signal type (SPGR/SSFP)";
			}
		}
		enum class OffResMode {
			Map = 0, Single, Multi, Bounded, MultiBounded
		};
	
	protected:
		const FieldStrength m_fieldStrength;
		const OffResMode m_offRes;
		const Scaling m_scaling;
		size_t m_nV;
		vector<Info> m_info;
		vector<ArrayXd> m_signals;
		vector<string> m_names; // Subclasses responsible for initialising this
		const bool m_debug;
	
		ArrayXXd offResBounds() {
			ArrayXXd b(nOffRes(), 2);
			bool symmetricB0 = true;
			for (auto &i : m_info) {
				if (fmod(i.phase, M_PI) > numeric_limits<double>::epsilon()) {
					symmetricB0 = false;
				}
			}
			for (size_t i = 0; i < nOffRes(); i++) {
				if (symmetricB0) b(i, 0) = 0;
				else b(i, 0) = -0.5 / m_info.at(i).TR;
				b(i, 1) = 0.5 / m_info.at(i).TR;
			}
			return b;
		}
	
		ArrayXXd PDBounds() {
			ArrayXXd b(nPD(), 2);
			for (size_t i = 0; i < nPD(); i++) {
				b(i, 0) = 1.e4;
				b(i, 1) = 5.e6;
			}
			return b;
		}
		
		virtual const ArrayXXd defaultBounds() = 0;
		
	public:
		const size_t nOffRes() const {
			switch (m_offRes) {
				case OffResMode::Map: return 0;
				case OffResMode::Single: return 1;
				case OffResMode::Multi: return m_info.size();
				case OffResMode::Bounded: return 1;
				case OffResMode::MultiBounded: return m_info.size();
			}
		}
		
		const size_t nPD() const {
			switch (m_scaling) {
				case Scaling::Global:        return 1;
				case Scaling::PerSignal:     return m_info.size();
				case Scaling::MeanPerSignal: return 0;
				case Scaling::MeanPerType:   return 0;
			}
		}
		
		virtual const size_t nP() const = 0;
		virtual const bool constraint(const VectorXd &params) const = 0;
		const long inputs() const override { return this->nP() + nOffRes() + nPD(); }
		const long values() const override { return m_nV; }
		
		DESPOTFunctor(vector<Info> &info_in,
					  const FieldStrength &tesla,
					  const OffResMode &offRes = OffResMode::Single,
					  const Scaling &s = Scaling::MeanPerSignal,
		              const bool &debug = false) :
			m_info(info_in),
			m_fieldStrength(tesla), m_offRes(offRes),
			m_scaling(s), m_debug(debug)
		{
			m_signals.reserve(m_info.size());
			m_nV = 0;
			for (size_t i = 0; i < m_info.size(); i++) {
				m_signals.emplace_back(m_info.at(i).nAngles());
				m_nV += m_info.at(i).nAngles();
			}
		}
		
		const vector<string> &names() { return m_names; }
		Info &info(const size_t i) { return m_info.at(i); }
		ArrayXd &signal(const size_t s) { return m_signals.at(s); }
		const ArrayXd signals() const {
			ArrayXd v(values());
			int index = 0;
			if (m_debug) cout << __PRETTY_FUNCTION__ << endl;
			for (size_t i = 0; i < m_signals.size(); i++) {
				v.segment(index, m_signals.at(i).size()) = m_signals.at(i);
				index += m_signals.at(i).size();
			}
			return v;
		}
		
		int operator()(const VectorXd &params, Ref<ArrayXd> diffs) override {
			eigen_assert(diffs.size() == values());
			ArrayXd t = theory(params);
			ArrayXd s = signals();
			diffs = t - s;
			if (m_debug) {
				cout << endl << __PRETTY_FUNCTION__ << endl;
				cout << "Diffs:  " << diffs.transpose() << endl;
				cout << "Sum:    " << diffs.square().sum() << endl;
			}
			return 0;
		}
};

//******************************************************************************
#pragma mark mcDESPOT Functor
//******************************************************************************
class mcDESPOT : public DESPOTFunctor {
	public:
		enum class Components {
			One, Two, Three
		};
		static const string to_string(const Components& c) {
			switch (c) {
				case Components::One: return "1";
				case Components::Two: return "2";
				case Components::Three: return "3";
			}
		};
	protected:
		const Components m_components;
	
	public:
		const size_t nP() const override {
			switch (m_components) {
				case Components::One: return 2;
				case Components::Two: return 6;
				case Components::Three: return 9;
			}
		}
		
		const ArrayXXd defaultBounds() override {
			ArrayXXd b(inputs(), 2);
			switch (m_fieldStrength) {
				case FieldStrength::Three:
					switch (m_components) {
						case Components::One:   b.block(0, 0, 2, 2) << 0.1, 4.0, 0.01, 1.5; break;
						case Components::Two:   b.block(0, 0, 6, 2) << 0.1, 0.25, 0.002, 0.03, 0.7, 4.0, 0.075, 0.145, 0.05, 0.3, 0.0, 0.95; break;
						case Components::Three: b.block(0, 0, 9, 2) << 0.1, 0.25, 0.002, 0.03, 0.7, 2.0, 0.075, 0.145, 3.5, 4.0, 0.8, 1.5, 0.05, 0.3, 0.001, 0.3, 0.0, 0.95; break;
					} break;
				case FieldStrength::Seven:
					switch (m_components) {
						case Components::One:   b.block(0, 0, 2, 2) << 0.1, 4.0, 0.01, 2.0; break;
						case Components::Two:   b.block(0, 0, 6, 2) << 0.1, 0.5, 0.001, 0.025, 1.0, 4.0, 0.04, 0.08, 0.01, 0.25, 0.001, 1.0; break;
						case Components::Three: b.block(0, 0, 9, 2) << 0.1, 0.5, 0.001, 0.025, 1.0, 2.5, 0.04, 0.08, 3., 4.5, 0.5, 2.0, 0.01, 0.25, 0.001, 0.4, 0.001, 1.0; break;
					} break;
				case FieldStrength::Unknown:
					switch (m_components) {
						case Components::One:   b.block(0, 0, 2, 2).setZero(); break;
						case Components::Two:   b.block(0, 0, 6, 2).setZero(); break;
						case Components::Three: b.block(0, 0, 9, 2).setZero(); break;
					} break;
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
			m.segment(nP(), nOffRes()).setConstant(0.1);
			m.tail(nPD()).setConstant(0.1);
			return m;
		}
		
		mcDESPOT(const Components &c, vector<Info> &data,
				 const FieldStrength &tesla, const OffResMode &offRes, const Scaling &s = Scaling::MeanPerSignal,
				 const bool &debug = false) :
			DESPOTFunctor(data, tesla, offRes, s, debug),
			m_components(c)
		{
			m_names.reserve(inputs());
			switch (c) {
				case Components::One: m_names = {"T1", "T2"}; break;
				case Components::Two: m_names = {"T1_a", "T2_a", "T1_b", "T2_b", "tau_a", "f_a"}; break;
				case Components::Three: m_names = {"T1_a", "T2_a", "T1_b", "T2_b", "T1_c", "T2_c", "tau_a", "f_a", "f_c"}; break;
			}
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
		
		const ArrayXd theory(const VectorXd &params) override {
			ArrayXd t(values());
			int index = 0;
			if (m_debug) cout << __PRETTY_FUNCTION__ << endl << "Params: " << params.transpose() << endl;
			for (size_t i = 0; i < m_info.size(); i++) {
				MagVector M(3, m_info[i].flip().size());
				if ((m_offRes == OffResMode::Single) || (m_offRes == OffResMode::Bounded))
					m_info.at(i).f0 = params[nP()];
				else if ((m_offRes == OffResMode::Multi) || (m_offRes == OffResMode::MultiBounded))
					m_info.at(i).f0 = params[nP() + i];
				double PD;
				switch (m_scaling) {
					case (Scaling::Global):        PD = params[nP() + nOffRes()]; break;
					case (Scaling::PerSignal):     PD = params[nP() + nOffRes() + i]; break;
					case (Scaling::MeanPerSignal): PD = 1; break;
					case (Scaling::MeanPerType):   PD = 1; break;
				}
				if (m_info[i].spoil == true) {
					switch (m_components) {
						case Components::One: M = One_SPGR(m_info[i], params.head(nP()), PD); break;
						case Components::Two: M = Two_SPGR(m_info[i], params.head(nP()), PD); break;
						case Components::Three: M = Three_SPGR(m_info[i], params.head(nP()), PD); break;
					}
				} else {
					switch (m_components) {
						case Components::One: M = One_SSFP(m_info[i], params.head(nP()), PD); break;
						case Components::Two: M = Two_SSFP(m_info[i], params.head(nP()), PD); break;
						case Components::Three: M = Three_SSFP(m_info[i], params.head(nP()), PD); break;
					}
				}
				ArrayXd theory = SigMag(M);
				if (m_scaling == Scaling::MeanPerSignal) {
					theory /= theory.mean();
				}
				t.segment(index, m_signals.at(i).size()) = theory;
				if (m_debug) cout << theory.transpose() << endl;
				index += m_signals.at(i).size();
			}
			return t;
		}
};

class mcFinite : public mcDESPOT {
	public:
		const size_t nP() const override {
			return mcDESPOT::nP() + 1;
		}
				
		const ArrayXXd defaultBounds() override {
			ArrayXXd b(inputs(), 2);
			switch (m_fieldStrength) {
				case FieldStrength::Three:
					switch (m_components) {
						case Components::One:   b.block(0, 0, 3, 2) << 0.1, 4.0, 0.01, 1.5, 0., 100.; break;
						case Components::Two:   b.block(0, 0, 7, 2) << 0.1, 0.25, 0.002, 0.03, 0.7, 4.0, 0.075, 0.145, 0.05, 0.3, 0.0, 0.95, 0., 100.; break;
						case Components::Three: b.block(0, 0, 10, 2) << 0.1, 0.25, 0.002, 0.03, 0.7, 2.0, 0.075, 0.145, 3.5, 4.0, 0.8, 1.5, 0.05, 0.3, 0.001, 0.3, 0.001, 0.95, 0., 100.; break;
					} break;
				case FieldStrength::Seven:
					switch (m_components) {
						case Components::One:   b.block(0, 0, 3, 2) << 0.1, 4.0, 0.01, 2.0, 0., 100.; break;
						case Components::Two:   b.block(0, 0, 7, 2) << 0.1, 0.5, 0.001, 0.025, 1.0, 2.5, 0.04, 0.08, 0.01, 0.25, 0.001, 1.0, 0., 100.; break;
						case Components::Three: b.block(0, 0, 10, 2) << 0.1, 0.5, 0.001, 0.025, 1.0, 2.5, 0.04, 0.08, 3., 4.5, 0.5, 2.0, 0.01, 0.25, 0.001, 0.4, 0.001, 1.0, 0., 100.; break;
					} break;
				case FieldStrength::Unknown:
					switch (m_components) {
						case Components::One:   b.block(0, 0, 3, 2).setZero(); break;
						case Components::Two:   b.block(0, 0, 7, 2).setZero(); break;
						case Components::Three: b.block(0, 0, 10, 2).setZero(); break;
					} break;
			}
			b.block(nP(), 0, nOffRes(), 2) = offResBounds();
			b.block(nP() + nOffRes(), 0, nPD(), 2) = PDBounds();
			return b;
		}
		
		const ArrayXd defaultThresholds() {
			ArrayXd m(inputs());
			switch (m_components) {
				case Components::One: m.head(nP()) << 0.05, 0.05, 1.0; break;
				case Components::Two: m.head(nP()) << 0.5, 0.5, 0.5, 0.5, 0.5, 0.05, 1.0; break;
				case Components::Three: m.head(nP()) << 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 0.5, 0.05, 0.05, 1.0; break;
			}
			m.segment(nP(), nOffRes()).setConstant(0.1);
			m.tail(nPD()).setConstant(0.1);
			return m;
		}
		
		mcFinite(const Components &c, vector<Info> &data,
				 const FieldStrength &tesla, const OffResMode &offRes, const Scaling &s = Scaling::MeanPerSignal,
				 const bool &debug = false) :
			mcDESPOT(c, data, tesla, offRes, s, debug)
		{
			m_names.insert(m_names.begin() + nP() - 1, "delta_f");
		}
		
		const ArrayXd theory(const VectorXd &params) override {
			ArrayXd t(values());
			int index = 0;
			if (m_debug) cout << __PRETTY_FUNCTION__ << endl << "Params: " << params.transpose() << endl;
			for (size_t i = 0; i < m_info.size(); i++) {
				MagVector M(3, m_info[i].flip().size());
				if ((m_offRes == OffResMode::Single) || (m_offRes == OffResMode::Bounded))
					m_info[i].f0 = params[nP()];
				else if ((m_offRes == OffResMode::Multi) || (m_offRes == OffResMode::MultiBounded))
					m_info[i].f0 = params[nP() + i];
				double PD;
				switch (m_scaling) {
					case (Scaling::Global):        PD = params[nP() + nOffRes()]; break;
					case (Scaling::PerSignal):     PD = params[nP() + nOffRes() + i]; break;
					case (Scaling::MeanPerSignal): PD = 1; break;
					case (Scaling::MeanPerType):   PD = 1; break;
				}
				switch (m_components) {
					case Components::One: M = One_SSFP_Finite(m_info[i], params.head(nP()), PD); break;
					case Components::Two: M = Two_SSFP_Finite(m_info[i], params.head(nP()), PD); break;
					case Components::Three: M = Three_SSFP_Finite(m_info[i], params.head(nP()), PD); break;
				}
				ArrayXd theory = SigMag(M);
				if (m_scaling == Scaling::MeanPerSignal) {
					theory /= theory.mean();
				}
				t.segment(index, m_signals.at(i).size()) = theory;
				index += m_signals.at(i).size();
				if (m_debug) cout << theory.transpose() << endl;
			}
			return t;
		}
};
//******************************************************************************
#pragma mark DESPOT2FM Functor
//******************************************************************************
class DESPOT2FM : public DESPOTFunctor {
	protected:
		double m_T1;
		bool m_finite;
		
	public:
		const size_t nP() const override {
			return 1;
		}
		
		const ArrayXXd defaultBounds() override {
			ArrayXXd b(inputs(), 2);
			switch (m_fieldStrength) {
				case FieldStrength::Three: b.block(0, 0, 1, 2) << 0.010, 1.5; break;
				case FieldStrength::Seven: b.block(0, 0, 1, 2) << 0.005, 2.0; break;
				case FieldStrength::Unknown: b.block(0, 0, 1, 2).setZero(); break;
			}
			b.block(nP(), 0, nOffRes(), 2) = offResBounds();
			b.block(nP() + nOffRes(), 0, nPD(), 2) = PDBounds();
			return b;
		}
		
		const ArrayXd defaultThresholds() {
			ArrayXd m(inputs());
			m.head(nP()) << 0.05, 0.05;
			m.segment(nP(), nOffRes()).setConstant(0.1);
			m.tail(nPD()).setConstant(0.1);
			return m;
		}
		
		const bool constraint(const VectorXd &params) const override {
			if (params[0] < 0.)
				return false;
			else
				return true;
		}
		
		DESPOT2FM(vector<Info> &data, const double T1,
				  const FieldStrength& tesla, const OffResMode &offRes, const Scaling &s = Scaling::MeanPerSignal,
				  const bool &finite = false, const bool &debug = false) :
			DESPOTFunctor(data, tesla, offRes, s, debug), m_T1(T1), m_finite(finite)
		{
			m_names.resize(inputs());
			m_names.at(0) = "T2";
			for (size_t i = 0; i < nOffRes(); i++)
				m_names.at(1 + i) = "f0_" + std::to_string(i);
			for (size_t i = 0; i < nPD(); i++)
				m_names.at(1 + nOffRes() + i) = "PD_" + std::to_string(i);
		}
		
		void setT1(const double T1) { m_T1 = T1; }
		
		const ArrayXd theory(const VectorXd &params) override {
			VectorXd T1T2(2);
			T1T2 << m_T1, params[0];
			
			ArrayXd t(values());
			int index = 0;
			for (size_t i = 0; i < m_info.size(); i++) {
				MagVector M(3, m_info[i].flip().size());
				if ((m_offRes == OffResMode::Single) || (m_offRes == OffResMode::Bounded))
					m_info[i].f0 = params[nP()];
				else if ((m_offRes == OffResMode::Multi) || (m_offRes == OffResMode::MultiBounded))
					m_info[i].f0 = params[nP() + i];
				double PD;
				switch (m_scaling) {
					case (Scaling::Global):        PD = params[nP() + nOffRes()]; break;
					case (Scaling::PerSignal):     PD = params[nP() + nOffRes() + i]; break;
					case (Scaling::MeanPerSignal): PD = 1; break;
					case (Scaling::MeanPerType):   PD = 1; break;
				}
				if (m_finite) {
					M = One_SSFP_Finite(m_info.at(i), T1T2, PD);
				} else {
					M = One_SSFP(m_info.at(i), T1T2, PD);
				}
				ArrayXd theory = SigMag(M);
				if (m_scaling == Scaling::MeanPerSignal) {
					theory /= theory.mean();
				}
				t.segment(index, m_info[i].flip().size()) = theory;
				index += m_info[i].flip().size();
			}
			return t;
		}
};

#endif
