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
		
		long inputs() const { return m_inputs; }
		long values() const { return m_values; }
		
		virtual int operator()(const VectorXd &params, ArrayXd &diffs) const = 0;
		
		virtual const ArrayXd theory(const VectorXd &params) const = 0;
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
		enum class PDMode {
			Normalise, Global, Individual
		};
		static const string to_string(const PDMode &p) {
			switch (p) {
				case PDMode::Normalise: return "normalise";
				case PDMode::Global: return "global";
				case PDMode::Individual: return "individual";
			}
		}
		enum class OffResMode {
			Map = 0, Single, Multi, Bounded, MultiBounded
		};
		
		static const size_t nPD(const PDMode& p, const size_t &nSignals) {
			switch (p) {
				case PDMode::Normalise: return 0;
				case PDMode::Global: return 1;
				case PDMode::Individual: return static_cast<size_t>(nSignals);
			}
		}
		
		static const size_t nOffRes(const OffResMode &m, const size_t &nSignals) {
			switch (m) {
				case OffResMode::Map: return 0;
				case OffResMode::Single: return 1;
				case OffResMode::Multi: return static_cast<size_t>(nSignals);
				case OffResMode::Bounded: return 1;
				case OffResMode::MultiBounded: return static_cast<size_t>(nSignals);
			}
		}
	
	protected:
		const OffResMode m_offRes;
		const PDMode m_PDMode;
		size_t m_nP, m_nV, m_nOffRes, m_nPD;
		vector<DESPOTData> &m_data;
		const bool m_debug;
	
	public:
		const long inputs() const { return m_nP + m_nOffRes + m_nPD; }
		const long values() const { return m_nV; }
		
		DESPOTFunctor(vector<DESPOTData> &data,
		              const OffResMode &offRes = OffResMode::Single, const PDMode &PD = PDMode::Normalise,
		              const bool &debug = false) :
			m_data(data), m_offRes(offRes), m_PDMode(PD), m_debug(debug)
		{
			m_nV = 0;
			for (int i = 0; i < data.size(); i++) {
				if (data[i].flip().size() != data[i].signal().size()) {
					cerr << "Angles and signals size mis-match for signal " << i << endl;
					cerr << "Angles = " << data[i].flip().size() << " signal = " << data[i].signal().size() << endl;
					exit(EXIT_FAILURE);
				}
				m_nV += m_data[i].flip().size();
			}
		}
		
		const ArrayXd signals() const {
			ArrayXd v(values());
			int index = 0;
			if (m_debug) cout << __PRETTY_FUNCTION__ << endl;
			for (int i = 0; i < m_data.size(); i++) {
				v.segment(index, m_data[i].signal().size()) = m_data[i].signal();
				index += m_data[i].signal().size();
			}
			return v;
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
				case Components::One: return "one";
				case Components::Two: return "two";
				case Components::Three: return "three";
			}
		};
		static const int nP(const Components &c) {
			switch (c) {
				case Components::One: return 2;
				case Components::Two: return 6;
				case Components::Three: return 9;
			}
		}
	protected:
		const Components m_components;
	
	public:
		static const vector<string> names(const Components &c) {
			switch (c) {
				case Components::One: return { "T1", "T2" };
				case Components::Two: return { "T1_a", "T2_a", "T1_b", "T2_b", "tau_a", "f_a"  };
				case Components::Three: return { "T1_a", "T2_a", "T1_b", "T2_b", "T1_c", "T2_c", "tau_a", "f_a", "f_c" };
			}
		}
		
		static const ArrayXXd defaultBounds(const Components &c, const FieldStrength &tesla) {
			ArrayXXd b(nP(c), 2);
			switch (tesla) {
				case FieldStrength::Three:
					switch (c) {
						case Components::One: b << 0.25, 3.0, 0.01, 0.25; break;
						case Components::Two: b << 0.25, 1.0, 0.01, 0.05, 0.75, 1.5, 0.01, 0.05, 0.01, 0.5, 0.001, 0.95; break;
						case Components::Three: b << 0.35, 0.55, 0.002, 0.016, 0.700, 2.0, 0.075, 0.145, 3.5, 7.5, 0.175, 0.5, 0.05, 0.3, 0.001, 0.3, 0.001, 0.95; break;
					}
					break;
				case FieldStrength::Seven:
					switch (c) {
						case Components::One: b << 0.25, 5.0, 0.01, 0.1; break;
						case Components::Two: b << 0.1, 0.5, 0.001, 0.025, 1.0, 2.5, 0.04, 0.08, 0.01, 0.25, 0.001, 1.0; break;
						case Components::Three: b << 0.1, 0.5, 0.001, 0.025, 1.0, 2.5, 0.04, 0.08, 3., 4.5, 0.5, 2.0, 0.01, 0.25, 0.001, 0.4, 0.001, 1.0; break;
					}
					break;
				case FieldStrength::Unknown:
					b.setZero(); break;
			}
			return b;
		}
	
		mcDESPOT(const Components &c, vector<DESPOTData> &data,
				 const OffResMode &offRes, const PDMode &PD = PDMode::Normalise,
				 const bool &debug = false) :
			DESPOTFunctor(data, offRes, PD, debug),
			m_components(c)
		{
			m_nP = nP(c);
			m_nOffRes = nOffRes(offRes, m_data.size());
			m_nPD = nPD(PD, m_data.size());
		}
		
		const bool constraint(const VectorXd &params) {
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
		
		const ArrayXd theory(const VectorXd &params) const {
			ArrayXd t(values());
			int index = 0;
			if (m_debug) cout << __PRETTY_FUNCTION__ << endl << "Params: " << params.transpose() << endl;
			for (int i = 0; i < m_data.size(); i++) {
				MagVector M(3, m_data[i].flip().size());
				if ((m_offRes == OffResMode::Single) || (m_offRes == OffResMode::Bounded))
					m_data[i].f0_off = params[m_nP];
				else if ((m_offRes == OffResMode::Multi) || (m_offRes == OffResMode::MultiBounded))
					m_data[i].f0_off = params[m_nP + i];
				double PD;
				switch (m_PDMode) {
					case (PDMode::Normalise): PD = 1.; break;
					case (PDMode::Global): PD = params[m_nP + m_nOffRes]; break;
					case (PDMode::Individual): PD = params[m_nP + m_nOffRes + i]; break;
				}
				if (m_data[i].spoil == true) {
					switch (m_components) {
						case Components::One: M = One_SPGR(m_data[i], params.head(m_nP), PD); break;
						case Components::Two: M = Two_SPGR(m_data[i], params.head(m_nP), PD); break;
						case Components::Three: M = Three_SPGR(m_data[i], params.head(m_nP), PD); break;
					}
				} else {
					switch (m_components) {
						case Components::One: M = One_SSFP(m_data[i], params.head(m_nP), PD); break;
						case Components::Two: M = Two_SSFP(m_data[i], params.head(m_nP), PD); break;
						case Components::Three: M = Three_SSFP(m_data[i], params.head(m_nP), PD); break;
					}
				}
				ArrayXd theory = SigMag(M);
				if (m_PDMode == PDMode::Normalise) {
					theory /= theory.mean();
				}
				t.segment(index, m_data[i].signal().size()) = theory;
				if (m_debug) cout << theory.transpose() << endl;
				index += m_data[i].signal().size();
			}
			return t;
		}
				
		int operator()(const VectorXd &params, ArrayXd &diffs) const {
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

class mcFinite : public mcDESPOT {
	public:
		static const int nP(const Components &c) {
			switch (c) {
				case Components::One: return 3;
				case Components::Two: return 7;
				case Components::Three: return 10;
			}
		}
		
		static const vector<string> names(const Components &c) {
			switch (c) {
				case Components::One: return { "T1", "T2", "delta_f" };
				case Components::Two: return { "T1_a", "T2_a", "T1_b", "T2_b", "tau_a", "f_a", "delta_f"  };
				case Components::Three: return { "T1_a", "T2_a", "T1_b", "T2_b", "T1_c", "T2_c", "tau_a", "f_a", "f_c", "delta_f" };
			}
		}
		
		static const ArrayXXd defaultBounds(const Components &c, const FieldStrength &tesla) {
			ArrayXXd b(nP(c), 2);
			switch (tesla) {
				case FieldStrength::Three:
					switch (c) {
						case Components::One: b << 0.25, 3.0, 0.01, 0.25, 0., 100.; break;
						case Components::Two: b << 0.25, 1.0, 0.01, 0.05, 0.75, 1.5, 0.01, 0.05, 0.01, 0.5, 0.001, 0.95, 0., 100.; break;
						case Components::Three: b << 0.35, 0.55, 0.002, 0.016, 0.700, 2.0, 0.075, 0.145, 3.5, 7.5, 0.175, 0.5, 0.05, 0.3, 0.001, 0.3, 0.001, 0.95, 0., 100.; break;
					}
					break;
				case FieldStrength::Seven:
					switch (c) {
						case Components::One: b << 0.25, 5.0, 0.01, 0.1, 0., 100.; break;
						case Components::Two: b << 0.1, 0.5, 0.001, 0.025, 1.0, 2.5, 0.04, 0.08, 0.01, 0.25, 0.001, 1.0, 0., 100.; break;
						case Components::Three: b << 0.1, 0.5, 0.001, 0.025, 1.0, 2.5, 0.04, 0.08, 3., 4.5, 0.5, 2.0, 0.01, 0.25, 0.001, 0.4, 0.001, 1.0, 0., 100.; break;
					}
					break;
				case FieldStrength::Unknown:
					b.setZero(); break;
			}
			return b;
		}
		
		mcFinite(const Components &c, vector<DESPOTData> &data,
				 const OffResMode &offRes, const PDMode &PD = PDMode::Normalise,
				 const bool &debug = false) :
			mcDESPOT(c, data, offRes, PD, debug)
		{
			m_nP = nP(c);
			m_nOffRes = nOffRes(offRes, m_data.size());
			m_nPD = nPD(PD, m_data.size());
		}
		
		const ArrayXd theory(const VectorXd &params) const {
			ArrayXd t(values());
			int index = 0;
			if (m_debug) cout << __PRETTY_FUNCTION__ << endl << "Params: " << params.transpose() << endl;
			for (int i = 0; i < m_data.size(); i++) {
				MagVector M(3, m_data[i].flip().size());
				if ((m_offRes == OffResMode::Single) || (m_offRes == OffResMode::Bounded))
					m_data[i].f0_off = params[m_nP];
				else if ((m_offRes == OffResMode::Multi) || (m_offRes == OffResMode::MultiBounded))
					m_data[i].f0_off = params[m_nP + i];
				double PD;
				switch (m_PDMode) {
					case (PDMode::Normalise): PD = 1.; break;
					case (PDMode::Global): PD = params[m_nP + m_nOffRes]; break;
					case (PDMode::Individual): PD = params[m_nP + m_nOffRes + i]; break;
				}
				switch (m_components) {
					case Components::One: M = One_SSFP_Finite(m_data[i], params.head(m_nP), PD); break;
					case Components::Two: M = Two_SSFP_Finite(m_data[i], params.head(m_nP), PD); break;
					case Components::Three: M = Three_SSFP_Finite(m_data[i], params.head(m_nP), PD); break;
				}
				ArrayXd theory = SigMag(M);
				if (m_PDMode == PDMode::Normalise) {
					theory /= theory.mean();
				}
				t.segment(index, m_data[i].signal().size()) = theory;
				index += m_data[i].signal().size();
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
	public:
		static const int nP() {
			return 1;
		}
		
		static const vector<string> &names() {
			static const vector<string> _names = { { "FM_T2" } };
			return _names;
		}
		
		static const ArrayXXd defaultBounds(const FieldStrength tesla) {
			ArrayXXd b(nP(), 2);
			switch (tesla) {
				case FieldStrength::Three: b << 0.010, 0.5; return b;
				case FieldStrength::Seven: b << 0.005, 0.25; return b;
				case FieldStrength::Unknown: b.setZero(); return b;
			}
		}
		
		const bool constraint(const VectorXd &params) {
			if (params[0] < 0.)
				return false;
			else
				return true;
		}
		
		DESPOT2FM(vector<DESPOTData> &data, const double T1,
				  const OffResMode &offRes, const PDMode &PD = PDMode::Global,
				  const bool &debug = false) :
			DESPOTFunctor(data, offRes, PD, debug), m_T1(T1)
		{
			m_nP = nP();
			m_nOffRes = nOffRes(offRes, m_data.size());
			m_nPD = nPD(PD, m_data.size());
		}
		DESPOT2FM(const DESPOT2FM &) = delete;
		
		const ArrayXd theory(const VectorXd &params) const {
			VectorXd T1T2(2);
			T1T2 << m_T1, params[0];
			
			ArrayXd t(values());
			int index = 0;
			for (int i = 0; i < m_data.size(); i++) {
				MagVector M(3, m_data[i].flip().size());
				if ((m_offRes == OffResMode::Single) || (m_offRes == OffResMode::Bounded))
					m_data[i].f0_off = params[m_nP];
				else if ((m_offRes == OffResMode::Multi) || (m_offRes == OffResMode::MultiBounded))
					m_data[i].f0_off = params[m_nP + i];
				double PD;
				switch (m_PDMode) {
					case (PDMode::Normalise): PD = 1.; break;
					case (PDMode::Global): PD = params[m_nP + m_nOffRes]; break;
					case (PDMode::Individual): PD = params[m_nP + m_nOffRes + i]; break;
				}
				M = One_SSFP(m_data.at(i), T1T2, PD);
				ArrayXd theory = SigMag(M);
				if (m_PDMode == PDMode::Normalise) {
					theory /= theory.mean();
				}
				t.segment(index, m_data[i].flip().size()) = theory;
				index += m_data[i].flip().size();
			}
			return t;
		}
				
		int operator()(const VectorXd &params, ArrayXd &diffs) const {
			eigen_assert(diffs.size() == values());
			ArrayXd t = theory(params);
			ArrayXd s = signals();
			diffs = t - s;
			return 0;
		}
};

#endif
