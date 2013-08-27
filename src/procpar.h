//
//  procpar.h
//  procparse
//
//  Created by Tobias Wood on 10/07/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#ifndef AGILENT_PROCPAR
#define AGILENT_PROCPAR

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include <map>
#include <algorithm>
#include <exception>

using namespace std;

namespace Agilent {
	
	class Parameter {
		public:
			enum class Type : size_t {
				Real = 1,
				String
			};
			
			enum class SubType : size_t {
				Junk = 0,
				Real = 1,
				String,
				Delay,
				Flag,
				Freq,
				Pulse,
				Int
			};
	
		protected:
			string m_name;
			Type m_type;
			SubType m_subtype;
			double  m_max, m_min, m_step;
			int  m_ggroup, m_dgroup, m_protection, m_active, m_intptr;
			vector<string> m_stringValues, m_stringAllowed;
			vector<double> m_realValues, m_realAllowed;
			
		public:
			Parameter();
			Parameter(const string &name, const SubType &st, const string &val);          //!< Construct a string parameter with one value
			Parameter(const string &name, const SubType &st, const double &val);          //!< Construct a string parameter with multiple values
			Parameter(const string &name, const SubType &st,
			          const vector<string> &vals, const vector<string> &allowed,
					  const int ggroup, const int dgroup,
					  const double max, const double min, const double step,
					  const int protection, const int active, const int intptr); //!< Construct a real parameter with one value
			Parameter(const string &name, const SubType &st,
			          const vector<double> &vals, const vector<double> &allowed,
					  const int ggroup, const int dgroup,
					  const double max, const double min, const double step,
					  const int protection, const int active, const int intptr); //!< Construct a real parameter with multiple values
			Parameter(const string &name, const SubType &st, const int n);                //!< Construct a blank parameter with space for several values

			const string &name() const;
			const size_t nvals() const;
			const size_t nallowed() const;
			const string &type_name() const;
			const string &subtype_name() const;
			
			const string &stringValue(const size_t i) const;
			const vector<string> &stringValues() const;
			const double &realValue(const size_t i) const;
			const vector<double> &realValues() const;
			const string print_values() const;
			const string print_allowed() const;
			const string print() const;
			
			const bool operator==(const Parameter &other);
			const bool operator!=(const Parameter &other);
			friend ostream& operator<<(ostream &os, const Parameter &p);
			friend istream& operator>>(istream &is, Parameter &p);
	};
	
	class ProcPar {
		protected:
			typedef map<string, Parameter> ParMap;
			ParMap m_parameters;
		
		public:
			friend ostream& operator<<(ostream &os, const ProcPar &p);
			friend istream& operator>>(istream &is, ProcPar &p);
			
			const bool contains(const string &name) const;
			void insert(const Parameter &p);
			void remove(const string &name);
			size_t count() const;
			
			const Parameter &parameter(const string &name) const;
			const vector<string> names() const;
			const double realValue(const string &name, const size_t index = 0) const;
			const vector<double> realValues(const string &name) const;
			const string &stringValue(const string &name, const size_t index = 0) const;
			const vector<string> stringValues(const string &name) const;
	};
} // End namespace Agilent

#endif
