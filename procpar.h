//
//  procpar.h
//  procparse
//
//  Created by Tobias Wood on 10/07/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#ifndef PROCPAR
#define PROCPAR

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

#define PROCPAR_ERROR( err ) do { std::cerr << __PRETTY_FUNCTION__ << ": " << ( err ) << std::flush << std::endl; } while(0)
#define PROCPAR_FAIL( err ) do { PROCPAR_ERROR( err ); exit(EXIT_FAILURE); } while(0)

namespace ProcPar {

	enum Type {
		TYPE_REAL = 1,
		TYPE_STRING
	};
	
	enum SubType {
		SUB_REAL = 1,
		SUB_STRING,
		SUB_DELAY,
		SUB_FLAG,
		SUB_FREQ,
		SUB_PULSE,
		SUB_INT
	};
	
	class Parameter {
		protected:
			string _name;
			int  _type, _subtype;
			double  _max, _min, _step;
			int  _ggroup, _dgroup, _protection, _active, _intptr;
			vector<string> _stringValues, _stringAllowed;
			vector<double> _realValues, _realAllowed;
		public:
			Parameter();
			Parameter(const string &name, const int &st, const string &val);          //!< Construct a string parameter with one value
			Parameter(const string &name, const int &st, const double &val);          //!< Construct a string parameter with multiple values
			Parameter(const string &name, const int &st,
			          const vector<string> &vals, const vector<string> &allowed,
					  const int ggroup, const int dgroup,
					  const double max, const double min, const double step,
					  const int protection, const int active, const int intptr); //!< Construct a real parameter with one value
			Parameter(const string &name, const int &st,
			          const vector<double> &vals, const vector<double> &allowed,
					  const int ggroup, const int dgroup,
					  const double max, const double min, const double step,
					  const int protection, const int active, const int intptr); //!< Construct a real parameter with multiple values
			Parameter(const string &name, const int &st, const int n);                //!< Construct a blank parameter with space for several values

			const string &name() const;
			const size_t nvals() const;
			const size_t nallowed() const;
			const string &type_name() const;
			const string &subtype_name() const;
			
			const string &stringValue(const int i) const;
			const vector<string> &stringValues() const;
			const double &realValue(const int i) const;
			const vector<double> &realValues() const;
			const string print_values() const;
			const string print_allowed() const;
			const string print() const;
			
			const bool operator==(const Parameter &other);
			friend ostream& operator<<(ostream &os, const Parameter &p);
			friend istream& operator>>(istream &is, const Parameter &p);
	};
	
	typedef map<string, Parameter> ParameterList;

	ParameterList ReadProcpar(const string &path);
	void PrintProcpar(const ParameterList &pp);

} // End namespace ProcPar

#endif
