//
//  procpar.c
//  procparse
//
//  Created by Tobias Wood on 10/07/2012.
//  Copyright (c) 2012 Tobias Wood. All rights reserved.
//

#include "procpar.h"

namespace Recon {

//******************************************************************************
#pragma mark Constructors
//******************************************************************************
Parameter::Parameter()
{


}

Parameter::Parameter(const string &name, const int &st, const string &val)
{
	switch (st) {
		case SUB_REAL: case SUB_DELAY: case SUB_FREQ: case SUB_PULSE: case SUB_INT:
			PROCPAR_FAIL("Tried to create string parameter " + _name + " with a real subtype");
			break;
		case SUB_STRING: case SUB_FLAG:
			break;
		case 0:
			// It's the stupid "junk" parameter
			break;
		default:
			PROCPAR_FAIL("Unknown subtype.");
	}
	_name = name;
	_type = TYPE_STRING;
	_subtype = st;
	_stringValues.push_back(val);
}

// The default value of 8 for max is taken from observations of procpar
Parameter::Parameter(const string &name, const int &st,
					 const vector<string> &vals, const vector<string> &allowed = vector<string>(),
					 const int ggroup = 0, const int dgroup = 0,
					 const double max = 8,
					 const double min = 0,
					 const double step = 0,
					 const int protection = 0, const int active = 0, const int intptr = 0) :
					 _ggroup(ggroup), _dgroup(dgroup), _protection(protection),
				     _active(active), _intptr(intptr), _max(max), _min(min), _step(step)
{
	switch (st) {
		case SUB_REAL: case SUB_DELAY: case SUB_FREQ: case SUB_PULSE: case SUB_INT:
			PROCPAR_FAIL("Tried to create string parameter " + _name + " with a real subtype");
			break;
		case SUB_STRING: case SUB_FLAG:
			break;
		case 0:
			// Junk subtype
			break;
		default:
			PROCPAR_FAIL("Unknown subtype.");
	}
	_name = name;
	_subtype = st;
	_type = TYPE_STRING;
	_stringValues = vals;
	_stringAllowed = allowed;
}

Parameter::Parameter(const string &name, const int &st, const double &val)
{
	switch (st) {
		case SUB_REAL: case SUB_DELAY: case SUB_FREQ: case SUB_PULSE: case SUB_INT:
			break;
		case SUB_STRING: case SUB_FLAG:
			PROCPAR_FAIL("Tried to create real parameter " + _name + " with a string subtype");
			break;
		case 0:
			break;
		default:
			PROCPAR_FAIL("Unknown subtype.");
	}
	_name = name;
	_type = TYPE_REAL;
	_subtype = st;
	_realValues.push_back(val);
}

Parameter::Parameter(const string &name, const int &st,
                     const vector<double> &vals, const vector<double> &allowed = vector<double>(),
					 const int ggroup = 0, const int dgroup = 0,
					 const double max = numeric_limits<double>::max(),
					 const double min = -numeric_limits<double>::max(),
					 const double step = 0,
					 const int protection = 0, const int active = 0, const int intptr = 0) :
					 _ggroup(ggroup), _dgroup(dgroup), _protection(protection),
				     _active(active), _intptr(intptr), _max(max), _min(min), _step(step)
{
	switch (st) {
		case SUB_REAL: case SUB_DELAY: case SUB_FREQ: case SUB_PULSE: case SUB_INT:
			break;
		case SUB_STRING: case SUB_FLAG:
			PROCPAR_FAIL("Tried to create real parameter " + _name + " with a string subtype");
			break;
		case 0:
			// Junk subtype
			break;
		default:
			PROCPAR_FAIL("Unknown subtype.");
	}
	_name = name;
	_type = TYPE_REAL;
	_subtype = st;
	_realValues = vals;
	_realAllowed = allowed;
}

Parameter::Parameter(const string &name, const int &st, const int n)
{
	_name = name;
	_subtype = st;
	switch (st) {
		case SUB_REAL: case SUB_DELAY: case SUB_FREQ: case SUB_PULSE: case SUB_INT:
			_type = TYPE_STRING;
			_stringValues.resize(n);
			break;
		case SUB_STRING: case SUB_FLAG:
			_type = TYPE_REAL;
			_realValues.resize(n);
			break;
		case 0:
			// Junk subtype, appears to be a string
			_type = TYPE_STRING;
			_stringValues.resize(n);
			break;
		default:
			PROCPAR_FAIL("Unknown subtype.");
	}
}

//******************************************************************************
#pragma mark Property accessors
//******************************************************************************
const string &Parameter::name() const { return _name; }

const size_t Parameter::nvals() const {
	if (_type == TYPE_STRING)
		return _stringValues.size();
	else
		return _realValues.size();
}

const size_t Parameter::nallowed() const {
	if (_type == TYPE_STRING)
		return _stringAllowed.size();
	else
		return _realAllowed.size();
}

const string &Parameter::stringValue(const int i = 0) const {
	if (_type == TYPE_STRING) {
		if (i < _stringValues.size())
			return _stringValues[i];
		else {
			stringstream ss;
			ss << "Tried to access element " << i << " of string parameter " << _name << " with only " << _stringValues.size() << " values.";
			PROCPAR_FAIL(ss.str());
		}
	}
	else
		PROCPAR_FAIL("Tried to read a string value from a real parameter.");
}

const double &Parameter::realValue(const int i = 0) const {
	if (_type == TYPE_REAL) {
		if (i < _realValues.size())
			return _realValues[i];
		else {
			stringstream ss;
			ss << "Tried to access element " << i << " of real parameter " << _name << " with only " << _realValues.size() << " values.";
			PROCPAR_FAIL(ss.str());
		}
	}
	else
		PROCPAR_FAIL("Tried to read a real value from string parameter " + _name + ".");
}

const vector<string> &Parameter::stringValues() const {
	if (_type == TYPE_STRING)
		return _stringValues;
	else
		PROCPAR_FAIL("Tried to read string values from real parameter " + _name + ".");
}

const vector<double> &Parameter::realValues() const {
	if (_type == TYPE_REAL)
		return _realValues;
	else
		PROCPAR_FAIL("Tried to read real values from string parameter " + _name + ".");
}

const string &Parameter::type_name() const
{
	static const string types[3] = { "", "Real", "String" };
	return types[_type];
}

const string &Parameter::subtype_name() const
{
	static const string subtypes[8] = { "", "Real", "String", "Delay",
									  "Flag", "Frequency", "Pulse", "Integer" };
	return subtypes[_subtype];
}

//******************************************************************************
#pragma mark Output
//******************************************************************************
const string Parameter::print_values() const {
	stringstream s;
	
	if (_type == TYPE_REAL) {
		for (vector<double>::const_iterator it = _realValues.begin(); it < _realValues.end(); it++)
			s << *it << " ";
	} else {
		vector<string>::const_iterator it = _stringValues.begin();
		s << "\"" << *it << "\"";
		for (it++; it < _stringValues.end(); it++)
			s << endl << "\"" << *it << "\"";
	}
	return s.str();
}

const string Parameter::print_allowed() const {
	stringstream s;
	
	if (_type == TYPE_REAL) {
		for (vector<double>::const_iterator it = _realAllowed.begin(); it < _realAllowed.end(); it++)
			s << *it;
	} else {
		for (vector<string>::const_iterator it = _stringAllowed.begin(); it < _stringAllowed.end(); it++)
			s << "\"" << *it << "\" ";
	}
	return s.str();
}

const string Parameter::print() const {
	stringstream s;

	s << _name << " " << _subtype << " " << _type << " "
	  << _max << " " << _min << " " << _step << " "
	  << _ggroup << " " << _dgroup << " " << _protection << " " << _active << " " << _intptr << endl;
	if (_type == TYPE_REAL) {
		s << _realValues.size() << " " << print_values() << endl
		  << _realAllowed.size() << " " << print_allowed();
	} else {
		s << _stringValues.size() << " " << print_values() << endl
		  << _stringAllowed.size() << " " << print_allowed();
	}
	return s.str();
}

ostream& operator<<(ostream &os, const Parameter &p) {
	os << p.print();
	return os;
}

//******************************************************************************
#pragma mark Input
//******************************************************************************
string readQuotedString(istream & is) {
	stringstream dummy, quoted;
	if (!is.get(*dummy.rdbuf(), '\"')) {
		PROCPAR_ERROR("Could not find an opening quote.");
		return "";
	}
	is.get(); // Consume the quote, get(...,delim) does not consume it
	if (is.peek() == '\"') {
		// We have an empty string
		is.get(); // Consume closing quote
		return "";
	} else if (!is.get(*quoted.rdbuf(), '\"')) {
		PROCPAR_ERROR("Could not find a closing quote.");
	} else
		is.get(); // Consume the quote
	return quoted.str();
}

istream& operator>>(istream &is, Parameter &p) {
	string name;
	int subtype, type, ggroup, dgroup, protection, active, intptr, nvals, nallowed;
	double max, min, step;
	vector<double> realVals, realAllowed;
	vector<string> stringVals, stringAllowed;
	if (is >> name >> subtype >> type
	       >> max >> min >> step
		   >> ggroup >> dgroup >> protection >> active >> intptr)
	{
		if (!(is >> nvals))
			PROCPAR_FAIL("Failed while reading number of values for parameter " + name);
		if (type == TYPE_REAL) {
			realVals.resize(nvals);
			for (int i = 0; i < nvals; i++)
				is >> realVals[i];
		} else if (type == TYPE_STRING) {
			stringVals.resize(nvals);
			for (int i = 0; i < nvals; i++)
				stringVals[i] = readQuotedString(is);
		} else
			PROCPAR_ERROR("Invalid type value for parameter " + name + ", no values read");
		if (!is)
			PROCPAR_FAIL("Failed while reading allowed values for parameter " + name + " from procpar file");
		if (!(is >> nallowed))
			PROCPAR_FAIL("Failed while reading number of allowed values for parameter " + name);
		if (type == TYPE_REAL) {
			realAllowed.resize(nallowed);
			for (int i = 0; i < nallowed; i++)
				is >> realAllowed[i];
		} else if (type == TYPE_STRING) {
			stringAllowed.resize(nallowed);
			for (int i = 0; i < nallowed; i++)
				stringAllowed[i] = readQuotedString(is);
		}
		if (!is)
			PROCPAR_FAIL("Failed while reading allowed values for parameter " + name + " from procpar file");
		
		if (type == TYPE_REAL)
			p = Parameter(name, subtype, realVals, realAllowed,
			              ggroup, dgroup, max, min, step, protection, active, intptr);
		else if (type == TYPE_STRING)
			p = Parameter(name, subtype, stringVals, stringAllowed,
			              ggroup, dgroup, max, min, step, protection, active, intptr);
		return is;
	} else {
		// An eof here is probably not an error, as procpar files have a trailing line
		// Anything else is an error.
		if (is.eof())
			return is;
		else
			PROCPAR_FAIL("Error while reading parameter definition line.");
	}
}

//******************************************************************************
#pragma mark Equality
//******************************************************************************
const bool Parameter::operator==(const Parameter &other)
{
	if (_type != other._type)
		return false;
	if (_subtype != other._subtype)
		return false;
	
	if (_type == TYPE_STRING) {
		if (_stringValues.size() != other._stringValues.size())
			return false;
		for (int i = 0; i < _stringValues.size(); i++)
			if (_stringValues[i] != other._stringValues[i]) return false;
	} else if (_type == TYPE_REAL) {
		if (_realValues.size() != other._realValues.size())
			return false;
		for (int i = 0; i < _realValues.size(); i++)
			if (_realValues[i] != other._realValues[i]) return false;
	}
	return true;
}

const bool Parameter::operator!=(const Parameter &other) { return !(operator==(other)); }
//******************************************************************************
#pragma mark Convenience
//******************************************************************************
/*
 *  Returns false if it can't find a procpar file, otherwise reads it into pp
 *
 */
bool ReadProcpar(const string &path, ParameterList &pp) {
	ifstream fpp;
	fpp.open(path.c_str());
	if (!fpp)
		return false;
	
	pp.clear();
	Parameter p;
	while(fpp >> p)
		pp.insert(pair<string, Parameter>(p.name(), p));
	
	fpp.close();
	return true;
}

bool WriteProcpar(const string &path, ParameterList &pp) {
	ofstream fpp;
	fpp.open(path.c_str());
	if (!fpp) {
		PROCPAR_ERROR("Cannot open procpar file for writing: " + path);
		return false;
	}
	
	ParameterList::const_iterator p;
	for (p = pp.begin(); p != pp.end(); p++)
		fpp << p->second << endl;
	fpp.close();
	return true;
}

const double RealValue(const ParameterList &pl, const string &name, const int index) {
	ParameterList::const_iterator p = pl.find(name);
	if (p == pl.end())
		PROCPAR_FAIL("Could not find parameter " + name);
	else
		return p->second.realValue(index);
}

const string &StringValue(const ParameterList &pl, const string &name, const int index) {
	ParameterList::const_iterator p = pl.find(name);
	if (p == pl.end())
		PROCPAR_FAIL("Could not find parameter " + name);
	else
		return p->second.stringValue(index);
}

}; // End namespace Recon