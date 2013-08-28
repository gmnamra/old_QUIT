//
//  fdfFile.cpp
//  Agilent
//
//  Created by Tobias Wood on 23/08/2013.
//
//

#include "fdfFile.h"

namespace Agilent {

//*************************
#pragma mark fdfValue
//*************************
fdfType to_Type(const string &s) {
	if (s == "int") {
		return fdfType::Integer;
	} else if (s == "float") {
		return fdfType::Float;
	} else if (s == "char") {
		return fdfType::String;
	} else {
		throw(invalid_argument("Unknown field type: " + s));
	}
}

string to_string(const fdfType &t) {
	switch (t) {
		case fdfType::Integer: return "int"; break;
		case fdfType::Float: return "float"; break;
		case fdfType::String: return "char"; break;
	}
}



istream &operator>>(istream &is, fdfType &t) {
	string type_string;
	is >> type_string;
	t = to_Type(type_string);
	if (t == fdfType::String) {
		// Eat the * for strings
		string temp;
		getline(is, temp, '*');
	}
	return is;
}
		
ostream &operator<<(ostream &os, const fdfType &t) {
	os << to_string(t);
	if (t == fdfType::String) {
		cout << " *";
	}
	return os;
}

fdfValue::fdfValue(const int &ival) { m_value.i_val = ival; }
fdfValue::fdfValue(const float &fval) { m_value.f_val = fval; }
fdfValue::fdfValue(const string &sval) { m_value.s_val = new string(sval); }
fdfValue::~fdfValue() {
	if (m_type == fdfType::String) {
		delete m_value.s_val;
	}
}

template <>
const string fdfValue::value() const {
	switch (m_type) {
		case fdfType::Integer: return std::to_string(m_value.i_val);
		case fdfType::Float:   return std::to_string(m_value.f_val);
		case fdfType::String: return *(m_value.s_val);
	}
}

//*************************
#pragma mark fdfField
//*************************
const string &fdfField::name() const { return m_name; }
istream &operator>> (istream &is, fdfField &f) {
	string name, values;
	is >> f.m_type;
	getline(is, name, '=');
	if (name.find("[]") == string::npos) {
		string temp, value;
		getline(is, value, ';');
		getline(is, temp);
		switch (f.m_type) {
			case fdfType::Integer: f.m_values.emplace_back(stoi(value)); break;
			case fdfType::Float: f.m_values.emplace_back(stof(value)); break;
			case fdfType::String: {
				size_t leftQuote = value.find("\"");
				size_t rightQuote = value.rfind("\"");
				f.m_values.emplace_back(value.substr(leftQuote, rightQuote - leftQuote));
			} break;
		}
	} else {
		// We have an array
		name.erase(name.find("[]"));
		string temp;
		getline(is, temp, '{'); // Get rid of the opening brace
		getline(is, values, '}');
		getline(is, temp); // Eat the remainder of the line
		size_t thisComma = 0, nextComma = values.find(",");
		do {
			string value = values.substr(thisComma, nextComma - thisComma);
			switch (f.m_type) {
				case fdfType::Integer: f.m_values.emplace_back(stoi(value)); break;
				case fdfType::Float: f.m_values.emplace_back(stof(value)); break;
				case fdfType::String: {
					size_t leftQuote = value.find("\"");
					size_t rightQuote = value.rfind("\"");
					f.m_values.emplace_back(value.substr(leftQuote, rightQuote - leftQuote));
				} break;
			}
			thisComma = nextComma;
			nextComma = values.find(",", nextComma);
			if (nextComma == string::npos)
				nextComma = values.size();
		} while (thisComma < values.size());
	}
	return is;
}

ostream &operator<<(ostream &os, const fdfField &f) {
	os << f.m_type << f.m_name << " = ";
	if (f.m_values.size() > 1) cout << "{";
	for (auto &v : f.m_values) {
		cout << v.value<string>();
		if (f.m_values.size() > 1) cout << ",";
	}
	if (f.m_values.size() > 1) cout << "}";
	cout << ";";
	return os;
}

//******************
#pragma mark fdfFile
//******************

fdfFile::fdfFile(const string &path) {
	open(path);
}

void fdfFile::open(const string &path) {	
	fstream m_file(path, ios::in);
	string nextLine;
		
	if (!m_file) {
		throw(runtime_error("Could not open file: " + path));
	}
	m_path = path;
	
	if (!getline(m_file, nextLine) || nextLine != "#!/usr/local/fdf/startup") {
		throw(runtime_error("Could not magic string in file: " + path));
	}

	m_header.clear();
	while (m_file.peek() != '\0') {
		fdfField temp;
		m_file >> temp;
		m_header.insert(pair<string, fdfField>(temp.name(), temp));
	}
	getline(m_file, nextLine, '\0'); // Just to eat the null
	m_hdrSize = m_file.tellg();
	m_dtype = headerValue<string>("storage");
	m_dims[0] = headerValue<size_t>("matrix", 0);
	m_dims[1] = headerValue<size_t>("matrix", 1);
	if (headerValue<size_t>("rank") == 3) {
		m_dims[2] = headerValue<size_t>("matrix", 2);
	} else {
		m_dims[2] = 1;
	}
}

const size_t fdfFile::rank() const { return m_rank; }
const size_t fdfFile::dim(const size_t d) const {
	if (d > 2) {
		throw(invalid_argument("Tried to access dimension " + std::to_string(d) + " of file: " + m_path));
	} else {
		return m_dims[d];
	}
}
const size_t fdfFile::dataSize() const {
	return m_dims[0] * m_dims[1] * m_dims[2];
}

} // End namespace Agilent