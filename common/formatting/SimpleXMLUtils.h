#ifndef SIMPLE_XML_UTILS_H_
#define SIMPLE_XML_UTILS_H_

#include <string>
#include <sstream>

using namespace std;

/*
 * WARNING !!! The functions here use TONS of temoprary variables and are 
 * very slow. Do not use for high-throughput xml generation.
 */

/*
 * none of this is really used... just the keyword value pair items.
 */

template<typename T_String, typename T_Value>
int OutputKeywordValuePair(T_String title, T_Value value, ostream out) {
	out << title << "=\"" << value << "\" ";
	return out.good();
}

template<typename T_String>
int OutputBeginElement(T_String title, ostream out) {
	out << "<" << title << " ";
	return 1;
}

template<typename T_String>
int OutputEndElement(T_String title, ostream out) {
	out << "/" << title;
	return 1;
}

template<typename T_String>
int OutputElement(T_String title, ostream out) {
	out << "<" << title << " ";
	return 1;
}

inline
void OutputEnd(ostream out) {
	out << "/>";
}


template<typename T_String, typename T_Value> 
T_String CreateKeywordValuePair(T_String title, T_Value value) {
	T_String keywordPair;
	stringstream sstrm;
	sstrm << title << "=\"" << value << "\"";
	keywordPair = sstrm.str();
	return keywordPair;
}

template<typename T_String>
T_String BeginDataEntry(T_String id, T_String data) {
	T_String entry;
	entry += "<";
	entry += id;
	entry += " ";
	entry += data;
	entry += ">";
	return entry;
}

template<typename T_String>
T_String EndDataEntry(T_String id) {
	T_String entry;
	entry = "<" + id + ">/";
	return entry;
}


template<typename T_String>
T_String CreateDataEntry(T_String id, T_String data) {
	T_String entry;
	entry += "<";
	entry += id;
	entry += " ";
	entry += data;
	entry += " />";
	return entry;
}


template<typename T_String> 
 T_String CreateStartEntry(T_String id, T_String data) {
	T_String entry;
	entry = "<" + id + " " + data + ">";
	return entry;
}

template<typename T_String>
T_String CreateEndEntry(T_String id) {
	T_String entry;
	entry = "<" + id + "/>";
	return entry;
}



	

#endif
