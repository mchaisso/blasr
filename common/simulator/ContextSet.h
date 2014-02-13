#ifndef SIMULATOR_CONTEXT_SET_H_
#define SIMULATOR_CONTEXT_SET_H_

#include "ContextSample.h"
#include <map>
#include <iostream>
using namespace std;
typedef map<string,ContextSample*> T_ContextSampleMap;

class ContextSampleMap : public T_ContextSampleMap {
 public:
 int contextLength;
 void Write(ofstream &out) {
	 T_ContextSampleMap::iterator mapIt, mapEnd;
	 int mapSize = this->size();
	 out.write((char*)&contextLength, sizeof(contextLength));
	 out.write((char*)&mapSize, sizeof(mapSize));
	 mapEnd = this->end();
	 for(mapIt = this->begin(); mapIt != mapEnd; ++mapIt) {
		 out.write(mapIt->first.c_str(), contextLength);
		 mapIt->second->Write(out);
	 }
 }
 void Read(ifstream &in) {
	 in.read((char*)&contextLength, sizeof(contextLength));
	 int numContext;
	 in.read((char*)&numContext, sizeof(numContext));
	 int i;
	 char *context = new char[contextLength+1];
	 context[contextLength] = '\0';
	 for (i = 0; i < numContext; i++) {
		 in.read(context, contextLength);
		 string contextString = context;
		 // Allocate the context
		 (*this)[contextString] = new ContextSample;
		 (*this)[contextString]->Read(in);
	 }
	 delete[] context;
 }
};
		 
#endif
