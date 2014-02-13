#ifndef SIMULATOR_CONTEXT_OUTPUT_LIST_H_
#define SIMULATOR_CONTEXT_OUTPUT_LIST_H_

#include "../../common/simulator/OutputList.h"
#include <stdlib.h>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include "utils.h"

using namespace std;

class ContextOutputList {
 public:
	map<string,OutputList*> outputMap;
	int contextLength;
	int ParsePair(string &output, string &outStr, int &outCount){ 
		string contextStr, outCountStr;
		int i;
		for (i = 0; i < output.size(); i++) {
			if (output[i] == '=') {
				int start = 0;
				while(start < i and (output[start] == ' ' or output[start] == '\t')) { start++;}
				outStr.assign(&output[start], i-start);
				outCountStr.assign(&output[i+1], output.size()-i-1);
				outCount = atoi(outCountStr.c_str());
				return 1;
			}
		}
		return 0;
	}
	
	int SampleRandomContext(string refContext, string &readContext) {
		//
		// Chec to see if there is a distribution for this ref context, if
		// not, just return the ref context so that things may proceed,
		// but return a fail code.
		//
		if (outputMap.find(refContext) == outputMap.end()) {
			readContext = refContext;
			return 0;
		}
		else {
			outputMap[refContext]->SelectRandomContect(readContext);
			return 1;
		}
	}


	void Read(string &inName) {
		ifstream in;
		CrucialOpen(inName, in, std::ios::in);
		Read(in);
	}

	void Read(ifstream &in) {
		int nLines = 0;
		contextLength = 0;
		while(in) {
			string context;
			if (!(in >> context)) break;
			
			contextLength = context.size();
			string outputsLine;
			getline(in,outputsLine);
			// parse ctx=num; pairs
			int i;
			int e;
			i = 0;
			if (outputMap.find(context) == outputMap.end()) {
				outputMap[context] = new OutputList;
			}
						
			while(i < outputsLine.size()) {
				e = i+1;
				while(e < outputsLine.size() and outputsLine[e] != ';') e++;
				string pairStr;
				pairStr.assign(&outputsLine[i], e-i);
				if (e > i + 1) {
					string outStr;
					int outCount;
					ParsePair(pairStr, outStr, outCount);
					outputMap[context]->AddOutput(outStr, outCount);
				}
				i = e + 1;
			}
			outputMap[context]->StoreCumulativeCounts();
		}
	}
  void Free() {
    map<string,OutputList*>::iterator mapIt;
    for (mapIt = outputMap.begin(); mapIt != outputMap.end(); ++mapIt) {
      delete mapIt->second;
    }
  }
};


#endif
