#ifndef SIMULATOR_OUTPUT_LIST_H_
#define SIMULATOR_OUTPUT_LIST_H_
#include <string>
#include <vector>
using namespace std;


#include "CDFMap.h"

class Output {
 public:
	string output;
	int    count;
};

class OutputList {
 public:
	int totalCount;
	vector<Output> outputs;
	CDFMap<vector<Output>::iterator> cdf;

	OutputList() {
		totalCount = 0;
	}

	int AddOutput(string str, int count) {
		outputs.resize(outputs.size()+1);
		outputs[outputs.size()-1].output = str;
		outputs[outputs.size()-1].count  = count;
		totalCount += count;
	}
	
	void StoreCumulativeCounts() {

		int outputIndex;
		int total = 0;
		for (outputIndex = 0; outputIndex < outputs.size(); outputIndex++) {
			total += outputs[outputIndex].count;
			cdf.cdf.push_back(total);
			cdf.data.push_back(outputs.begin() + outputIndex);
		}
		assert(total == totalCount);
	}

	int ReturnUniformRandomContextIt(vector<Output>::iterator &outputIt) {
		return cdf.SelectRandomValue(outputIt);
	}
	int SelectRandomContect(string &outputContext) {
		vector<Output>::iterator outIt;
		ReturnUniformRandomContextIt(outIt);
		outputContext = outIt->output;
		return outIt - outputs.begin();
	}
		
};


#endif
