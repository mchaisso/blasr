#ifndef SIMULATOR_LENGTH_HISTOGRAM_H_
#define SIMULATOR_LENGTH_HISTOGRAM_H_

#include "simulator/CDFMap.h"
#include "datastructures/alignment/CmpAlignment.h"
#include "utils.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class LengthHistogram {
 public:
	bool interpolate;
	CDFMap<int> lengthHistogram;

	LengthHistogram() {
		interpolate = false;
	}
	
	int Read(string &inName) {
		ifstream in;
		CrucialOpen(inName, in, std::ios::in);
		return Read(in);
	}

	int Read(ifstream &in) {
		while(in) {
			int length, count;
			in >> length;
			in >> count;
			lengthHistogram.data.push_back(length);
			if (lengthHistogram.cdf.size() == 0) {
				lengthHistogram.cdf.push_back(count);
			}
			else {
				lengthHistogram.cdf.push_back(lengthHistogram.cdf[lengthHistogram.cdf.size()-1] + count);
			}
		}
	}

	void GetRandomLength(int &length) {
		int index = lengthHistogram.SelectRandomValue(length);
		if (interpolate == true) {
			if (index < lengthHistogram.data.size() - 1) {
				int nextLength = lengthHistogram.data[index+1];
				length += RandomInt(nextLength - length);
			}
		}
	}

  void BuildFromAlignmentLengths(vector<int> &lengths) {
    int i;
    sort(lengths.begin(), lengths.end());
    int f;
		//
		// For now, just use all lengths.
		//
    for (f = 0, i = 1; i < lengths.size(); i++, f++) {
			lengthHistogram.data.push_back(lengths[f]);
			lengthHistogram.cdf.push_back(i);
    }
    /* Tests:
     * indices                0 1 2 3 4  5  6
     * lengths:               1 3 5 9 10 10 11
     * lengthHistogram.data:  1 3 5 9 10 11
     * lengthHistogram.cdf :  1 2 3 4 6  7
     *
     * lengths:               1 3 5 9 10 11
     * lengthHistogram.data:  1 3 5 9 10 11
     * lengthHistogram.cdf :  1 2 3 4 5  6
     *
     * lengths:               10
     * lengthHistogram.data:  10
     * lengthHistogram.cdf :  1 
     */ 
  }
};

#endif
