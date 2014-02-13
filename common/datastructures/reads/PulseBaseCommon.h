#ifndef DATASTRUCTURES_READS_PULSE_BASE_COMMON_H_
#define DATASTRUCTURES_READS_PULSE_BASE_COMMON_H_

// 
// This includes values that both pulse and base files must have.
//

#include "ScanData.h"
class PulseBaseCommon {
 public:
	ScanData scanData;
	vector<uint32_t> holeNumbers;

	float GetFrameRate() {
		return scanData.frameRate;
	}
	unsigned int GetNumFrames() {
		return scanData.numFrames;
	}
	string GetMovieName() {
		return scanData.movieName;
	}

    map<char, int> GetBaseMap() {
        return scanData.baseMap;
    }

	bool LookupReadIndexByHoleNumber(uint32_t holeNumber, int &readIndex) {
		vector<uint32_t>::iterator holeIt;
		if (holeNumbers.size() == 0) {
			return false;
		}
		holeIt = lower_bound(holeNumbers.begin(), holeNumbers.end(), holeNumber);
		if (holeIt == holeNumbers.end()) {
			return false;
		}
		if (*holeIt == holeNumber) {
			readIndex = holeIt - holeNumbers.begin();
			return true;
		}
		else {
			return false;
		}
	}

};

#endif
