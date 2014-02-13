#ifndef DATASTRUCTURES_READS_BAS_FILE_H_
#define DATASTRUCTURES_READS_BAS_FILE_H_

#include <vector>
#include <algorithm>
#include "HoleXY.h"
#include "PulseBaseCommon.h"
#include "ScanData.h"
#include "Enumerations.h"
#include "SMRTSequence.h"

class CompareHoleXY {
 public:
	bool operator()(const HoleXY &lhs, const HoleXY &rhs) const {
		if (lhs.xy[0] == rhs.xy[0]) {
			return lhs.xy[1] < rhs.xy[1];
		}
		else {
			return lhs.xy[0] < rhs.xy[0];
		}
	}
};

class BaseFile : public PulseBaseCommon {
 public:
	vector<unsigned char> baseCalls;
	vector<uint8_t> holeStatus;
	vector<HoleXY>   holeXY;
	vector<uint16_t> basWidthInFrames;
	vector<uint16_t> preBaseFrames;
	vector<int> pulseIndex;
	vector<unsigned char> qualityValues;
	vector<unsigned char> deletionQV;
	vector<unsigned char> deletionTag;
	vector<unsigned char> insertionQV;
	vector<unsigned char> substitutionQV;
	vector<unsigned char> substitutionTag;
	vector<unsigned char> mergeQV;
	vector<int>  readLengths;
  vector<int>  readStartPositions;
	int nReads;
	int nBases;

	bool LookupReadIndexByXY(uint16_t x, uint16_t y, int &index) {
		int16_t xy[2];
		xy[0] = x; xy[1] = y;
		vector<HoleXY>::iterator holeIt;
		holeIt = lower_bound(holeXY.begin(), holeXY.end(), xy);
		if ((*holeIt).xy[0] == xy[0] and (*holeIt).xy[1] == xy[1]) {
			index = holeIt - holeXY.begin();
			return true;
		}
		else {
			return false;
		}
	}

  
  template<typename T>
    void CopyArray(vector<T> &fullArray, int pos, int length, T*dest) {
    memcpy(dest, &fullArray[pos], sizeof(T) * length);
  }

  void CopyReadAt(uint32_t readIndex, SMRTSequence &read) {
    assert(holeNumbers.size() > readIndex);
    read.zmwData.holeNumber = holeNumbers[readIndex];
    if (holeXY.size() > 0) {
      assert(holeXY.size() > readIndex);
      read.zmwData.x = holeXY[readIndex].xy[0];
      read.zmwData.y = holeXY[readIndex].xy[1];
    }
    
    int startPos = readStartPositions[readIndex];
    int readLength = readLengths[readIndex];
    read.length = readLength;
    read.Allocate(readLength);
    if (baseCalls.size() > 0) {
      assert(baseCalls.size() >= readLength + startPos);
      CopyArray(baseCalls, startPos, readLength, read.seq);
    }
    if (qualityValues.size() > 0) {
      assert(qualityValues.size() >= readLength + startPos);
      CopyArray(qualityValues, startPos, readLength, read.qual.data);
    }
    if (basWidthInFrames.size() > 0) {
      assert(basWidthInFrames.size() >= readLength + startPos);
      CopyArray(basWidthInFrames, startPos, readLength, read.widthInFrames);
    }
    if (deletionQV.size() > 0) {
      assert(deletionQV.size() >= readLength + startPos);
      CopyArray(deletionQV, startPos, readLength, read.deletionQV.data);
    }
    if (deletionTag.size() > 0) {
      assert(deletionTag.size() >= readLength + startPos);
      CopyArray(deletionTag, startPos, readLength, read.deletionTag);
    }
    if (insertionQV.size() > 0) {
      assert(insertionQV.size()  >= readLength + startPos);
      CopyArray(insertionQV, startPos, readLength, read.insertionQV.data);
    }
    if (substitutionQV.size() > 0) {
      assert(substitutionQV.size() >= readLength + startPos);
      CopyArray(substitutionQV, startPos, readLength, read.substitutionQV.data);
    }
    if (mergeQV.size() > 0) {
      assert(mergeQV.size() >= readLength + startPos);
      CopyArray(mergeQV, startPos, readLength, read.mergeQV.data);
    }
    if (substitutionTag.size() > 0) {
      assert(substitutionTag.size() >= readLength + startPos);
      CopyArray(substitutionTag, startPos, readLength, read.substitutionTag);
    }
    if (preBaseFrames.size() > 0) {
      assert(preBaseFrames.size() >= readLength + startPos);
      CopyArray(preBaseFrames, startPos, readLength, read.preBaseFrames);
    }
    
  }

};

#endif
