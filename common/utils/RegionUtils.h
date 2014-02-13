#ifndef UTILS_REGION_UTILS_H_
#define UTILS_REGION_UTILS_H_
#include <algorithm>
#include "SMRTSequence.h"
#include "datastructures/reads/ReadInterval.h"
#include "datastructures/reads/RegionTable.h"

bool LookupHQRegion(int holeNumber, RegionTable &regionTable, int &start, int &end, int &score) {
	int regionLowIndex, regionHighIndex;
	regionLowIndex = regionHighIndex = 0;
	regionTable.LookupRegionsByHoleNumber(holeNumber, regionLowIndex, regionHighIndex);
	bool readHasGoodRegion = true;
	int  regionIndex = regionLowIndex;
	while (regionIndex < regionHighIndex and 
				 regionTable.GetType(regionIndex) != HQRegion) {
		regionIndex++;
	}
	
	if (regionIndex == regionHighIndex) {
    start = end = score = 0;
		return false;
	}
	else {
		start = regionTable.GetStart(regionIndex);
		end   = regionTable.GetEnd(regionIndex);
    score = regionTable.GetScore(regionIndex);
		return true;
	}
}


template<typename T_Sequence>
bool MaskRead(T_Sequence &fastaRead, ZMWGroupEntry &zmwData, RegionTable &regionTable) {
	int regionIndex;						 
	int regionLowIndex, regionHighIndex;
	regionLowIndex = regionHighIndex = 0;
	regionTable.LookupRegionsByHoleNumber(zmwData.holeNumber, regionLowIndex, regionHighIndex);
	bool readHasGoodRegion = true;

	DNALength readPos;

	regionIndex = regionLowIndex;
	int lastHQRegionIndex;
	
	int hqRegionStart=0, hqRegionEnd=0, hqRegionScore = 0;
	readHasGoodRegion = LookupHQRegion(zmwData.holeNumber, regionTable, hqRegionStart, hqRegionEnd, hqRegionScore);

	//
	// Mask off the low quality portion of this read.
	//
	for (readPos = 0; (readPos < hqRegionStart and
											 readPos < fastaRead.length); readPos++) {
		fastaRead.seq[readPos] = 'N';
	}
	for (readPos = hqRegionEnd; readPos < fastaRead.length; readPos++) {
		fastaRead.seq[readPos] = 'N';
	}

	//
	// Look to see if there is region information provided, but the entire read is bad.
	//
	if (hqRegionEnd == hqRegionStart) {
		//
		// This read is entirely bad, flag that.
		//
		readHasGoodRegion = false;
	}

	return readHasGoodRegion;
}


template<typename T_Sequence>
bool GetReadTrimCoordinates(T_Sequence &fastaRead,
														ZMWGroupEntry &zmwData,
														RegionTable &regionTable,
														DNALength &readStart,
														DNALength &readEnd,
                            int &score) {

	int regionIndex;						 
	int regionLowIndex, regionHighIndex;
	regionLowIndex = regionHighIndex = 0;
	regionTable.LookupRegionsByHoleNumber(zmwData.holeNumber, regionLowIndex, regionHighIndex);
	bool readHasGoodRegion = true;

	DNALength readPos;

	regionIndex = regionLowIndex;
	int lastHQRegionIndex;
	
	while (regionIndex < regionHighIndex and 
				 regionTable.GetType(regionIndex) != HQRegion) {
		regionIndex++;
	}
	
	if (regionIndex < regionHighIndex ) {
		readStart = regionTable.GetStart(regionIndex);
		readEnd   = regionTable.GetEnd(regionIndex);
    score     = regionTable.GetScore(regionIndex);
		return true;
	}
	else {
		readStart = 0;
		readEnd   = fastaRead.length;
		return false;
	}
}

template<typename T_Sequence>
bool TrimRead(T_Sequence &fastaRead, ZMWGroupEntry &zmwData, RegionTable &regionTable, T_Sequence &trimmedRead) {

	DNALength readStart, readEnd;
	GetReadTrimCoordinates(fastaRead, zmwData, regionTable, readStart, readEnd);
	if (readEnd - readStart > 0) {
		trimmedRead.CopySubsequence((FASTQSequence&)fastaRead, 
																readStart, readEnd);
		// signal that the read has a good region.
		return true;
	}
	else {

		//
		// There is no information for this read. Make it skipped.
		//
		trimmedRead.seq = NULL;
		trimmedRead.CopyTitle(fastaRead.title);
		// signal this read has no good region.
		return false;
	}
}

class CompareRegionIndicesByStart {
public:
	RegionTable *regionTablePtr;
	int operator()(const int a, const int b) const {
		if (regionTablePtr->GetStart(a) == regionTablePtr->GetStart(b)) {
			return (regionTablePtr->GetEnd(a) < regionTablePtr->GetEnd(b));
		}
		else {
			return (regionTablePtr->GetStart(a) < regionTablePtr->GetStart(b));
		}
	}
};

		
int SortRegionIndicesByStart(RegionTable &regionTable, vector<int> &indices) {
	CompareRegionIndicesByStart cmpFctr;
	cmpFctr.regionTablePtr = &regionTable;
	std::sort(indices.begin(), indices.end(), cmpFctr);
  return indices.size();
}

class OrderRegionsByReadStart {
 public:
	int operator()(const ReadInterval &lhs, const ReadInterval &rhs) const {
		return lhs.start < rhs.start;
	}
};

int FindRegionIndices(unsigned int holeNumber, RegionTable *regionTablePtr, int &regionLowIndex, int &regionHighIndex) {
	int regionIndex;						 
  regionLowIndex = regionHighIndex = 0;
	regionTablePtr->LookupRegionsByHoleNumber(holeNumber, regionLowIndex, regionHighIndex);  
  return regionHighIndex - regionLowIndex;
}

int FindRegionIndices(SMRTSequence &read, RegionTable *regionTablePtr, int &regionLowIndex, int &regionHighIndex) {
  return FindRegionIndices(read.zmwData.holeNumber, regionTablePtr, regionLowIndex, regionHighIndex);
}

//
// Collect region indices for either all region types, or just a few specific region types.
//

int CollectRegionIndices(SMRTSequence &read, RegionTable &regionTable, vector<int> &regionIndices,
                         RegionType *regionTypes=NULL, int numRegionTypes = 0) {
  int regionLow, regionHigh;
  int prevNumRegionIndices = regionIndices.size();
  if (FindRegionIndices(read, &regionTable, regionLow, regionHigh)) {
    int i;
    for (i = regionLow; i < regionHigh; i++) {
      if (regionTypes == NULL) {
        regionIndices.push_back(i);
      }
      else {
        int t;
        for (t = 0; t < numRegionTypes; t++) {
          if (regionTable.GetType(i) == regionTypes[t]) {
            regionIndices.push_back(i);
            break;
          }
        }
      }
    }
  }
  return regionIndices.size() - prevNumRegionIndices;
}




template<typename T_Sequence>
void CollectSubreadIntervals(T_Sequence &read, RegionTable *regionTablePtr, vector<ReadInterval> &subreadIntervals, bool byAdapter=false) {
	int regionIndex;						 
	int regionLowIndex, regionHighIndex;
	regionLowIndex = regionHighIndex = 0;
	regionTablePtr->LookupRegionsByHoleNumber(read.zmwData.holeNumber, regionLowIndex, regionHighIndex);
	if (byAdapter == false) {
		for (regionIndex = regionLowIndex; regionIndex < regionHighIndex; regionIndex++) {
			if (regionTablePtr->GetType(regionIndex) ==  Insert) {
				subreadIntervals.push_back(ReadInterval(regionTablePtr->table[regionIndex].row[RegionAnnotation::RegionStart],
																								regionTablePtr->table[regionIndex].row[RegionAnnotation::RegionEnd],
                                                regionTablePtr->table[regionIndex].row[RegionAnnotation::RegionScore]));
			}
		}
	}
	else {
		vector<int> adapterIntervalIndices;
		for (regionIndex = regionLowIndex; regionIndex < regionHighIndex; regionIndex++) {
			if (regionTablePtr->GetType(regionIndex) == Adapter) {
				adapterIntervalIndices.push_back(regionIndex);
			}
		}
		// Sort indices so that the intervals appear in order on the read.
		SortRegionIndicesByStart(*regionTablePtr, adapterIntervalIndices);
		int curIntervalStart = 0;
		int i;
		if (adapterIntervalIndices.size() == 0) {
			subreadIntervals.push_back(ReadInterval(0, read.length));
		}
		else {
			subreadIntervals.push_back(ReadInterval(0, regionTablePtr->table[adapterIntervalIndices[0]].row[RegionAnnotation::RegionStart]));
			for (i = 0; i + 1 < adapterIntervalIndices.size() ; i++) {
				subreadIntervals.push_back(ReadInterval(regionTablePtr->table[adapterIntervalIndices[i  ]].row[RegionAnnotation::RegionEnd],
																								regionTablePtr->table[adapterIntervalIndices[i+1]].row[RegionAnnotation::RegionStart]));
			}
			subreadIntervals.push_back(ReadInterval(regionTablePtr->table[adapterIntervalIndices[adapterIntervalIndices.size()-1]].row[RegionAnnotation::RegionEnd],
																							read.length));
		}
	}
	sort(subreadIntervals.begin(), subreadIntervals.end(), OrderRegionsByReadStart());
}


// Get all adapter intervals of a ZMW.
// Input:
//   read - read.zmwData.holeNumber specifies the zmw.
//   regionTablePtr - a pointer to a region table.
// Output:
//   adapterIntervals - where to assign all adapter intervals of the zmw
template<typename T_Sequence>
void CollectAdapterIntervals(T_Sequence &read, RegionTable *regionTablePtr, vector<ReadInterval> &adapterIntervals) {
  assert(regionTablePtr != NULL);
  int regionLowIndex = 0, regionHighIndex = 0;
  regionTablePtr->LookupRegionsByHoleNumber(read.zmwData.holeNumber, regionLowIndex, regionHighIndex);
  for (int regionIndex = regionLowIndex; regionIndex < regionHighIndex; regionIndex++) {
    if (regionTablePtr->GetType(regionIndex) ==  Adapter) {
      adapterIntervals.push_back(ReadInterval(
          regionTablePtr->table[regionIndex].row[RegionAnnotation::RegionStart],
          regionTablePtr->table[regionIndex].row[RegionAnnotation::RegionEnd],
          regionTablePtr->table[regionIndex].row[RegionAnnotation::RegionScore]));
    }
  }
}


// Given a vecotr of ReadInterval objects and their corresponding 
// directions, intersect each object with an interval 
// [hqStart, hqEnd), if there is no intersection or the intersected
// interval is less than minIntervalLength, remove this object and
// their corresponding directions; otherwise, replace this object 
// with the intersected interval and keep their directions. 
// Return index of the (left-most) longest subread interval in the
// updated vector.
int GetHighQualitySubreadsIntervals(vector<ReadInterval> & subreadIntervals, 
                                    vector<int> & subreadDirections,
                                    int hqStart, int hqEnd, 
                                    int minIntervalLength = 0) {
    // Avoid using vector.erase() when possible, as it is slow.
    int ret = -1;
    int maxLength = 0;
    assert(subreadIntervals.size() == subreadDirections.size());
    vector<ReadInterval> subreadIntervals2; 
    vector<int> subreadDirections2;
    for(int i = 0; i < int(subreadIntervals.size()); i++) {
        int & thisStart = subreadIntervals[i].start;
        int & thisEnd   = subreadIntervals[i].end;
        if (thisStart >= hqEnd or thisEnd <= hqStart) {
            continue;
        } 
        if (thisStart < hqStart and thisEnd > hqStart) {
            thisStart = hqStart;
        }
        if (thisStart < hqEnd   and thisEnd > hqEnd  ) {
            thisEnd   = hqEnd;
        }
        if (thisEnd - thisStart >= minIntervalLength) {
            if (maxLength < thisEnd - thisStart) {
                ret = subreadIntervals2.size();
                maxLength = thisEnd - thisStart;
            }
            subreadIntervals2.push_back(subreadIntervals[i]);
            subreadDirections2.push_back(subreadDirections[i]);
        }
    }
    subreadIntervals  = subreadIntervals2;
    subreadDirections = subreadDirections2;
    return ret;
}

// Given a vector of subreads and a vector of adapters, return
// index of the (left-most) longest subread which has both
// adapters before & after itself.
int GetLongestFullSubreadIndex(vector<ReadInterval> & subreadIntervals,
                               vector<ReadInterval> & adapterIntervals) {
  int longestLength = 0;
  int index = -1; // Index of the longest fullpass subread.
  for(int i = 0; i < subreadIntervals.size(); i++) {
    ReadInterval & subread = subreadIntervals[i];
    bool ladapter = false, radapter = false;
    for(int j = 0; j < adapterIntervals.size(); j++) {
      ReadInterval & adapter = adapterIntervals[j];
      if (abs(subread.start - adapter.end) < 10) {
          ladapter = true;
      } else if(abs(subread.end - adapter.start) < 10) {
          radapter = true;
      }
      if (ladapter && radapter) {
        if (longestLength < subread.end - subread.start) {
          longestLength = subread.end - subread.start;
          index = i;
        } else {
          break;
        }
      }
    }
  }
  return index;
}


// Create a vector of n directions consisting of interleaved 0 and 1s.
void CreateDirections(vector<int> & directions, const int & n) {
    directions.clear();
    directions.resize(n);
    for(int i = 0; i < n; i++) {
        directions[i] = i % 2;
    }
}

// Flop all directions in the given vector, if flop is true.
void UpdateDirections(vector<int> & directions, bool flop = false) {
  if (not flop) return;
  for (int i = 0; i < int(directions.size()); i++) {
    assert(directions[i] == 0 or directions[i] == 1);
    directions[i] = (directions[i] == 0)?1:0;
  }
}


#endif
