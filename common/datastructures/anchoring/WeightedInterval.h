#ifndef WEIGHTED_INTERVAL_H_
#define WEIGHTED_INTERVAL_H_

#include <vector>
#include <queue>
#include "MatchPos.h"
#include "../../DNASequence.h"


using namespace std;

class WeightedInterval {
public:
	DNALength size; // not necessarily end - start + 1
	DNALength start;
	DNALength end;
	DNALength qStart, qEnd;
	int readIndex;
	float pValue;
	vector<int> positions;
	vector<ChainedMatchPos> matches;
  float pValueVariance, pValueNStdDev;
  float sizeVariance, sizeNStdDev;
  int nAnchors;
  int totalAnchorSize;
  bool isOverlapping;
  int GetStrandIndex() const {
    return readIndex;
  }

  void SetPValueVariance(float v) {
    pValueVariance = v;
  }

  void SetPValueNStdDev(float v) {
    pValueNStdDev = v;
  }

  void SetSizeVariance(float v) {
    sizeVariance = v;
  }

  void SetSizeNStdDev(float v) {
    sizeNStdDev = v;
  }

	int operator<(const WeightedInterval &intv) const {
		if (size == intv.size) {
			return start > intv.start;
		}
		else {
			return size < intv.size;
		}
	}
	int operator==(const WeightedInterval &intv) const {
		return size == intv.size;
	}
	WeightedInterval(){}
	void Init(int _size, int _start, int _end, int _readIndex, float _pValue) {
		size      = _size; 
		start     = _start; 
		end       = _end; 
		readIndex = _readIndex;
		pValue    = _pValue;
		qStart = 0;
		qEnd   = 0;
    nAnchors = 0;
    totalAnchorSize = 0;
    pValueVariance =
      pValueNStdDev = 
      sizeVariance = 
      sizeNStdDev = 0;
	}

	WeightedInterval(int _size, int _start, int _end, int _readIndex, float _pValue =0.0) {
		Init(_size, _start, _end, _readIndex, _pValue);
	}
	
	WeightedInterval(int _size, int _start, int _end, int _readIndex, float _pValue, int _qStart, int _qEnd){
		Init(_size, _start, _end, _readIndex, _pValue);
		qStart    = _qStart;
		qEnd      = _qEnd;
	}

	WeightedInterval(int _size, unsigned int _nAnchors, unsigned int _totalAnchorSize, int _start, int _end, int _readIndex, float _pValue, int _qStart, int _qEnd, vector<ChainedMatchPos> &_matches) {
		Init(_size, _start, _end, _readIndex, _pValue);
		qStart    = _qStart;
		qEnd      = _qEnd;
		matches   = _matches;
    nAnchors  = _nAnchors;
    totalAnchorSize = _totalAnchorSize;
	}

	float PValue() const {
		return pValue;
	}
	int Size() const {
		return size;
	}
};

class CompareWeightedIntervalByPValue {
 public:
	int operator()(const WeightedInterval &a, const WeightedInterval &b) const {
    if (a.PValue() != b.PValue()) {
			return a.PValue() < b.PValue();
		}
		else {
			return a.start < b.start;
    }
	}
};

typedef vector<WeightedInterval> WeightedIntervalVector;
typedef multiset<WeightedInterval, CompareWeightedIntervalByPValue> T_WeightedIntervalMultiSet;

class WeightedIntervalSet : public T_WeightedIntervalMultiSet {
 public:
  int  maxSize;

  WeightedIntervalSet() {
    maxSize = 0;
  }

  WeightedIntervalSet(int maxSizeP) : maxSize(maxSizeP) {} 
 
	bool insert(WeightedInterval &intv) {

    intv.isOverlapping = false;

		//
		// Make sure this interval is not contained inside any other
		// weighted intervals.  
		//

		WeightedIntervalSet::iterator it = (*this).begin();
		WeightedIntervalSet::iterator endit = (*this).end();
		bool isContained = false;
		while (it != endit and isContained == false) {

			if (intv.qStart >= (*it).qStart and intv.qEnd <= (*it).qEnd and
					intv.start  >= (*it).start and intv.end <= (*it).end and 
          intv.readIndex == (*it).readIndex and
          intv.pValue >= (*it).pValue) {
				//
				// This already overlaps an existing interval, don't bother
				// trying to add it.
				//
				isContained = true;
        intv.isOverlapping = true;
			}
			else if ((*it).start >= intv.start and (*it).end <= intv.end and
							 (*it).qStart >= intv.qStart and (*it).qEnd <= intv.qEnd and
               (*it).readIndex == intv.readIndex and
               (*it).pValue >= intv.pValue) {
				WeightedIntervalSet::iterator next = it;
				++next;
				this->erase(it);
				it = next;
			}
			else {
				++it;
			}
		}

    //
    // Take a peek to see if this interval is too low of a score to
    // bother attempting to add at all. 
    //
    if (size() >= maxSize and maxSize > 0) {
      WeightedIntervalSet::iterator last = (*this).end();
      last--;
 
      if (last->pValue < intv.pValue) {

        return false;
      }
    }

		if (isContained == false) {
      bool addInsert = false;
      if (size() == 0) {
        addInsert = true;
      }
      else {
        it = end();
        --it;
        if (size() < maxSize or (*it).pValue > intv.pValue) {
          addInsert = true;
          //
          // Keep the size of the stack the same if it is at the limit.
          //
          if (maxSize != 0 and size() >= maxSize and size() > 0) { 
            erase(it);
          }
        }
      }
      if (addInsert) {
        ((T_WeightedIntervalMultiSet*)this)->insert(intv);
      }
      return true;
		}
    return false;
	}
};


#endif
