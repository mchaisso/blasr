#ifndef WEIGHTED_INTERVAL_H_
#define WEIGHTED_INTERVAL_H_

#include <vector>
#include <queue>
#include "MatchPos.h"
#include "DNASequence.h"


using namespace std;
template<typename T_ChainedMatch> 
class WeightedInterval {
public:
	DNALength size; // not necessarily end - start + 1
	DNALength start;
	DNALength end;
	DNALength qStart, qEnd, qLength;
	mutable int readIndex;
	float pValue;
	vector<int> positions;
	mutable vector<T_ChainedMatch> matches;
	float pValueVariance, pValueNStdDev;
	float sizeVariance, sizeNStdDev;
	 int nAnchors;
	int totalAnchorSize;
	mutable bool isOverlapping;
	string tName;
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

	int operator<(const WeightedInterval<T_ChainedMatch> &intv) const {
		if (size == intv.size) {
			return start > intv.start;
		}
		else {
			return size < intv.size;
		}
	}
	int operator==(const WeightedInterval<T_ChainedMatch> &intv) const {
		return size == intv.size;
	}
	WeightedInterval<T_ChainedMatch>(){}
	void Init(int _size, int _start, int _end, int _readIndex, float _pValue) {
		size      = _size; 
		start     = _start; 
		end       = _end; 
		readIndex = _readIndex;
		pValue    = _pValue;
		qStart = 0;
		qEnd   = 0;
		qLength = 0;
    nAnchors = 0;
    totalAnchorSize = 0;
    pValueVariance =
      pValueNStdDev = 
      sizeVariance = 
      sizeNStdDev = 0;
	}

	WeightedInterval<T_ChainedMatch>(int _size, int _start, int _end, int _readIndex, float _pValue =0.0) {
		Init(_size, _start, _end, _readIndex, _pValue);
	}
	
	WeightedInterval<T_ChainedMatch>(int _size, int _start, int _end, int _readIndex, float _pValue, int _qStart, int _qEnd, int _qLength){
		Init(_size, _start, _end, _readIndex, _pValue);
		qStart    = _qStart;
		qEnd      = _qEnd;
		qLength   = _qLength;
	}

	WeightedInterval(int _size, unsigned int _nAnchors, unsigned int _totalAnchorSize, 
									 int _start, int _end, int _readIndex, float _pValue, 
									 int _qStart, int _qEnd, int _qLength, 
									 vector<T_ChainedMatch> &_matches) {
		Init(_size, _start, _end, _readIndex, _pValue);
		qStart    = _qStart;
		qEnd      = _qEnd;
		qLength   = _qLength;
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

template<typename T_ChainedMatch>
class CompareWeightedIntervalByPValue {
 public:
	int operator()(const WeightedInterval<T_ChainedMatch> &a, const WeightedInterval<T_ChainedMatch> &b) const {
    if (a.PValue() != b.PValue()) {
			return a.PValue() < b.PValue();
		}
		else {
			return a.start < b.start;
    }
	}
};

template<typename T_Chained_Anchor>
class WeightedIntervalSet :
  public multiset<WeightedInterval<T_Chained_Anchor>, CompareWeightedIntervalByPValue<T_Chained_Anchor> >  {
 public:
  int  maxSize;
	
	void RemoveContained(float maxRatio = 0.9) {
		typename WeightedIntervalSet<T_Chained_Anchor>::iterator it = (*this).begin();
		typename WeightedIntervalSet<T_Chained_Anchor>::iterator endIt = (*this).end();
		bool isContained = false;
		while (it != endIt and isContained == false) {
			DNALength curStart = (*it).qStart;
			DNALength curEnd   = (*it).qEnd;
			int curAnchorSize = (*it).totalAnchorSize;
			DNALength curReadLength = (*it).qLength;

			if ((*it).GetStrandIndex() == 1) {
				curStart = (*it).qLength - (*it).qEnd;
				curEnd   = (*it).qLength - (*it).qStart;
			}
			typename WeightedIntervalSet<T_Chained_Anchor>::iterator nextIt = it;
			++nextIt;
			while (nextIt != endIt) {
				float overlap = 0;
				float nextStart = (*nextIt).qStart;
				float nextEnd   = (*nextIt).qEnd;

				if ((*nextIt).GetStrandIndex() == 1) {
					nextEnd   = (*nextIt).qLength - (*nextIt).qStart;
					nextStart = (*nextIt).qLength - (*nextIt).qEnd;
				}	
				float overlapRatio = 0;
				if ((nextStart >= curStart and nextEnd < curEnd) or 
						(nextStart <= curStart and nextEnd >= curEnd)) {
					overlapRatio = 1;
				}
				else if (nextEnd < curStart or nextStart > curEnd) {
					overlapRatio = 0;
				}
				else if (nextStart < curStart and nextEnd > curStart) {
					overlapRatio = (nextEnd - curStart)*2/(nextEnd - nextStart + curEnd - curStart);
				}
				else if (nextStart < curEnd and nextEnd > curEnd) {
					overlapRatio = (curEnd - nextStart)*2/(nextEnd - nextStart + curEnd - curStart);
				}
				
				float anchorSizeRatio = ((float)(*nextIt).totalAnchorSize)/curAnchorSize;
				if (nextStart >= curStart and nextEnd <= curEnd and (overlapRatio > maxRatio or anchorSizeRatio < 0.25)) {
					typename WeightedIntervalSet<T_Chained_Anchor>::iterator skip = nextIt;
					++skip;
					this->erase(nextIt);
					nextIt = skip;
				}
				else {
					++nextIt;
				}
			}
			++it;
		}
	}

  WeightedIntervalSet<T_Chained_Anchor>() {
    maxSize = 0;
  }

  WeightedIntervalSet<T_Chained_Anchor>(int maxSizeP) : maxSize(maxSizeP) {} 
 
	bool insert( WeightedInterval<T_Chained_Anchor> &intv) {

    intv.isOverlapping = false;

		//
		// Make sure this interval is not contained inside any other
		// weighted intervals.  
		//

		typename WeightedIntervalSet<T_Chained_Anchor>::iterator it = (*this).begin();
		typename WeightedIntervalSet<T_Chained_Anchor>::iterator endIt = (*this).end();
		bool isContained = false;


		while (it != endIt and isContained == false) {
			
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
				typename WeightedIntervalSet<T_Chained_Anchor>::iterator next = it;
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
    if (this->size() >= maxSize and maxSize > 0) {
      typename WeightedIntervalSet<T_Chained_Anchor>::iterator last = (*this).end();
      last--;
 
      if (last->pValue < intv.pValue) {

        return false;
      }
    }

		if (isContained == false) {
      bool addInsert = false;
      if (this->size() == 0) {
        addInsert = true;
      }
      else {
        it = this->end();
        --it;
        if (this->size() < maxSize or (*it).pValue > intv.pValue) {
          addInsert = true;
          //
          // Keep the size of the stack the same if it is at the limit.
          //
          if (maxSize != 0 and this->size() >= maxSize and this->size() > 0) { 
            this->erase(it);
          }
        }
      }
      if (addInsert) {
        ((multiset<WeightedInterval<T_Chained_Anchor>, CompareWeightedIntervalByPValue<T_Chained_Anchor> >*)this)->insert(intv);
      }
      return true;
		}
    return false;
	}
};


#endif
