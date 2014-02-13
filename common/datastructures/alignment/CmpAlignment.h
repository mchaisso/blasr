#ifndef DATASTRUCTURES_ALIGNMENT_CMP_ALIGNMENT_H_
#define DATASTRUCTURES_ALIGNMENT_CMP_ALIGNMENT_H_
#include "datastructures/alignment/Alignment.h"
#include "Enumerations.h"
#include <vector>
#include <iostream>
#include <map>
#include <assert.h>
#include <algorithm>

using namespace std;
class CmpAlignmentBase  {
 public:
	//
	// For use in referencing alignment sets. TODO: subclass.
	//
	PlatformType platformType;
	int Z;
	unsigned int index, readGroupId, movieId, refSeqId; 
	unsigned int expId, runId, panel;
	unsigned int x, y;
	unsigned int rcRefStrand;
	unsigned int holeNumber;
	unsigned int offsetBegin, offsetEnd;
	unsigned int setNumber, strobeNumber, mapQV, nBackRead, nReadOverlap;
	unsigned int subreadId;
	unsigned int nMatch, nMismatch, nIns, nDel;
	vector<unsigned char> alignmentArray;
	vector<unsigned int> alignmentIndex;
	static map<string,int> columnNameToIndex;
  static bool initializedColumnNameToIndex;
  
  map<string, vector<UChar> > fields;
  
  
      
	unsigned int *GetAlignmentIndex() {
		return &alignmentIndex[0];
	}
	int GetAlignmentIndexSize() {
		return alignmentIndex.size();
	}
	unsigned int GetAlignedStrand() {
		return LookupColumnValue("AlignedStrand");
	}
  
	unsigned int GetRCRefStrand() {
		return LookupColumnValue("RCRefStrand");
	}

  // synonym 
  unsigned int GetTStrand() {
    return GetRCRefStrand();
  }

	bool GetX(int &xp) {
		if (alignmentIndex.size() > 0) {
			xp = alignmentIndex[columnNameToIndex["x"]];
			return true;
		}
		else {
			xp = -1;
			return false;
		}
	}

	unsigned int GetAlignmentId() {
		return LookupColumnValue("AlnID");
	}

	unsigned int GetX() {
		return LookupColumnValue("x");
	}
	
	unsigned int GetY() {
		return LookupColumnValue("y");
	}

	unsigned int GetMovieId() {
		return LookupColumnValue("MovieId");
	}

	unsigned int GetAlnGroupId() {
		return LookupColumnValue("AlnGroupId");
	}

	unsigned int GetReadGroupId() {
		return LookupColumnValue("ReadGroupId");
	}
	

	unsigned int LookupColumnValue(const char * columnName) {
		if (columnNameToIndex.find(columnName) != columnNameToIndex.end()) {
			int columnIndex = columnNameToIndex[columnName];
			return alignmentIndex[columnIndex];
		}
		else {
			cout << "ERROR, For now cmp files must contain a column " << columnName << endl;
			cout << "size of columnNameToIndex: " << columnNameToIndex.size() << endl;
			assert(0);
		}
	}		


	void InitializeColumnNameToIndex(vector<string> &columnNames) {
		int i;
		for (i = 0; i < columnNames.size(); i++ ){
			columnNameToIndex[columnNames[i]] = i;
		}
	}


	unsigned int GetHoleNumber() {
		return LookupColumnValue("HoleNumber");
	}

	unsigned int GetRefGroupId() {
		return LookupColumnValue("RefGroupId");
	}

	unsigned int GetRefSeqId() {
		return LookupColumnValue("RefSeqId");
	}

	unsigned int GetOffsetBegin() {
		return LookupColumnValue("offset_begin");
	}

	unsigned int GetOffsetEnd() {
		return LookupColumnValue("offset_end");
	}
	
	unsigned int GetQueryStart() {
		return LookupColumnValue("rStart");
	}

	unsigned int GetQueryEnd() {
		return LookupColumnValue("rEnd");
	}

	unsigned int GetRefStart() {
		return LookupColumnValue("tStart");
	}

	unsigned int GetRefEnd() {
		return LookupColumnValue("tEnd");
	}

  unsigned int GetNMatch() {
    return LookupColumnValue("nM");
  }

  unsigned int GetNMismatch() {
    return LookupColumnValue("nMM");
  }

  unsigned int GetNInsertions() {
    return LookupColumnValue("nIns");
  }

  unsigned int GetNDeletions() {
    return LookupColumnValue("nDel");
  }

  unsigned int GetMapQV() {
    return LookupColumnValue("MapQV");
  }

  unsigned int GetSubreadId() {
    return LookupColumnValue("SubreadId");
  }

  unsigned int GetStrobeNumber() {
    return LookupColumnValue("StrobeNumber");
  }

  unsigned int GetSetNumber() {
    return LookupColumnValue("SetNumber");
  }

  

  
  

	CmpAlignmentBase(PlatformType platformTypeP=Springfield) {
		platformType = platformTypeP;
	}

	void SetPlatformType(PlatformType platformTypeP) {
		platformType = platformTypeP;
	}
};

map<string,int> CmpAlignmentBase::columnNameToIndex;
bool initializedColumnNameToIndex = false;

class CmpAlignment : public CmpAlignmentBase {
 public:
	int qStrand, tStrand;
	int qStart, qLength;
	int tStart, tLength;
	//
	// Default constructor just calls the base constructor to initialize platoformType
  CmpAlignment(PlatformType ptype=Springfield) : CmpAlignmentBase(ptype) {
	}

	void StoreAlignmentIndex(unsigned int *alignmentIndexPtr, int alignmentIndexLength) {
		alignmentIndex.clear();
		alignmentIndex.insert(alignmentIndex.begin(), &alignmentIndexPtr[0], &alignmentIndexPtr[alignmentIndexLength]);
	}

	void StoreAlignmentArray(unsigned char* alignmentArrayPtr, int alignmentArrayLength) {
		alignmentArray.resize(alignmentArrayLength);
		unsigned int a;
		for (a = 0; a < alignmentArrayLength; a++ ){
			alignmentArray[a] = alignmentArrayPtr[a];
		}
	}

  template<typename T_Field>
    void StoreField(string fieldName, T_Field* fieldValues, int length) {
    fields[fieldName].resize(length);
    memcpy(&fields[fieldName][0], fieldValues, length * sizeof(T_Field));
  }    

  CmpAlignment &operator=(const CmpAlignment &rhs) {
    // deep copy the alignment index
    alignmentIndex.resize(rhs.alignmentIndex.size());
    copy(rhs.alignmentIndex.begin(), rhs.alignmentIndex.end(), alignmentIndex.begin());
    // deep copy the alignment array
    alignmentArray.resize(rhs.alignmentIndex.size());
    copy(rhs.alignmentArray.begin(), rhs.alignmentArray.end(), alignmentArray.begin());
    // copy fields
    Z = rhs.Z;
    index = rhs.index; readGroupId = rhs.readGroupId; movieId = rhs.movieId;
    refSeqId = rhs.refSeqId;
    expId = rhs.expId; runId = rhs.runId; panel = rhs.panel;
    x = rhs.x; y = rhs.y;
    rcRefStrand = rhs.rcRefStrand;
    holeNumber = rhs.holeNumber;
    offsetBegin = rhs.offsetBegin; offsetEnd = rhs.offsetEnd;
    setNumber = rhs.setNumber, strobeNumber = rhs.strobeNumber, mapQV = rhs.mapQV, nBackRead = rhs.nBackRead, nReadOverlap = rhs.nReadOverlap;
    subreadId = rhs.subreadId;
    nMatch = rhs.nMatch, nMismatch = rhs.nMismatch, nIns = rhs.nIns, nDel = rhs.nDel;
    return *this;
  }

  int operator<(const CmpAlignment &rhs) const {
    if (alignmentArray[1] == rhs.alignmentArray[1]) {
      if (alignmentArray[2] == rhs.alignmentArray[2]) {
        if (alignmentArray[10] == rhs.alignmentArray[10]) {
          return (alignmentArray[4] < rhs.alignmentArray[4]);
        }
        else {
          return alignmentArray[10] < rhs.alignmentArray[10];
        }
      }
      else {
        return alignmentArray[2] < rhs.alignmentArray[2];
      }
    }
    else {
      return alignmentArray[1] < rhs.alignmentArray[1];
    }
  }


};


class CmpFullAlignment : public Alignment {
 public:
	void StoreAlignmentArray(unsigned char *alignmentArrayPtr, int alignmentArrayLength) {
		//
		// This version of this template should eventually populate the
		// internal non-8 bit alignment structure so that it may be used elsewhere.
		// 
		cout << "NOT YET IMPLEMENTED." << endl;
		assert(0);
	}
};

#endif
