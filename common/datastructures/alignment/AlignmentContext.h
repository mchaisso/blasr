#ifndef DATASTRUCTURES_ALIGNMENT_ALIGNMENT_CONTEXT_H_
#define DATASTRUCTURES_ALIGNMENT_ALIGNMENT_CONTEXT_H_

#include "../../Enumerations.h"

class AlignmentContext {
public:
  bool isPrimary;
  int  subreadIndex;
  int  numProperlyAlignedSubreads;
  bool allSubreadsProperlyAligned;
  bool isFinal;
  int  nextSubreadPos;
  int  nextSubreadDir;
  bool hasNextSubreadPos;
  int  nSubreads;
  string rNext;
  string readGroupId;
  string chipId;
  AlignMode alignMode;
  int editDist;
  AlignmentContext() {
    isPrimary = true;
    subreadIndex = 0;
    isFinal   = true;
    nextSubreadPos = 0;
    hasNextSubreadPos = false;
    numProperlyAlignedSubreads = 0;
    allSubreadsProperlyAligned = false;
    nSubreads = 0;
    nextSubreadDir = 0;
    rNext = "";
    readGroupId = "";
    chipId = "";
    alignMode = NoAlignMode;
    editDist = 0;
  }
  bool IsFirst() {
    return subreadIndex == 0;
  }
  bool IsLast() {
    return subreadIndex == nSubreads-1;
  }
  bool AllSubreadsAligned() {
    if (numProperlyAlignedSubreads == nSubreads) {
      return true;
    }
    else {
      return false;
    }
  }
};

#endif
