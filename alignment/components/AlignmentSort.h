#ifndef ALIGNMENT_COMPONENTS_ALIGNMENT_SORT_H_
#define ALIGNMENT_COMPONENTS_ALIGNMENT_SORT_H_

#include "datastructures/alignment/AlignmentCandidate.h"

class SortAlignmentPointersByScore {
public:
  int operator()(T_AlignmentCandidate *lhs, T_AlignmentCandidate* rhs) {
    if (lhs->score == rhs->score) {
      return lhs->tPos + lhs->tAlignedSeqPos < rhs->tPos + rhs->tAlignedSeqPos;
    }
    else {
      return lhs->score < rhs->score;
    }
  }
};

class SortAlignmentPointersByMapQV {
public:
  int operator()(T_AlignmentCandidate *lhs, T_AlignmentCandidate* rhs) {
    if (lhs->mapQV == rhs->mapQV) {
      if (lhs->score == rhs->score) {
        return lhs->tPos + lhs->tAlignedSeqPos < rhs->tPos + rhs->tAlignedSeqPos;
      }
      else {
        return lhs->score < rhs->score;
      }
    }
    else {
      return lhs->mapQV > rhs->mapQV;
    }
  }
};


#endif
