#ifndef ALIGNMENT_COMPONENTS_QV_PURITY_H_
#define ALIGNMENT_COMPONENTS_QV_PURITY_H_

#include "FASTQSequence.h"

int CountZero(unsigned char *ptr, int length) {
  int i;
  int nZero = 0;
  for (i = 0; i < length; i++) {
    if (ptr[i] == 0) { ++nZero; }
  }
  return nZero;
}


bool ReadHasMeaningfulQualityValues(FASTQSequence &sequence) {
  if (sequence.qual.Empty() == true) {
    return 0;
  }
  else {
    int q;
    int numZero=0, numNonZero=0;
    if (sequence.qual.data == NULL) {
      return false;
    }
    numZero = CountZero(sequence.qual.data, sequence.length);
    numNonZero = sequence.length - numZero;
    int subNumZero = 0, subNonZero = 0;

    if (sequence.substitutionQV.data == NULL) {
      return false;
    }
    subNumZero = CountZero(sequence.substitutionQV.data, sequence.length);
    subNonZero = sequence.length - subNumZero;

    if (numZero < 0.5*numNonZero and subNumZero < 0.5 * subNonZero) {
       return true;
     }
    else {
      return false;
    }
  }
}


#endif
