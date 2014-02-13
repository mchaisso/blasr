#ifndef DATASTRUCTURES_ALIGNMENT_GAP_LIST_H_
#define DATASTRUCTURES_ALIGNMENT_GAP_LIST_H_

#include <vector>
using namespace std;



class Gap {
 public:
  enum GapSeq {Query, Target};
  GapSeq seq;
  int length;
  Gap() {
    seq    = Query;
    length = 0;
  }
  Gap(GapSeq seqP, int lengthP) {
    seq = seqP;
    length = lengthP;
  }
};

typedef vector<Gap> GapList;

#endif
