#ifndef ADAPTERS_SEQUENCE_READER_H_
#define ADAPTERS_SEQUENCE_READER_H_

#include "SMRTSequence.h"

class SequenceReader {
 public:
  virtual void Initialize(string &inFileName);
  virtual bool GetNext(SMRTSequence &seq);
  virtual bool Advance(int numSeq);
};
  

class DefaultAdvancePolicy {
 public:
  int GetNumberToAdvance() {
    return 1;
  }
}

template<typename T_Reader, typename T_Sequence=SMRTSequence, typename T_AdvancePolicy=DefaultAdvancePolicy>
class SequenceReaderBase {
 public:
};

#endif
