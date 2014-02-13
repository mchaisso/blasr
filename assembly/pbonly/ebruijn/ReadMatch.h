#ifndef EBRUIJN_READ_WORD_MATCH_H_
#define EBRUIJN_READ_WORD_MATCH_H_

#include <vector>
#include "DNASequence.h"

class ReadWordMatch {
 public:
  unsigned char *seq;
  vector<bool> pos;
  ReadWordMatch() {
    seq = NULL;
  }

  void Initialize(DNASequence &read) {
    pos.resize(read.length);
    pos.fill(pos.begin(), pos.end(), false);
  }

  ReadWordMatch(DNASequence &read) {
    Initialize(read);
  }

  void SetMatch(int p) {
    assert(p < pos.size());
    pos[p] = true;
  }
};

typedef vector<ReadWordMatch> ReadWordVector;

void InitializeFromReads(vector<DNASequence> &reads, ReadWordVector &readWordMatches) {
  int r;
  readWordMatches.resize(reads.size());
  for (r = 0; r < reads.size(); r++) {
    readWordMatches[r].Initialize(reads[r]);
  }
}





#endif
