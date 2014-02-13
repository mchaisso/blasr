#ifndef EBRUIJN_READ_WORD_MATCH_H_
#define EBRUIJN_READ_WORD_MATCH_H_

#include <vector>
#include <bitset>

#include "DNASequence.h"

#include <ostream>

class bitvect: public vector<bool> {
 public:
  int count() {
    int i;
    int n = 0;
    for (i = 0; i < this->size(); i++) {
      (*this)[i] && n++;
    }
    return n;
  }
};
class ReadWordMatch {
 public:
  unsigned char *seq;
  bitvect pos;

  void PrintPos(ostream &out) {
    int i;
    for (i = 0; i < pos.size(); i++) {
      if (pos[i]) { out << "1"; }
      else { out <<"0";}
    }
    out << endl;
  }


  vector<int> parents;
  void PrintParents(ostream &out) {
    int i;
    for (i = 0; i < parents.size(); i++) {
      out << parents[i] << " ";
    }
    out << endl;
  }

  int CreateParents() {
    parents.resize(pos.size());
    fill(parents.begin(), parents.end(), 0);
    return pos.count();
  }

  ReadWordMatch() {
    seq = NULL;
  }

  int size() {
    return pos.size();
  }

  void Initialize(DNASequence &read) {
    pos.resize(read.length);
    fill(pos.begin(), pos.end(), false);
  }

  ReadWordMatch(DNASequence &read) {
    Initialize(read);
  }

  void SetMatch(int p) {
    assert(p < pos.size());
    pos[p] = true;
  }

  int CountMatches() {
    return pos.count();
  }
  
};

typedef vector<ReadWordMatch> ReadWordMatchVector;

template<typename T_Read>
void InitializeFromReads(vector<T_Read> &reads, ReadWordMatchVector &readWordMatches) {
  int r;
  readWordMatches.resize(reads.size());
  for (r = 0; r < reads.size(); r++) {
    readWordMatches[r].Initialize(reads[r]);
  }
}





#endif
