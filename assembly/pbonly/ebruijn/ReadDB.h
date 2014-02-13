#ifndef EBRUIJN_READ_DB_H_
#define EBRUIJN_READ_DB_H_
#include "FASTASequence.h"
#include <map>
#include <string>
void BuildReadNameToIndexMap(vector<FASTASequence > &reads,
                             map<string, int> &readNameToIndex) {
  int i;
  for (i = 0; i < reads.size(); i++) {
    readNameToIndex[reads[i].title] = i;
  }
}

#endif
