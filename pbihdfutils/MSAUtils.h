#ifndef PBIHDFUTILS_MSA_UTILS_H_
#define PBIHDFUTILS_MSA_UTILS_H_

#include "InsertedString.h"

void ComputePositions(ByteAlignment &aln, char charMap[], vector<int> &pos) {
  pos.resize(aln.size());
  fill(pos.begin(), pos.end(), -1);
  int i;
  int curPos = 0;
  for (i = 0; i < aln.size(); i++) {
    if (charMap[aln[i]] != ' ') {
      pos[i] = curPos;
      ++curPos;
    }
  }
}

void ComputeQueryPositions(ByteAlignment &aln, vector<int> &queryPos) {
  ComputePositions(aln, QueryChar, queryPos);
}

void ComputeRefPositions(ByteAlignment &aln, vector<int> &refPos) {
  ComputePositions(aln, RefChar, refPos);
}


void StoreInsertedStrings(ByteAlignment &aln, vector<int> &refPositions, vector<int> &queryPositions,
                          InsertedStringList &insertions, int start, int end) {
  int p = start;
  //
  // Do not store insertions at the beginning of the string since they
  // start to the left of where printing begins.
  // 
  while (p < end and RefChar[aln[p]] == ' ') { p++; }
  while (p < end) {
    // look to see if an insertion starts at p
    if (RefChar[aln[p]] == ' ') {
      int insEnd = p;
      string insertedString = "";
      while (insEnd < end and RefChar[aln[insEnd]] == ' ') {
        insertedString.push_back(QueryChar[aln[insEnd]]);
        insEnd++;
      }
      //
      // Make sure the gap doesn't start on a gap.
      //
      if (p > 0 and refPositions[p-1] != -1) {
        insertions.push_back(InsertedString(insertedString, refPositions[p-1], queryPositions[p-1], p) );
      }
      p = insEnd;
    }
    else {
      p++;
    }
  }
}


#endif
