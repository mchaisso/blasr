#ifndef REPMASK_ROW_H_
#define REPMASK_ROW_H_

#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using namespace std;

class RepmaskRow {
 public:

  //   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat
  //score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID
  // 1004  17.5  0.5  0.1  chr21     9719769 9721892 (37222431) +  ALR/Alpha      Satellite/centr          1 2132   (12)      1
  
  int swScore;
  float pctDiv, pctDel, pctIns;
  string querySeq;
  long genomeStartPos, genomeEndPos, genomeLeft;
  string matchingRepeat;
  string repeatClass;
  int  repeatStartPos, repeatEndPos, repeatStrand;
  int id;
  char strand;

  template<typename int_type>
  int ParseParentheticalNumber(string &parenNumber, int_type &result) {
    if (parenNumber.size() < 2) {
      result = atoi(parenNumber.c_str());
      return 0;
    }
    else {
      if (parenNumber[0] == '(' and parenNumber[parenNumber.size()-1] == ')') {
        result = atoi(parenNumber.c_str() + 1);
        return 1;
      }
      else {
        result = atoi(parenNumber.c_str());
        return 0;
      }
    }
  }

};

istream &operator>>(istream &in, RepmaskRow &r) {
  string genomeLeftStr, repPosBeginStr, repPosEndStr, repLeftStr;
    
  if (! (in >> r.swScore >> r.pctDiv >> r.pctDel >> r.pctIns)) return in;
  if (! (in >> r.querySeq >> r.genomeStartPos >> r.genomeEndPos >> genomeLeftStr)) return in;
  if (! (in >> r.strand) ) return in;
  if (! (in >> r.matchingRepeat >> r.repeatClass)) return in;
  if (! (in >> repPosBeginStr >> repPosEndStr >> repLeftStr)) return in;
  if (! (in >> r.id)) return in;
    
  r.ParseParentheticalNumber(genomeLeftStr, r.genomeLeft);
  r.ParseParentheticalNumber(repPosBeginStr, r.repeatStartPos);
  r.ParseParentheticalNumber(repPosEndStr, r.repeatEndPos);
  r.ParseParentheticalNumber(repLeftStr, r.repeatStrand);
    
  return in;
}


#endif
