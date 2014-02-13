#ifndef SAME_MOLECULE_H_
#define SAME_MOLECULE_H_

#include <set>
#include "Types.h"

using namespace std;

class Read;

class SameMolecule {
 public:
  typedef set<Read*> ReadSet;
  ReadSet subreads;
  UInt length;

  void Write(ofstream &outFile) {
    outFile.write((const char*) &length, sizeof(length));
    UInt nSubreads = subreads.size();
    outFile.write((const char*) &nSubreads, sizeof(nSubreads));
    ReadSet::iterator subreadIt, subreadEnd;
    subreadEnd = subreads.end();
    for (subreadIt = subreads.begin();
         subreadIt != subreadEnd;
         ++subreadIt) {
      
    }
  }
    
};


#endif
