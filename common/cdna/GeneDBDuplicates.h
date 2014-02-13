#ifndef GENE_DB_DUPLICATES_H_
#define GENE_DB_DUPLICATES_H_
#include <map>
#include <vector>
#include "GeneDBLocus.h"

using namespace std;

typedef map<string, vector<GeneDBLocus>* >  DuplicateMap;
class GeneDBDuplicates {
 public:
  DuplicateMap duplicates;
  vector<GeneDBLocus>* FindDuplicate(string key) {
    DuplicateMap::iterator it;
    it = duplicates.find(key);
    if (it != duplicates.end()){
      return it->second;
    }
    else {
      return NULL;
    }
  }
  
  void AddDuplicate(string key, GeneDBLocus &locus) {
    vector<GeneDBLocus>* duplicateList = FindDuplicate(key);
    if (duplicateList == NULL) {
      duplicateList = new vector<GeneDBLocus>;
      duplicates[key] = duplicateList;
    }
    duplicateList->push_back(locus);
  }

};

#endif
