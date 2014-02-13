#ifndef ALIGNMENT_SAM_HEADER_H_
#define ALIGNMENT_SAM_HEADER_H_

#include "SAMKeywordValuePair.h"
#include <string>

class SAMHeader {
 public:
  string formatVersion;
  enum SortingOrder {unknown, sorted , queryname, coordinate};
  SortingOrder sortingOrder;

  void StoreValues(vector<SAMKeywordValuePair> &kvPairs, int lineNumber = 0) {
    int i;
    for ( i = 0; i < kvPairs.size(); i++) {
      if (kvPairs[i].key == "VN") {
        formatVersion = kvPairs[i].value;
      }
      else if (kvPairs[i].key == "SO") {
        if (kvPairs[i].value == "unknown" ||
            kvPairs[i].value == "unsorted") {
          sortingOrder = unknown;
        }
        else if (kvPairs[i].value == "sorted") {
          sortingOrder = sorted;
        }
        else if (kvPairs[i].value == "queryname") {
          sortingOrder =queryname;
        }
        else if (kvPairs[i].value == "coordinate") {
          sortingOrder = coordinate;
        }
        else {
          cout << "Invalid sorting order " << kvPairs[i].value << " at line " << lineNumber;
        }
      }
    }
  }
};
#endif
