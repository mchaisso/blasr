#ifndef UTILS_CHANGELIST_ID_H_
#define UTILS_CHANGELIST_ID_H_

#include <string>
#include <vector>
#include <stdlib.h>
#include "StringUtils.h"
using namespace std;

class ChangeListID {
 public:
  string idString;
  vector<string> strVer;
  vector<int>    intVer;
  ChangeListID() {}
  ChangeListID(string &idStringP) {
    StoreString(idStringP);
  }
  
  void StoreString(string &idStringP) {
    idString = idStringP;
    Tokenize(idString, ".", strVer);
    int i;
    intVer.resize(strVer.size());
    for (i = 0; i < strVer.size(); i++) {
      intVer[i] = atoi(strVer[i].c_str());
    }
  }
    
  int LessThan(ChangeListID &rhs, int depth = 0) {
    if (depth == 0) {
      depth = min(intVer.size(), rhs.intVer.size());
    }
    int i;
    for (i = 0; i < depth; i++) {
      if (intVer[i] != rhs.intVer[i]) {
        return intVer[i] < rhs.intVer[i];
      }
    }
    return 0; // making it here they are equal
  }
};

void AppendPerforceChangelist(string perforceVersionString, string &version) {
  if (perforceVersionString.size() > 12) {
    version.insert(version.size(), ".");
    version.insert(version.size(), perforceVersionString, 9, perforceVersionString.size() - 11);
  }
}

#endif
