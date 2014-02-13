#ifndef ALIGNMENT_SET_PROGRAM_H_
#define ALIGNMENT_SET_PROGRAM_H_

#include <string>
class Program {
 public:
  string id;
  string programName;
  string commandLine;
  string previousId;
  string programVersion;

  void StoreValues(vector<SAMKeywordValuePair> &kvPairs, int lineNumber = 0) {
    // no op for now, do this later.

  }
};


#endif
