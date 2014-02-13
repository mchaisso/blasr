#ifndef PBReadNameParser_H_
#define PBReadNameParser_H_

#include <string>
#include <sstream>

using namespace std;
class PBReadNameParser {
 public:
  // Parse read names in the format 
  // m120308_205818_42142_cunknown0_s1_p0/94/0_419
  static bool GetMoleculeNameFromReadName(string &readName, string &moleculeName) {
    int firstPos = readName.find("/");
    if (firstPos == readName.npos) {
      return false;
    }
    else {
      int secondPos = readName.find("/", firstPos+1);
      if (secondPos == readName.npos) return false;
      moleculeName.assign(readName, 0, secondPos);
      return true;
    }
  }

  static bool GetReadCoordinatesFromReadName(string &readName, 
                                             unsigned int &begin, 
                                             unsigned int &end) {
    int lastPos = readName.rfind("/");
    if (lastPos == readName.npos) {
      return false;
    }
    string coordinates;
    coordinates.assign(readName.c_str(), lastPos + 1, readName.size() - lastPos - 1);
    stringstream strm(coordinates);
    if (! (strm >> begin) ) {
      return false;
    }
    strm.get();// get _
    if (! (strm >> end) ) {
      return false;
    }
    return true;
  }

};

#endif
