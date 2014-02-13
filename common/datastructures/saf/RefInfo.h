#ifndef DATASTRUCTURES_SAF_REF_INFO
#define DATASTRUCTURES_SAF_REF_INFO

#include <vector>

using namespace std;

class OneRefInfo {
 public:
  string fullName;
  unsigned int id;
  unsigned int length;
  string md5;
  OneRefInfo() {
    fullName = md5 = "";
    id = length = 0;
  }
};


class RefInfo {
 public:
	vector<string> fullName;
	vector<uint32_t> id;
	vector<uint32_t> length;
	vector<string> md5;
  bool RefIdToIndex(uint32_t qid, int &index) {
    /*
     * Translate from external id to position in arry.
     */
    int i;
    for (i = 0; i < id.size(); i++) {
      if (id[i] == qid) {
        index = i;
        return true;
      }
    }
    return false;
  }
};
#endif
