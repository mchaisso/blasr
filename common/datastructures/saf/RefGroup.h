#ifndef DATASTRUCTURES_SAF_REF_GROUP_H_
#define DATASTRUCTURES_SAF_REF_GROUP_H_

#include <string>
#include <vector>

#include <stdint.h>

using namespace std;
class RefGroup {
 public:
	vector<uint32_t> id;
	vector<string> path;
	vector<string> refGroupName;
	vector<uint32_t> refInfoId;
  bool IdToIndex(int idKey, int &idIndex) {
    int i;
    for (i = 0; i < refInfoId.size(); i++) {
      if (refInfoId[i] == idKey) {
        idIndex = i; return true;
      }
    }
    return false;
  }

	int FindPath(int idKey, string &pathVal, string &groupNameVal) {
		int i;
		for (i = 0; i < id.size(); i++) {
			if (id[i] == idKey) {
				pathVal = path[i];
				groupNameVal = refGroupName[i];
				return 1;
			}
		}
		return 0;
	}
};

#endif
