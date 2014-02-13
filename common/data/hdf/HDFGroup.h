#ifndef DATA_HDF_HDF_GROUP_H_
#define DATA_HDF_HDF_GROUP_H_

#include <vector>
#include <iostream>
#include <stdlib.h>
#include "H5Cpp.h"
#include "HDFAttributable.h"
#include "../../utils/StringUtils.h"

using namespace H5;
using namespace std;

class HDFGroup : public HDFAttributable {
 public:
	vector<string> objectNames;
	string objectName;
	Group group;
	bool  groupIsInitialized;

 HDFGroup() : HDFAttributable() {
		groupIsInitialized = false;
	}
  
	void AddGroup(string groupName) {
		group.createGroup(groupName);
		return;
	}

  H5Object* GetObject() {
    return &group;
  }

	int Initialize(CommonFG &fg, string groupName){ 
		try {
			group   = fg.openGroup(groupName.c_str());
			groupIsInitialized = true;
			return 1;
		}
		catch(FileIException &e) {
			return 0;
		}
		catch(GroupIException &e) {
			return 0;
		}
		catch(Exception &e ) {
			return 0;
		}
		return 1;
	}

  int Initialize(HDFGroup & parentGroup, string groupName) {
    return Initialize(parentGroup.group, groupName);
  }

	bool ContainsObject(string queryObjectName) {
    hsize_t objIdx;
    int numGroupObjs = group.getNumObjs();
    for (objIdx = 0; objIdx < numGroupObjs; objIdx++) {
      H5std_string groupObjectName;
      size_t objNameSize;
      groupObjectName = group.getObjnameByIdx(objIdx);
      if (groupObjectName == queryObjectName) {
        return true;
      }
    }
    return false;
  }

	void Close() {
		if (groupIsInitialized) {
			group.close();
		}
	}
};


#endif
