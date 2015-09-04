#ifndef DATA_HDF_HDF_ATTRIBUTABLE_H_
#define DATA_HDF_HDF_ATTRIBUTABLE_H_

#include "H5Cpp.h"
#include <string>
#include <vector>
#include <assert.h>

using namespace std;
using namespace H5;

template<typename T_H5Object>
void CallStoreAttributeName(T_H5Object &obj, string attrName, void *attrList){ 
	((vector<string>*)attrList)->push_back(attrName);
 }

class HDFAttributable {
 public:
	vector<string> attributeNameList;

	void StoreAttributeNames(H5Object &thisobject, vector<string> &attributeNames) {
		void *destAndData[2];
		int nAttr = thisobject.getNumAttrs();
		unsigned int bounds[2];
		bounds[0] = 0;
		bounds[1] = nAttr;
		attributeNameList.clear();
		thisobject.iterateAttrs(&CallStoreAttributeName, 
														bounds, (void*) &attributeNames);
	}

  virtual H5Object* GetObject() {
    return NULL;
  }

	int ContainsAttribute(string attributeName) {
		int i;
    vector<string> tmpAttributeNames;
    H5Object *obj = GetObject();
    assert(obj != NULL);
    StoreAttributeNames(*obj, tmpAttributeNames);
		for (i = 0; i < tmpAttributeNames.size(); i++) {
			if (tmpAttributeNames[i] == attributeName) return true;
		}
		return false;
	}

};


#endif
