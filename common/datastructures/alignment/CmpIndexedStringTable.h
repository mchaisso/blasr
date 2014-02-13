#ifndef DATASTRUCTURES_ALIGNMENT_CMP_INDEXED_STRING_TABLE_H_
#define DATASTRUCTURES_ALIGNMENT_CMP_INDEXED_STRING_TABLE_H_

#include <vector>
#include <string>
#include <map>
using namespace std;

class CmpIndexedStringTable {
 public:

	void resize(int size) {
		names.resize(size);
		ids.resize(size);
	}

	void StoreArrayIndexMap() {
		int i;
		for (i = 0; i < ids.size(); i++) {
			idToArrayIndex[ids[i]] = i;
		}
	}

    //
    // The terminology of Index here is confusing.
    // Actually 'Index' is equivalent to 'id'. 
    // Each id represents index of an indexed string.
    // That's why an id is called an 'Index' in this function.
    // GetNameAtIndex returns name of an indexed string, 
    // whose index is the given value
    //
	bool GetNameAtIndex(int index, string &name) {
		map<int,int>::iterator mapIt;
		mapIt = idToArrayIndex.find(index);
		if (mapIt != idToArrayIndex.end()) {
			name = names[mapIt->second];
			return true;
		}
		else {
			return false;
		}
	}

    //
    // Here 'Id' means indexes of indexed strings, 
    // 'Index' means index of an 'Id' in ids
    //
	bool GetIndexOfId(int id, int &index) {
		map<int,int>::iterator mapIt;
		mapIt = idToArrayIndex.find(id);
		if (mapIt != idToArrayIndex.end()) {
			index = mapIt->second;
			return true;
		}
		else {
			return false;
		}
	}
	vector<int> ids;
	vector<string> names;
	map<int,int> idToArrayIndex;
};




#endif
