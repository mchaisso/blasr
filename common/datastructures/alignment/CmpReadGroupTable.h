#ifndef DATASTRUCTURES_ALIGNMENT_CMP_READ_GROUP_TABLE_H_
#define DATASTRUCTURES_ALIGNMENT_CMP_READ_GROUP_TABLE_H_

#include <vector>
#include <string>
using namespace std;

class CmpReadGroupTable {
 public:
	void resize(int size) {
		readGroupNameIds.resize(size);
		readGroupNames.resize(size);
	}
	vector<int>    readGroupNameIds;
	vector<string> readGroupNames;
	int lastRow;
};




#endif
