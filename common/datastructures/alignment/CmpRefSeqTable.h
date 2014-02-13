#ifndef DATASTRUCTURES_ALIGNMENT_CMP_REF_SEQ_TABLE_H_
#define DATASTRUCTURES_ALIGNMENT_CMP_REF_SEQ_TABLE_H_

#include <vector>
#include <string>
using namespace std;

class CmpRefSeqTable {
 public:
	void resize(int size) {
		refSeqNameIds.resize(size);
		refSeqNames.resize(size);
	}
	vector<int> refSeqNameIds;
	vector<string> refSeqNames;
	int lastRow;
};




#endif
