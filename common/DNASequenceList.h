#ifndef DNA_SEQUENCE_LIST_H_
#define DNA_SEQUENCE_LIST_H_

#include <vector>
#include <DNASequence>

class DNASequenceList {
 public:
	vector<DNASequence*> seqPtrList;
	vector<int> cumLenghts;
	void InitFromList(vector<DNASequence> *seqList) {
		int i;
		if (seqList->size() == 0)
			return;

		for (i = 0; i < seqList.size(); i++) {
			seqPtrList[i] = &(*seqList)[i];
		}

		//
		// Store n+1 cumulative lengths, from 0 to N, where
		// N is the sum of all sequence lengths.
		//
		cumLengths.push_back(0);
		for (i = 1; i < seqPtrList.size(); i++) {
			cumLengths.push_back(cumLengths[i-1] + seqList[i-1].length);
		}
		cumLengths.push_back(cumLengths[i-1] + seqList[i-1].length);
	}

	int ListIndexToSeqIndex(int index, int &seq, int &seqIndex) {
		int i;
		for (i =1; i < cumLengths.size(); i++) {
			if (index >= cumLengths[i-1] && index < cumLengths[i]) {
				seq = i;
				seqIndex = index - cumLengths[i-1];
				return 1;
			}
		}
		// signal that there was no valid mapping.
		return 0;
	}
	
	
	
};
			
			

#endif
