#ifndef ALGORITHMS_METAGENOMICS_RANDOM_POS_GENERATOR_H_
#define ALGORITHMS_METAGENOMICS_RANDOM_POS_GENERATOR_H_

#include "DNASequence.h"
#include "statistics/statutils.h"

#include <vector>

template<typename T_Sequence>
int FindRandomPos(vector<T_Sequence> &sequences, DNALength &seqIndex, DNALength &seqPos, DNALength seqLength=0 ) {
	vector<UInt> cumulativeLengths;
	cumulativeLengths.resize(sequences.size());
	int i;
	if (sequences.size() == 0) {
		return 0;
	}
	DNALength cumulativeLength;
	cumulativeLengths[0] = sequences[0].length;
	cumulativeLength = cumulativeLengths[0];
	for (i = 1; i < sequences.size(); i++) {
		cumulativeLengths[i] = cumulativeLength = cumulativeLengths[i-1] + sequences[i].length;
	}
	bool validPosFound = false;
	int iter = 0;
	int max_iter = 100000;
	while (validPosFound == false and iter < max_iter) {
		++iter;
    if (seqLength > cumulativeLength) {
      validPosFound = false;
      iter = max_iter;
      break;
    }
		DNALength pos = RandomUnsignedInt(cumulativeLength - seqLength);
		// Make sure this sequence fits 
		for (seqIndex = 0; seqIndex < sequences.size(); seqIndex++) {
			if (cumulativeLengths[seqIndex] > pos) break;
		}
		if (cumulativeLengths[seqIndex] - pos < seqLength) {
			continue;
		}
		UInt pi;
		if (seqIndex == 0) {
			seqPos = pos;
		}
		else {
			seqPos = pos - cumulativeLengths[seqIndex-1];
		}
		bool seqContainsN = false;
		for (pi = seqPos; pi < seqPos + seqLength; pi++) {
			if (toupper(sequences[seqIndex].seq[pi]) == 'N') {
				seqContainsN = true;
				break;
			}
		}
		if (seqContainsN) {
			continue;
		}
		else {
			validPosFound = true;
		}
	}
	if (iter == max_iter) {
		return 0;
	}
	else {
		return 1;
	}
}


#endif
