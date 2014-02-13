#ifndef FILES_BASE_CCS_ITERATOR_H_
#define FILES_BASE_CCS_ITERATOR_H_

#include "../CCSSequence.h"

class CCSIterator {
 public:
	CCSSequence *seqPtr;
	int curPass;
	int numPasses;

  
	virtual void Initialize(CCSSequence *_seqPtr) {
		seqPtr = _seqPtr;
		curPass = 0;
		numPasses = seqPtr->passDirection.size();
	}

	virtual int GetNext(int &direction, int &startBase, int &numBases) {
		if (curPass >= numPasses) {
			return 0;
		}
		else {
			direction = seqPtr->passDirection[curPass];
			startBase = seqPtr->passStartBase[curPass];
			numBases  = seqPtr->passNumBases[curPass];
			++curPass;
			return 1;
		}
	}

	void Reset() {
		curPass = 0;
	}
  
  int GetNumPasses() {
    return numPasses;
  }
};

#endif
