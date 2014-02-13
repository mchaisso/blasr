#ifndef CCS_SEQUENCE_H_
#define CCS_SEQUENCE_H_

#include "SMRTSequence.h"
#include "Enumerations.h"
#include "utils/VectorUtils.h"
//
// A CCS Sequence is both a SMRTSequence itself, and contains a list of SMRTSequences.
//
class CCSSequence : public SMRTSequence {
 public:
	UInt numPasses;
	UInt numConsensusBases;
	vector<DNALength> passStartPulse, passNumPulses, passStartBase, passNumBases;
	vector<Byte> passDirection;
	vector<Byte>      adapterHitBefore, adapterHitAfter, adapterHitConfidence;
	//
	// The CCS Sequence originates from a full length read.  That read
	// is stored here for reference later on.  The ccs read is stored in
	// the inherited fields from SMRTSequence so that it may be worked
	// with as if it were a normal non-ccs sequence.
	//
	SMRTSequence      unrolledRead;
	void Free() {
        numPasses = 0;
        numConsensusBases = 0;
		SMRTSequence::Free();
		unrolledRead.Free();
        /*
        ClearMemory(passStartPulse);
        ClearMemory(passNumPulses);
        ClearMemory(passStartBase);
        ClearMemory(passNumBases);
        ClearMemory(passDirection);
        ClearMemory(adapterHitBefore);
        ClearMemory(adapterHitAfter);
        ClearMemory(adapterHitConfidence);
        */
	}
	int GetStorageSize() {
		return SMRTSequence::GetStorageSize() + unrolledRead.GetStorageSize();
	}
	//
	// In the first iteration, Explode simply pulls the subreads out
	// that are used in the ccs.   Eventually, it will pull out all
	// high-quality subreads.
	// 
	void Explode(vector<SMRTSequence> &subreads) {
		subreads.resize(numPasses);
		int subreadIndex;
		for (subreadIndex = 0; subreadIndex < numPasses; subreadIndex++) {
			subreads[subreadIndex].ReferenceSubstring(this->unrolledRead, passStartBase[subreadIndex], passNumBases[subreadIndex]);
		}
	}
};

#endif
