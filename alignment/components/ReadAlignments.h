#ifndef COMPONENTS_ALIGNMENT_READ_ALIGNMENTS_H_
#define COMPONENTS_ALIGNMENT_READ_ALIGNMENTS_H_

#include "Enumerations.h"
#include "SMRTSequence.h"

class ReadAlignments {
public:
  /*
    This class stores the alignments from a read.  A read may be
    aligned in several different modes:
    1. noSplitSureads - Treat the read as a unit from start to end
    2. subreads       - Align each subread independently
    3. denovo         - Only align the CCS sequence from a read
    4. allpass        - Align the de novo ccs sequences and then the
    subreads to where the denovo ccs aligned.
    5. fullpass       - Same as allpass, except using only complete
    subreads.
   
    The alignments are a raggad array of n sequences; n is 1 for cases 
    1 and 3, the number of subreads for cases 2 and 4, and the number
    of full length passes for case 5.

    A ReadAligments class must only have alignments for a single type
    of read in it.

  */

  vector<vector<T_AlignmentCandidate*> > subreadAlignments;
  vector<SMRTSequence> subreads;
  AlignMode alignMode;
  SMRTSequence read;
  int GetNAlignedSeq() {
    return subreadAlignments.size();
  }

  bool AllSubreadsHaveAlignments() {
    int i, nAlignedSeq;
    nAlignedSeq = subreadAlignments.size();
    for (i = 0; i < nAlignedSeq; i++) {
      if (subreadAlignments[i].size() == 0) {
        return false;
      }
    }
    return true;
  }

  void Clear() {
    int i;
    int nAlignedSeq;
    for (i = 0, nAlignedSeq = subreadAlignments.size(); i < nAlignedSeq; i++) {
      int nAlignments;
      int a;
      for (a = 0, nAlignments = subreadAlignments[i].size(); a < nAlignments; a++) {
        delete subreadAlignments[i][a];
      }
      subreadAlignments[i].clear();
    }

    for (i = 0, nAlignedSeq = subreads.size(); i< nAlignedSeq; i++) {
      subreads[i].FreeIfControlled();
      if (subreads[i].title != NULL) {
        delete[] subreads[i].title;
        subreads[i].title = NULL;
      }
    }

    subreadAlignments.clear();
    read.Free();
  }

  void Resize(int nSeq) {
    subreadAlignments.resize(nSeq);
    subreads.resize(nSeq);
  }
  
  void CheckSeqIndex(int seqIndex) {
    if ( seqIndex < 0 or seqIndex >= subreads.size() ) {
        cout << "ERROR, adding a sequence to an unallocated position." << endl;
        assert(0);
    }
  }

  void SetSequence(int seqIndex, SMRTSequence &seq) {
    CheckSeqIndex(seqIndex);
    subreads[seqIndex] = seq;
  }

  void AddAlignmentForSeq(int seqIndex, T_AlignmentCandidate *alignmentPtr) {
    CheckSeqIndex(seqIndex);
    subreadAlignments[seqIndex].push_back(alignmentPtr);
  }

  void AddAlignmentsForSeq(int seqIndex, vector<T_AlignmentCandidate*> &seqAlignmentPtrs) {
    CheckSeqIndex(seqIndex);
    subreadAlignments[seqIndex].insert(subreadAlignments[seqIndex].end(), seqAlignmentPtrs.begin(), seqAlignmentPtrs.end());
  }

  ~ReadAlignments() {
    read.FreeIfControlled();
  }
};

#endif
