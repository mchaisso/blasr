#ifndef ALIGNMENT_COMPONENTS_ASSIGN_REF_CONTIG_LOCATION_H_
#define ALIGNMENT_COMPONENTS_ASSIGN_REF_CONTIG_LOCATION_H_

#include "datastructures/alignment/AlignmentCandidate.h"
#include "FASTQSequence.h"
#include "DNASequence.h"
#include "datastructures/metagenome/SequenceIndexDatabase.h"
#include "Types.h"

void AssignRefContigLocation(T_AlignmentCandidate &alignment, SequenceIndexDatabase<FASTQSequence> &seqdb, DNASequence &genome) {
    //
    // If the sequence database is used, the start position of
    // the alignment is relative to the start of the chromosome,
    // not the entire index.  Subtract off the start position of
    // the chromosome to get the true position.
    //
  DNALength forwardTPos;
  int seqDBIndex;
  if (alignment.tStrand == 0) {
    forwardTPos = alignment.tAlignedSeqPos;
    seqDBIndex = seqdb.SearchForIndex(forwardTPos);
    alignment.tAlignedSeqPos -= seqdb.seqStartPos[seqDBIndex];
  }
  else {
    //
    // Flip coordinates into forward strand in order to find the boundaries 
    // of the contig, then reverse them in order to find offset.
    //

    // Find the reverse complement coordinate of the index of the last aligned base.
    assert(alignment.tAlignedSeqLength > 0);
    forwardTPos = genome.MakeRCCoordinate(alignment.tAlignedSeqPos + alignment.tAlignedSeqLength - 1);
    seqDBIndex  = seqdb.SearchForIndex(forwardTPos);

    
    //
    // Find the reverse comlement coordinate of the last base of this
    // sequence.  This would normally be the start of the next contig
    // -1 to get the length, but since an 'N' is added between every
    // pair of sequences, this is -2.
    //
    DNALength reverseTOffset;
    reverseTOffset = genome.MakeRCCoordinate(seqdb.seqStartPos[seqDBIndex+1]-2);
    alignment.tAlignedSeqPos -= reverseTOffset;
  }
}

void AssignRefContigLocations(vector<T_AlignmentCandidate*> &alignmentPtrs, SequenceIndexDatabase<FASTQSequence> &seqdb, DNASequence &genome) {
  
  UInt i;
  for (i = 0; i < alignmentPtrs.size(); i++) {
    T_AlignmentCandidate *aref = alignmentPtrs[i];
    AssignRefContigLocation(*aref, seqdb, genome);
  }
}


#endif
