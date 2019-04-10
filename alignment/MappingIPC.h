#ifndef MAPPING_IPC_H_
#define MAPPING_IPC_H_

#include <pthread.h>

#include "MappingParameters.h"

#include "FASTASequence.h"
#include "FASTQSequence.h"
#include "tuples/DNATuple.h"
#include "tuples/CompressedDNATuple.h"
#include "files/ReaderAgglomerate.h"
#include "datastructures/mapping/MappingMetrics.h"
#include "datastructures/tuplelists/TupleCountTable.h"
#include "datastructures/suffixarray/SuffixArrayTypes.h"
#include "datastructures/metagenome/SequenceIndexDatabase.h"
#include "datastructures/reads/RegionTable.h"
#include "datastructures/bwt/BWT.h"
/*
 * This structure contains pointers to all required data structures
 * for mapping reads to a suffix array and evaluating the significance
 * of the matches.
 */

template<typename T_SuffixArray, typename T_GenomeSequence, typename T_Tuple>
class MappingData {
 public:
	T_SuffixArray        *suffixArrayPtr;
	BWT                  *bwtPtr;
	T_GenomeSequence     *referenceSeqPtr;
	SequenceIndexDatabase<FASTASequence> *seqDBPtr;
	TupleCountTable<T_GenomeSequence, T_Tuple> *ctabPtr;
	MappingParameters     params;
	MappingMetrics        metrics;
	RegionTable          *regionTablePtr;
	ReaderAgglomerate    *reader;
	ostream *outFilePtr;
	ostream *unalignedFilePtr;
  ostream *anchorFilePtr;
  ostream *clusterFilePtr;
  ostream *lcpBoundsOutPtr;
  
  // Declare a semaphore for blocking on reading from the same hdhf file.
	
  void ShallowCopySuffixArray(T_SuffixArray &dest) {
		dest.index              = suffixArrayPtr->index;
		dest.length             = suffixArrayPtr->length;
		dest.target             = suffixArrayPtr->target;
		dest.startPosTable      = suffixArrayPtr->startPosTable;
		dest.endPosTable        = suffixArrayPtr->endPosTable;
		dest.lookupTableLength  = suffixArrayPtr->lookupTableLength;
		dest.lookupPrefixLength = suffixArrayPtr->lookupPrefixLength;
		dest.tm                 = suffixArrayPtr->tm;
		dest.deleteStructures   = false;
		//		dest.useLCPTable        = suffixArrayPtr->useLCPTable;
	}
	
	void ShallowCopySequenceIndexDatabase(SequenceIndexDatabase<FASTQSequence> &dest) {
		dest.nSeqPos     = seqDBPtr->nSeqPos;
		dest.seqStartPos = seqDBPtr->seqStartPos;
		dest.nameLengths = seqDBPtr->nameLengths;
		dest.names       = seqDBPtr->names;
		dest.startPos    = seqDBPtr->startPos;
		dest.endPos      = seqDBPtr->endPos;
		dest.deleteStructures = false;
	}

	void ShallowCopyTupleCountTable(	TupleCountTable<T_GenomeSequence, T_Tuple> &dest) {
		dest.countTable       = ctabPtr->countTable;
		dest.countTableLength = ctabPtr->countTableLength;
		dest.nTuples          = ctabPtr->nTuples;
		dest.tm               = ctabPtr->tm;
		dest.deleteStructures = false;
	}

	void ShallowCopyReferenceSequence(T_GenomeSequence &refSeq) {
		refSeq.ShallowCopy(*referenceSeqPtr);
		refSeq.deleteOnExit = false;
	}

	void Initialize(T_SuffixArray *saP, T_GenomeSequence *refP, 
									SequenceIndexDatabase<FASTASequence> *seqDBP,
									TupleCountTable<T_GenomeSequence, T_Tuple> *ctabP,
									ReverseCompressIndex *rciP,
									MappingParameters &paramsP,
									ReaderAgglomerate *readerP,
									RegionTable *regionTableP,
									ostream *outFileP,
									ostream *unalignedFileP,
                  ostream *anchorFilePtrP, 
                  ostream *clusterFilePtrP=NULL) {
		suffixArrayPtr     = saP;
		referenceSeqPtr    = refP;
		seqDBPtr           = seqDBP;
		ctabPtr            = ctabP;
		regionTablePtr     = regionTableP;
		params             = paramsP;
		reader             = readerP;
		outFilePtr         = outFileP;
		unalignedFilePtr   = unalignedFileP;
    anchorFilePtr      = anchorFilePtrP;
    clusterFilePtr= clusterFilePtrP;
	}
};

#endif
