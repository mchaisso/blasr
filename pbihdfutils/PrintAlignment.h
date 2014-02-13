#ifndef PRINT_ALIGNMENT_H_
#define PRINT_ALIGNMENT_H_

#include <sstream>
#include "datastructures/alignment/ByteAlignment.h"

void PrintAlignment(vector<unsigned char> &byteAlignment,
                    ostream &out=cout,
                    DNALength qStart = 0,
                    DNALength tStart = 0) {
  string   refSequence;
  string   readSequence;
                    
                          
  readSequence.resize(byteAlignment.size());
  refSequence.resize(byteAlignment.size());

				
  stringstream sstrm;
  //sstrm << alnHoleNumber << "/" << qStart << "_" << cmpFile.alnInfo.alignments[alignmentIndex].GetQueryEnd();
  sstrm << "alignment";

  readSequence.resize(byteAlignment.size());
  refSequence.resize(byteAlignment.size());

  ByteAlignmentToQueryString(&byteAlignment[0], byteAlignment.size(), &readSequence[0]);
  ByteAlignmentToRefString(&byteAlignment[0], byteAlignment.size(), &refSequence[0]);				
  string ungappedRead, ungappedRef;
  RemoveGaps(readSequence, ungappedRead);
  RemoveGaps(refSequence, ungappedRef);
  Alignment alignment;
  GappedStringsToAlignment(readSequence, refSequence, alignment);
  DNASequence qAlignedSeq, rAlignedSeq;
  qAlignedSeq.seq = (Nucleotide*) &ungappedRead[0];
  qAlignedSeq.length = ungappedRead.size();
  rAlignedSeq.seq = (Nucleotide*) &ungappedRef[0];
  rAlignedSeq.length = ungappedRef.size();
				
  sstrm << "alignment";
  alignment.qName = sstrm.str();
  StickPrintAlignment(alignment, qAlignedSeq, rAlignedSeq, out, qStart, tStart);
}

#endif
