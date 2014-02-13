#ifndef SUMMARY_ALIGNMENT_PRINTER_H_
#define SUMMARY_ALIGNMENT_PRINTER_H_

#include "../../../datastructures/alignment/AlignmentCandidate.h"
#include "../../../FASTQSequence.h"
class SummaryAlignmentPrinter {
 public:
  static void Print(AlignmentCandidate<DNASequence,FASTQSequence> &alignment, ostream &outFile) {
    int lastBlock = alignment.blocks.size()-1;
    outFile << alignment.qName << " "
            << alignment.tName << " " 
            << alignment.qStrand << " " 
            << alignment.tStrand << " " 
            << alignment.score << " " 
            << alignment.pctSimilarity << " "
            << alignment.tAlignedSeqPos + alignment.blocks[0].tPos << " " 
            << alignment.tAlignedSeqPos + alignment.blocks[lastBlock].tPos + alignment.blocks[lastBlock].length << " " 
            << alignment.tLength << " "
            << alignment.qAlignedSeqPos + alignment.blocks[0].qPos << " " 
            << alignment.qAlignedSeqPos + alignment.blocks[lastBlock].qPos + alignment.blocks[lastBlock].length << " " 
            << alignment.qLength << " " << alignment.clusterWeight << endl;
  }
  static void PrintHeader(ostream &out) {
      out << "qname tname qstrand tstrand score pctsimilarity tstart tend tlength qstart qend qlength ncells" << endl;
  }
};


#endif
