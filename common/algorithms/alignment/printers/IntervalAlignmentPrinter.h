#ifndef INTERVAL_ALIGNMENT_PRINTER_H_
#define INTERVAL_ALIGNMENT_PRINTER_H_

#include "../../../datastructures/alignment/AlignmentCandidate.h"
#include "../../../FASTQSequence.h"

class IntervalAlignmentPrinter {
public:
    static void Print(T_AlignmentCandidate &alignment, 
                      ostream &outFile) {
        int mapQV = alignment.mapQV;
        int lastBlock = alignment.blocks.size()-1;
        if (lastBlock < 0) return;

        outFile << alignment.qName << " " << alignment.tName << " " 
                << alignment.score << " " 
                << alignment.pctSimilarity << " " 
                << alignment.qStrand << " " 
                << alignment.QAlignStart() << " "
                << alignment.QAlignEnd() << " "
                //<< alignment.qAlignedSeqPos + alignment.blocks[0].qPos << " "
                //<< (alignment.qAlignedSeqPos 
                //    + alignment.blocks[lastBlock].qPos 
                //    + alignment.blocks[lastBlock].length) << " "
                << alignment.qLength << " "
                << alignment.tStrand << " "
                << alignment.TAlignStart() << " "
                //<< alignment.tAlignedSeqPos + alignment.blocks[0].tPos << " "
                << (alignment.tAlignedSeqPos + alignment.tPos
                    + alignment.blocks[lastBlock].tPos 
                    + alignment.blocks[lastBlock].length) << " "
                << alignment.tLength << " " 
                << mapQV << endl;
                //Remove the last four fields from m4 format.
                //<< " " << alignment.nCells << " " << alignment.clusterScore 
                //<< " " << alignment.probScore << " " 
                //<< alignment.numSignificantClusters 
    } 

    // Print an alignment from a sam file in Interval (m 4) format.
    static void PrintFromSAM(AlignmentCandidate<> &alignment, 
                             ostream &outFile) {
        int mapQV = alignment.mapQV;
        int lastBlock = alignment.blocks.size()-1;
        if (lastBlock < 0) return;

        outFile << alignment.qName << " " << alignment.tName << " " 
                << alignment.score << " " 
                << alignment.pctSimilarity << " " 
                << alignment.qStrand << " " 
                << alignment.QAlignStart() << " "
                << alignment.QAlignEnd() << " "
                << alignment.qLength << " "
                << alignment.tStrand << " ";

        DNALength tS = alignment.TAlignStart(); 
        DNALength tE = alignment.tAlignedSeqPos + alignment.tPos
                       + alignment.blocks[lastBlock].tPos 
                       + alignment.blocks[lastBlock].length;

        if (alignment.tStrand == 1) {
            // Since the alignment is from a SAM file and the reference 
            // is reverse, compute tS and tE in the coordinate of the reverse
            // complementary sequence
            DNALength tmp = tS;
            tS = alignment.tLength - tE;
            tE = alignment.tLength - tmp;
        }
        outFile << tS << " "
                << tE << " "
                << alignment.tLength << " " 
                << mapQV << endl;
    }

    static void PrintHeader(ostream &out) {
        out << "qname tname score pctsimilarity qstrand qstart qend qseqlength tstrand tstart tend tseqlength mapqv" << endl;
        //ncells clusterScore probscore numSigClusters" << endl;
    }
};


#endif
