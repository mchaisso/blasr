#ifndef ALIGNMENT_TYPE_H_
#define ALIGNMENT_TYPE_H_

#include <vector>
#include <string>
#include <stdlib.h>
using namespace std;
#include "DNASequence.h" // to define DNALength
#include "AlignmentMap.h"
#include "AlignmentBlock.h"
#include "AlignmentGapList.h"
#include "Path.h"
#include "AlignmentStats.h"



class Alignment : public AlignmentStats {
 public:

	// the FASTA titles of each sequence
	string qName, tName;

	// Strands represented in the alignment, 0=forward, 1=reverse
	int qStrand, tStrand;

	// The starting pos in the text and query of the start of the 
	// alignment, in the window that is matched.
	DNALength qPos, tPos;
	DNALength  qAlignLength;
  DNALength tAlignLength;
	DNALength qLength;
	DNALength tLength;

	double probability;
	float zScore;
	float probScore;
	int   sumQVScore; 
	int   nCells;
  int   nSampledPaths;
	vector<Block> blocks;
  vector<GapList> gaps;

  void CopyStats(Alignment &rhs) {
    AlignmentStats::CopyStats((AlignmentStats) rhs);
    probability = rhs.probability;
    probScore   = rhs.probScore;
    sumQVScore  = rhs.sumQVScore;
    nCells      = rhs.nCells;
    nSampledPaths = rhs.nSampledPaths;
  }
	// 
	// The position in the query is qPos + block[i].qPos
	// and the position in the text is tPos + block[i].tPos
	//
	void Clear() {
		qName = "";
		tName = "";
		blocks.clear();
    gaps.clear();
	}

	Alignment& operator=(const Alignment &rhs) {
		qName = rhs.qName;
		tName = rhs.tName;
		qStrand = rhs.qStrand; tStrand = rhs.tStrand;
		qPos = rhs.qPos; tPos = rhs.tPos;
		qAlignLength = rhs.qAlignLength; tAlignLength = rhs.tAlignLength;
		qLength = rhs.qLength;
    tLength = rhs.tLength;
		zScore = rhs.zScore;
		blocks.clear();
		blocks = rhs.blocks;
    gaps.clear();
    gaps   = rhs.gaps;
		nCells = rhs.nCells;
    nSampledPaths = rhs.nSampledPaths;
    AlignmentStats::Assign(rhs);
		return *this;
	}

 Alignment() : AlignmentStats() {
		qName = ""; tName = "";
		qStrand = tStrand = 0;
		qPos = tPos = 0;
		qLength = tLength = 0;

		probability   = 0;
		nCells        = 0;
    nSampledPaths = 0;
    zScore        = 0;
    score = probScore = sumQVScore = 0;
    qAlignLength  = tAlignLength = 0;
	}
	
	unsigned int size() {
		return blocks.size();
	}

	void Assign(Alignment &rhs) {
		((AlignmentStats*)this)->Assign(rhs);
		qPos = rhs.qPos; tPos = rhs.tPos;
		qAlignLength = rhs.qAlignLength;
		tAlignLength = rhs.tAlignLength;
		qLength = rhs.qLength;
		zScore  = rhs.zScore;
		qName   = rhs.qName;
		tName   = rhs.tName;
		qStrand = rhs.qStrand;
		tStrand = rhs.tStrand;
		nCells  = rhs.nCells;
		vector<Block> empty;
		blocks.swap(empty);
		blocks.resize(rhs.size());
		int b;
		for (b = 0; b < rhs.blocks.size(); b++) {
			blocks[b].Assign(rhs.blocks[b]);
		}
	}
  int ComputeNumAnchors(int minAnchorSize, int &nAnchors, int &nAnchorBases) {
    int i;
    nAnchors = 0;
    nAnchorBases = 0;
    for (i = 0; i < blocks.size(); i++) {
      if (blocks[i].length >= minAnchorSize) {
        nAnchors++;
        nAnchorBases += blocks[i].length;
      }
    }
    return nAnchors;
  }

	void AllocateBlocks(int nBlocks) {
		blocks.resize(nBlocks);
	}

  void AppendAlignmentGaps(Alignment &next, bool mergeFirst=false) {
    //
    // append all gaps belonging to ext to the gap list in this 
    // alignment.  The logic for determining just what the gap should
    // be between two separate alignments is somewhat complicated, and
    // should probably be handled outside this function.  The default
    // to do this is to skip the first gap in the next alignment, and
    // assume that the code has handled computing the gap stored in
    // the last gap in this alignment correctly. 
    //
    // At your own risk, you can simply merge the first gap in next
    // with the last gap in the current alignment.

    assert(gaps.size() > 0);
    
    // Ditto for the next alignment.
    assert(next.gaps.size() > 0);

    vector<GapList>::iterator secondGapIt = next.gaps.begin();
    if (mergeFirst) {
      // Merge the gap list in the first gap in next sequence to the
      // last gap in the current sequence.
      gaps[gaps.size()-1].insert(gaps[gaps.size()-1].end(), secondGapIt->begin(), secondGapIt->end());
    }
    // Append all other gaps to the current one.
    secondGapIt++;
    gaps.insert(gaps.end(), secondGapIt, next.gaps.end());
  }

  void AppendAlignmentBlocks(Alignment &next, int qOffset = 0, int tOffset = 0) {
		VectorIndex n;
		Block tempBlock;
		for (n = 0; n < next.blocks.size(); n++ ) {
			tempBlock = next.blocks[n];
			tempBlock.qPos += qOffset;
			tempBlock.tPos += tOffset;
			blocks.push_back(tempBlock);
		}
  }

	void AppendAlignment(Alignment &next) {
		int qOffset = next.qPos - qPos;
		int tOffset = next.tPos - tPos;
    AppendAlignmentBlocks(next, qOffset, tOffset);
	}
  
  /*
    Transform the series of operations in an optimal dynamic
    programming path to a block representation of the alignment.

    Since it is possible to have an adjacent insertion and deletion,
    the gap blocks are tracked in addition to the match blocks.
  */

	void ArrowPathToAlignment(vector<Arrow> &optPath) {
		int q, t;
		VectorIndex a = 1;
		q = 0; t = 0;
		Block b;
		a = 0;
    bool beforeFirstBlock = true;
		while (a < optPath.size()) {
      //
      // Allow for there to be a block at the beginning
      // of the alignment, so process gap characters first.
      //
      if (!beforeFirstBlock) {
        if (optPath[a] == Diagonal) {
          // Start of a block;
          b.qPos = q;
          b.tPos = t;
          b.length = 0;
          while(a < optPath.size() and optPath[a] == Diagonal) {
            b.length++;
            a++;
            t++;
            q++;
          }
          blocks.push_back(b);
        }
      }
      gaps.push_back(GapList());
      int curGapList = gaps.size() - 1;
      //
      // Add gaps as condensed blocks of insertions or deletions.  It
      // is possible there are multiple stretches of ins/del/ins/del
      // patterns, so have to loop over all of these.
      //
      while (a < optPath.size() and 
             (optPath[a] == Left or 
              optPath[a] == Up)) {
        if (a < optPath.size() and optPath[a] == Left) {
          int gapStart = a;
          while(a < optPath.size() and optPath[a] == Left) {
            t++;
            a++;
          }
          gaps[curGapList].push_back(Gap(Gap::Query, a - gapStart));
          continue;
        }
        else if (a < optPath.size() and optPath[a] == Up) {
          int gapStart = a;
          while(a < optPath.size() and optPath[a] == Up) {
            q++;
            a++;
          }
          gaps[curGapList].push_back(Gap(Gap::Target, a - gapStart));
          continue;
        }
      }
      if (a == optPath.size()) {
        if (gaps.size() > 0) {
          gaps[curGapList].clear();
        }
      }
      assert(a == optPath.size() or gaps[curGapList].size() != 0 or beforeFirstBlock == true);
      beforeFirstBlock = false;
		}
	}

  //
  // samtools / picard do not like the pattern
  // insertion/deletion/insertion (or the opposite).  To get around
  // this, reorder the idi patterns to iid (or did to idd).  This
  // produces the same scoring alignment, however it is reordered so
  // that Picard / samtools accepts the alignments.
  //
  void OrderGapsByType() {
    int g;
    //
    // Gaps at the beginning and end of the sequence are hard to deal
    // with. Just get rid of them. 
    //
    RemoveEndGaps();
    //
    // Start at 1 since the gaps at the beginning of the sequence are
    // removed.
    //
    for (g = 1; g < gaps.size(); g++ ) {
      if (gaps[g].size() <= 1) {
        continue;
      }
      Gap queryGap, targetGap;
      GapList condensedGapList;
      int gi;
      targetGap.seq = Gap::Target;
      queryGap.seq  = Gap::Query;
      for (gi = 0; gi < gaps[g].size(); gi++) {
        if (gaps[g][gi].seq == Gap::Target) {
          targetGap.length += gaps[g][gi].length;
        }
        else {
          queryGap.length += gaps[g][gi].length;
        }
      }
      int nTypes = 0;
      gaps[g].clear();
      int matchExtend = 0;
      if (targetGap.length > queryGap.length) {
        targetGap.length -= queryGap.length;
        gaps[g].push_back(targetGap);
        matchExtend = queryGap.length;
      }
      else if (queryGap.length > targetGap.length) {
        queryGap.length -= targetGap.length;
        gaps[g].push_back(queryGap);
        matchExtend = targetGap.length;
      }
      else {
        matchExtend = targetGap.length;
      }

      if (matchExtend > 0) {
        assert(g>0);
        blocks[g-1].length += matchExtend;
      }
      //
      // When targetGap.length == queryGap.length, there is no gap, so
      // just leave the gap list cleared.
      //
    }
  }

  //
  // Transform an alignment that has up to one long gap in it to a
  // block based alignment.

	void LongGapArrowPathToAlignment(vector<Arrow> &optPath, DNALength lengthOfLongGap)  {

    DNALength i;
    int numLongGaps = 0;

    // Input checking.
    // Only one long gap is allowed per alignment.  Make sure this is
    // the case on the input.
    //
    for (i = 0; i < optPath.size(); i++) {
      if (optPath[i] == AffineLongDelLeft or
          optPath[i] == AffineLongDelClose) {
        numLongGaps++;
      }
    }
    
    if (numLongGaps > 1) {
      cout << "ERROR. Only one long gap per alignment is allowed." << endl;
      exit(1);
    }

    //
    // First locate the position of the gap.
    // Also, change the gap to a normal arrow.
    //
    DNALength indexOfLongGap; // undefined until one is found
    bool aLongGapWasFound = false;
    Arrow longGapArrow; // will hold the type of the long gap
    int numBlocksBeforeGap = 0;
    int indexOfLastMatchBeforeGap = 0;
    // Now locate both the type of type of the arrow, and the 
    // position 
    for (i = 0; i < optPath.size(); i++) {
      //
      // Count the number of blocks. This will tell us where to insert
      // the gap.
      //
      if ( i > 0 and 
           optPath[i-1] == Diagonal and 
           optPath[i] != Diagonal ) {
        numBlocksBeforeGap++;
        indexOfLastMatchBeforeGap = i;
      }
      //
      // Look for the gap.
      //
      if (optPath[i] == AffineLongDelLeft or
          optPath[i] == AffineLongDelClose) {
        aLongGapWasFound = true;
        longGapArrow = optPath[i];
        indexOfLongGap = i;
        optPath[i] = Left;
        break;
      }
    }

    //
    // Next transform the path into an alignment that lacks the gap. 
    //
    ArrowPathToAlignment(optPath);

    //
    // Finally, insert the gap into the block form of the alignment.
    //
    
    if (aLongGapWasFound and numBlocksBeforeGap < blocks.size()) {
      // Found a gap, add it.
      
      // First, find which gap corresponds to the long gap.  This is
      // the hardest step.  First, find how many arrow instructions
      // there are between the last block and the arrow.

      int numGapChars = indexOfLongGap - indexOfLastMatchBeforeGap + 1;
      int gi;
      
      // Define some variables for readability.
      int indexOfBlockBeforeGap;
      int gapIndex;
      gapIndex = numBlocksBeforeGap;

      assert(gapIndex < gaps.size());

      // There must be at least one gap here (the long deletion)
      assert(gaps[gapIndex].size() > 0);
      
      int cumulativeGapLength = 0;
      bool indexOfGapFound = false;
      for (gi = 0; gi < gaps[gapIndex].size(); gi++) {
        cumulativeGapLength += gaps[gapIndex][gi].length;
        if (cumulativeGapLength >= numGapChars) {
          // Found the gap where the long deletion happened.
          // Make sure this is on a deletion.
          assert(gaps[gapIndex][gi].seq == Gap::Query);
          indexOfGapFound = true;
          break;
        }
      }

      assert(indexOfGapFound == true);
      //
      // Found the gap corresponding to the long deletion.
      // Now, add in the length of the long gap, taking into account
      // the fact that there is already one base used in the previous
      // accounting of the gap.
      gaps[gapIndex][gi].length += lengthOfLongGap - 1;

      // 
      // Now fix the offsets of the positions of the rest of the
      // blocks in the sequence.
      //
      UInt b;
      for (b = numBlocksBeforeGap; b < blocks.size(); b++) {
        blocks[b].tPos += lengthOfLongGap - 1;
      }
    }
  }

  //
  // The length of the aligned sequence in the query.
  //
	DNALength QEnd() {
		if (blocks.size() > 0) {
			return blocks[blocks.size()-1].QEnd();
		}
		else { return 0; }
	}
  
  //
  // The lenght of the aligned sequence in the target.
  //
	DNALength TEnd() {
		if (blocks.size() > 0) {
			return blocks[blocks.size()-1].TEnd();
		}
		else { return 0; }
	}

  DNALength GenomicTBegin() {
    return tPos;
  }
  
  DNALength GenomicTEnd() {
    return tPos + TEnd();
  }
  
  //
  // Some programs do not accept alignments that have gaps at their
  // ends.  This is used to trim gaps at the ends of alignments (even
  // if the structure represents an acceptable alignment).
  //

  void RemoveEndGaps() {
    if (gaps.size() > 0 and gaps[0].size() > 0) {
      int i;
      for (i = 0; i < gaps[0].size(); i++) {
        if (gaps[0][i].seq == Gap::Target) {
          qPos += gaps[0][i].length;
        }
        else {
          tPos += gaps[0][i].length;
        }
      }
      gaps[0].clear();
    }

    if (gaps.size() > 1 ) {
      int lastGap = gaps.size() - 1;
      gaps[lastGap].clear();
    }
  }

  //
  // Make a map that takes an alignment and allows one to query for
  // what the position in a sequence is given a position in the alignment.
  //
	void AlignmentToMap(AlignmentMap &map) {
		map.tPos = tPos;
		map.qPos = qPos;
		if (blocks.size() == 0) {
			map.alignPos.resize(0);
			return;
		}

		int lastQueryPos = blocks[blocks.size()- 1].qPos + blocks[blocks.size() - 1].length;
		map.alignPos.resize(lastQueryPos - qPos + 1);
		std::fill(map.alignPos.begin(),
							map.alignPos.end(), 
							-1);

    //
    // Not done.
    cout << "AlignmentToMap is not yet impelmented." << endl;
    assert(0);
	}

	/*
		COMMENTED OUT -- To be uncommented when this is needed. It is not
	guaranteed to work, either.
	void CreateReverseComplement(Alignment &rev, DNALength queryLength) {
		int b;
		rev.AllocateBlocks(blocks.size());
		int revBlock = 0;
		int qGap, tGap;
		if (blocks.size() == 0) {
			// 
			// Do some work for this empty alignment. For now, just flag an error.
			//
			cout << "ERROR, taking the complement of an empty alignment." << endl;
			assert(0);
		}
		b = blocks.size() - 1;
		qGap = blocks[b].qPos - blocks[b-1].qPos;
		tGap = blocks[b].tPos - blocks[b-1].tPos;
		blocks[revBlock].tPos = 0;
		blocks[revBlock].qPos = 0;
		for (b = blocks.size() - 2 ; b >= 0; b--) {
			rev.blocks[revBlock].length = blocks[b].length;
			qGap = blocks[b].qPos - blocks[b-1].qPos;
			tGap = blocks[b].tPos - blocks[b-1].tPos;
			rev.blocks[revBlock].

		}
	}
	*/
};




class MatchedAlignment : public Alignment {
	//
	// This data structure holds two things: alignments, of course, and in addition
	// the coordinates of sequences that are successively refined in order to produce 
	// the alignment.  This is somewhat tricky when the target genome has been
	// transformed by some noise-reducing function phi(t). 
	//
	// Before aligning a read, it is first mapped to the genome, or transformed 
	// then mapped to the transformed genome.  Because the mapping is inexact, the
	// region a read is mapped to is typically much larger than the read.  
	// The coordinates of the mapped region are stored in tStart and tEnd
	// Alignments are performed in native nucleotide space (not transformed).
	
	// For now, the query is always qStart=0, qEnd = queryLength.


	//
	// When mapping a read to a set of concatenated chromosomes, each chromosome
	// has an offset into the file.  Therefore though a sequence may be aligned 
	// to a region starting at tStart, the relative offset into the chromosome
	// is tStart - tChromOffset.  This is used when printing the coordinates of a match.
	//
 public: 
	int refIndex;
	int readIndex;
	DNALength tStart, tEnd, qStart, qEnd;
	int tChromOffset;

	MatchedAlignment &Assign(MatchedAlignment &rhs) {
		((Alignment*)(this))->Assign(rhs);
		refIndex = rhs.refIndex;
		readIndex = rhs.readIndex;
		tStart   = rhs.tStart;
		tEnd     = rhs.tEnd;
		qStart   = rhs.qStart;
		qEnd     = rhs.qEnd;
		tChromOffset = rhs.tChromOffset;
		return *this;
	}
};


/*
 *  Create a structure for storing the information output by compare sequences.
 *  Namely, the two string representations of the alignment.
 */
class CompSeqAlignment : public Alignment {
 public:
	string tString, qString, alignString;
};

#endif
