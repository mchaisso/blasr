#ifndef ALGORITHMS_ALIGNMENT_PRINTERS_SAMPRINTER_H_
#define ALGORITHMS_ALIGNMENT_PRINTERS_SAMPRINTER_H_

#include <sstream>
#include <stdint.h>
#include "SMRTSequence.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "datastructures/alignment/AlignmentContext.h"
#include "datastructures/alignmentset/SAMSupplementalQVList.h"

#define MULTI_SEGMENTS 0x1
#define ALL_SEGMENTS_ALIGNED 0x2
#define SEGMENT_UNMAPPED 0x4
#define NEXT_SEGMENT_UNMAPPED 0x8
#define SEQ_REVERSED 0x10
#define SEQ_NEXT_REVERSED 0x20
#define FIRST_SEGMENT 0x40
#define LAST_SEGMENT 0x80
#define SECONDARY_ALIGNMENT 0x100
#define NO_PASS_QUALITY 0x200
#define PCR_OR_OPTICAL_DUPLICATE 0x400


namespace SAMOutput {

  enum Clipping {hard, soft, subread, none};
	static string SAMVersion;
	
	template<typename T_Alignment>
  void BuildFlag(T_Alignment &alignment, AlignmentContext &context, uint16_t &flag) {

    /*
     *  Much of the flags are commented out for now since they do not
     *  generate PICARD compliant SAM.  This needs to be worked out. 
     */


    flag = 0;

    if (alignment.tStrand != alignment.qStrand) {
      flag |= SEQ_REVERSED;
    }
    else if (context.nSubreads > 1) {
      /*
       * Remember, if you're not first, you're last.
       */
    }
    if (context.isPrimary == false) {
     flag |= SECONDARY_ALIGNMENT;
    }
  }

  
  // 
  // The aligned sequence is either the sequence from the first
  // aligned base to the last (hard and no clipping), or first high
  // quality base to the last high quality base (soft clipping).
  //
  template<typename T_Alignment, typename T_Sequence>
  void SetAlignedSequence(T_Alignment &alignment, T_Sequence &read,
                          T_Sequence &alignedSeq,
                          Clipping clipping = none) {
    //
    // In both no, and hard clipping, the dna sequence that is output
    // solely corresponds to the aligned sequence.
    //
    DNALength clippedReadLength = 0;
    DNALength clippedStartPos   = 0;


    if (clipping == none or clipping == hard) {
			DNALength qStart = alignment.QAlignStart();
			DNALength qEnd   = alignment.QAlignEnd();
      clippedReadLength = qEnd - qStart;
      clippedStartPos   = qStart;
		}
		else if (clipping == soft) {
      clippedReadLength = read.length - read.lowQualityPrefix - read.lowQualitySuffix;
			clippedStartPos   = read.lowQualityPrefix;
		}
		else if (clipping == subread) {
			clippedReadLength = read.subreadEnd - read.subreadStart;
			clippedStartPos = read.subreadStart;
		}			
    else {
      cout <<" ERROR! The clipping must be none, hard, or soft when setting the aligned sequence." << endl;
      assert(0);
    }

		//
		// Set the aligned sequence according to the clipping boundaries.
		//
		if (alignment.tStrand == 0) {
			alignedSeq.ReferenceSubstring(read, clippedStartPos, clippedReadLength);
		}
		else {
			T_Sequence subSeq;
			subSeq.ReferenceSubstring(read, clippedStartPos, clippedReadLength);
			subSeq.MakeRC(alignedSeq);
			alignedSeq.deleteOnExit = true;
		}
  }



 void CreateDNAString(DNASequence &seq, 
                      DNASequence &clippedSeq,
                      //
                      // Trimming is used for both hard non-clipping
                      // so it is called trim instead of clip.
                      //
                      int trimFront=0,
                      int trimEnd=0) {
   assert(seq.length - trimEnd >= trimFront);

   clippedSeq.seq    = &seq.seq[trimFront];
   clippedSeq.length = seq.length - trimEnd - trimFront;
 }
 
 template<typename T_Alignment>
 void AddGaps(T_Alignment &alignment, int gapIndex,
              vector<int> &opSize, vector<char> &opChar, int &qPos, int &tPos) {
   int g;
   for (g = 0; g < alignment.gaps[gapIndex].size(); g++) {
     if (alignment.gaps[gapIndex][g].seq == Gap::Query) {
       opSize.push_back(alignment.gaps[gapIndex][g].length);
       opChar.push_back('D');
			 tPos += alignment.gaps[gapIndex][g].length;
     }
     else if (alignment.gaps[gapIndex][g].seq == Gap::Target) {
       opSize.push_back(alignment.gaps[gapIndex][g].length);
       opChar.push_back('I');
			 qPos += alignment.gaps[gapIndex][g].length;
     }
   }
 }
 template<typename T_Alignment>
 void AddUngappedOperations(T_Alignment &alignment, 
														int blockIndex,
														int qPos,
														int tPos,
														vector<int> &opSize, 
														vector<char> &opChar) {
	 int i;
	 int opStart;
	 i = 0;
	 string qstr, tstr;
	 qstr.insert(0, (char*) &alignment.qAlignedSeq.seq[qPos], alignment.blocks[blockIndex].length);
	 tstr.insert(0, (char*) &alignment.tAlignedSeq.seq[tPos], alignment.blocks[blockIndex].length);
	 while (i < alignment.blocks[blockIndex].length) {
		 opStart = i;
		 while (i < alignment.blocks[blockIndex].length and 
						alignment.qAlignedSeq.seq[qPos+i] != alignment.tAlignedSeq.seq[tPos+i]) {
			 i+=1;
		 }
		 if (i > opStart) {
			 opSize.push_back(i - opStart);
			 opChar.push_back('X');
		 }
		 opStart = i;
		 while (i < alignment.blocks[blockIndex].length and 
						alignment.qAlignedSeq.seq[qPos+i] == alignment.tAlignedSeq.seq[tPos+i]) {
			 i+=1;
		 }
		 if (i > opStart) {
			 opSize.push_back(i - opStart);
			 opChar.push_back('=');
		 }
	 }		 
 }
 template<typename T_Alignment>
 void AddUnmatchedOperations(T_Alignment &alignment, 
														 int qPos,
														 int tPos,
														 int alnLength,
														 vector<int> &opSize, 
														 vector<char> &opChar) {
	 int i;
	 int opStart;
	 i = 0;
	 while (i < alnLength) {
		 opStart = i;
		 while (i < alnLength and 
						alignment.qAlignedSeq.seq[qPos+i] != alignment.tAlignedSeq.seq[tPos+i]) {
			 i+=1;
		 }
		 if (i > opStart) {
			 opSize.push_back(i - opStart);
			 opChar.push_back('X');
		 }
		 opStart = i;
		 while (i < alnLength and 
						alignment.qAlignedSeq.seq[qPos+i] == alignment.tAlignedSeq.seq[tPos+i]) {
			 i+=1;
		 }
		 if (i > opStart) {
			 opSize.push_back(i - opStart);
			 opChar.push_back('=');
		 }
	 }		 
 }

 template<typename T_Alignment>
 void CreateNoClippingCigarOps(T_Alignment &alignment, 
															 int qPos,
															 int tPos,
															 vector<int> &opSize, 
															 vector<char> &opChar) {
    //
    // Create the cigar string for the aligned region of a read,
    // excluding the clipping.
    //
    int b;
    // Each block creates a match NM (N=block.length)
    int nBlocks = alignment.blocks.size();
    int nGaps   = alignment.gaps.size();

    //
    // Add gaps at the beginning of the alignment.
    //
    if (nGaps > 0) {
      AddGaps(alignment, 0, opSize, opChar, qPos, tPos);
    }
    for (b = 0; b < nBlocks; b++) {
      //
      // Determine the gap before printing the match, since it is
      // possible that the qurey and target are gapped at the same
      // time, which merges into a mismatch.
      //
      int qGap=0, tGap=0, commonGap=0;

      if (nGaps == 0) {
				qGap = 0;
				tGap = 0;
				//
				// For blocks in the middle of gaps, add insertion or deletion characters.
				//
        if (b + 1 < nBlocks) {
          qGap = alignment.blocks[b+1].qPos - alignment.blocks[b].qPos - alignment.blocks[b].length;
          tGap = alignment.blocks[b+1].tPos - alignment.blocks[b].tPos - alignment.blocks[b].length;
				}
				//
				// Add the current block to the alignment.
				//
				AddUngappedOperations(alignment, b, qPos, tPos, opSize, opChar);
				qPos += alignment.blocks[b].length;
				tPos += alignment.blocks[b].length;

				//
				// In sloppy regions, there may be a set of bases that are not
				// aligned in both the reference and query.  It is necessary
				// to account for all bases, but SAM does not allow
				// overlapping gaps.  To account for this the largest of the
				// two (insertion or deletion) is printed, minus the
				// difference.
				// Next, the extra characters must be printed in an alignment

				int unmatchedLength = 0;
				if (qGap > 0 and tGap > 0) {
					int commonGap;
					commonGap = min(qGap, tGap);
					qGap -= commonGap;
					tGap -= commonGap;
					unmatchedLength = commonGap;
				}

				if (qGap > 0 or tGap > 0) {
					if (qGap > 0) {
						opSize.push_back(qGap);
						opChar.push_back('I');
						qPos += qGap;
					}
					if (tGap > 0) {
						opSize.push_back(tGap);
						opChar.push_back('D'); 
						tPos += tGap;
					}
				}

				AddUnmatchedOperations(alignment, qPos, tPos, unmatchedLength, opSize, opChar);
				qPos += unmatchedLength;
				tPos += unmatchedLength;

      }
      else {
				AddUngappedOperations(alignment, b, qPos, tPos, opSize, opChar);
				qPos += alignment.blocks[b].length;
				tPos += alignment.blocks[b].length;
        AddGaps(alignment, b+1, opSize, opChar, qPos, tPos);
      }
    }
  }

 template<typename T_Alignment, typename T_Sequence>
  void SetSoftClip(T_Alignment &alignment,
                   T_Sequence &read,
									 DNALength hardClipPrefix,
									 DNALength hardClipSuffix,
                   DNALength &softClipPrefix, 
                   DNALength &softClipSuffix) {
    DNALength qStart, qEnd;
    qStart = alignment.QAlignStart() + alignment.blocks[0].qPos;
    qEnd   = alignment.QAlignEnd();

    assert(qStart >= hardClipPrefix);
    softClipPrefix = qStart - hardClipPrefix;
    assert(alignment.QAlignEnd() + hardClipSuffix <= read.length);
    softClipSuffix = read.length - hardClipSuffix - qEnd;
  }
 
 template<typename T_Alignment, typename T_Sequence>
  void SetHardClip(T_Alignment &alignment, 
                   T_Sequence &read,
                   DNALength &prefixClip,
                   DNALength &suffixClip) {
    //
    // Set the hard clipping assuming the read is in the forward
    // direction.
    //
    prefixClip = alignment.QAlignStart();
    suffixClip = read.length - alignment.QAlignEnd();
    if (alignment.tStrand == 1) {
      //
      // If the read is instead reverse, swap the clipping lengths.
      //
      swap(prefixClip, suffixClip);
    }
  }

  void CigarOpsToString(vector<int> &opSize,
                        vector<char> &opChar,
                        string &cigarString) {
    stringstream sstrm;
    int i, nElem;
    for (i = 0, nElem = opSize.size(); i < nElem; i++) {
      sstrm << opSize[i] << opChar[i];
    }
    cigarString = sstrm.str();
  }

  //
  // Straight forward: create the cigar string allowing some clipping
  // The read is provided to give length and hq information.
  //
  template<typename T_Alignment, typename T_Sequence>
  void CreateCIGARString(T_Alignment &alignment,
                         T_Sequence &read,
                         string &cigarString, 
												 Clipping clipping, 
												 DNALength &prefixSoftClip, DNALength &suffixSoftClip,
												 DNALength &prefixHardClip, DNALength &suffixHardClip
                         ) {
    cigarString = "";
    // All cigarString use the no clipping core
    vector<int> opSize;
    vector<char> opChar;
		prefixSoftClip = prefixHardClip = 0;
		suffixSoftClip = suffixHardClip = 0;
    if (clipping == hard) {
      SetHardClip(alignment, read, prefixHardClip, suffixHardClip);
		}
		
    if (clipping == soft or clipping == subread) {
      //
      // Even if clipping is soft, the hard clipping removes the 
      // low quality regions
      //
			if (clipping == soft) {
				prefixHardClip = read.lowQualityPrefix;
				suffixHardClip = read.lowQualitySuffix;
			}
			else if (clipping == subread) {
				prefixHardClip = max((DNALength) read.subreadStart, read.lowQualityPrefix);
				suffixHardClip = max((DNALength)(read.length - read.subreadEnd), read.lowQualitySuffix);
			}

			SetSoftClip(alignment, read, prefixHardClip, suffixHardClip, prefixSoftClip, suffixSoftClip);
		}

		//
		// The position of the alignment in the query and target.
		//
		int qPos = alignment.qPos + alignment.blocks[0].qPos;
		int tPos = alignment.tPos + alignment.blocks[0].tPos;

		if (prefixHardClip > 0) {
			opSize.push_back(prefixHardClip);
			opChar.push_back('H');
		}
		if (prefixSoftClip > 0) {
			opSize.push_back(prefixSoftClip);
			opChar.push_back('S');
		}
		
    CreateNoClippingCigarOps(alignment, qPos, tPos, opSize, opChar);
      
		//
		// Append the hard and soft clipping so they are in the order S
		// then H. 
		//
		if (suffixSoftClip > 0) {
			opSize.push_back(suffixSoftClip);
			opChar.push_back('S');
		}
		if (suffixHardClip > 0) {
			opSize.push_back(suffixHardClip);
			opChar.push_back('H');
		}

    if (alignment.tStrand == 1) {
      reverse(opSize.begin(), opSize.end());
      reverse(opChar.begin(), opChar.end());
    }

    CigarOpsToString(opSize, opChar, cigarString);
  }


  template<typename T_Alignment, typename T_Sequence>
  void PrintAlignment(T_Alignment &alignment,
                      T_Sequence &read,
                      ostream &samFile,
                      AlignmentContext &context,
											SupplementalQVList &qvlist,
                      Clipping clipping = none,
                      int subreadIndex = 0,
                      int nSubreads = 0) {
    
    string cigarString;
    uint16_t flag;
    T_Sequence alignedSequence;
    DNALength prefixSoftClip=0, suffixSoftClip=0;
    DNALength prefixHardClip=0, suffixHardClip=0;
		
    CreateCIGARString(alignment, read, cigarString, clipping, prefixSoftClip, suffixSoftClip, prefixHardClip, suffixHardClip);
    SetAlignedSequence(alignment, read, alignedSequence, clipping);
    BuildFlag(alignment, context, flag);
    samFile << alignment.qName << "\t" 
            << flag << "\t" 
            << alignment.tName << "\t";   // RNAME
		DNALength tAlignStart = alignment.TAlignStart() + alignment.blocks[0].tPos;
    if (alignment.tStrand == 0) {
      samFile << tAlignStart + 1 << "\t"; // POS, add 1 to
                                             // get 1 based
                                             // coordinate
                                             // system   
    }
    else {
      samFile << alignment.tLength - (tAlignStart + alignment.TEnd()) + 1 << "\t"; // includes - 1 for rev-comp,  +1 for one-based
    }
    samFile << (int) alignment.mapQV << "\t"// MAPQ
            << cigarString << "\t"; // CIGAR
      
      //
      // Determine RNEXT

    string rNext;
    rNext = "*";
    samFile << rNext << "\t"; // RNEXT
    
    DNALength nextSubreadPos = 0;
    samFile << nextSubreadPos << "\t"; // RNEXT, add 1 for 1 based
                                           // indexing


    DNALength tLen = alignment.GenomicTEnd() - alignment.GenomicTBegin();
    samFile << tLen << "\t"; // TLEN
    // Print the sequence on one line, and suppress printing the
    // newline (by setting the line length to alignedSequence.length
    ((DNASequence)alignedSequence).PrintSeq(samFile, 0);  // SEQ
    samFile << "\t";

    if (alignedSequence.qual.data != NULL) {
      alignedSequence.PrintAsciiQual(samFile, 0);  // QUAL
    }
    else {
      samFile <<"*";
    }
    samFile << "\t";
    //
    // Add optional fields
    //
    samFile << "RG:Z:" << context.readGroupId << "\t";
    samFile << "AS:i:" << alignment.score << "\t";

		// must recompute the 

    //
    // "RG" read group Id
    // "AS" alignment score
    // "XS" read alignment start position without counting previous soft clips (1 based) 
    // "XE" read alignment end position without counting previous soft clips (1 based) 
    // "XL" aligned read length 
    // "XQ" query sequence length
    // "XT" # of continues reads, always 1 for blasr 
    // "NM" edit distance 
    // "FI" read alignment start position (1 based) 
    //

		DNALength qAlignStart = alignment.QAlignStart();
		DNALength qAlignEnd   = alignment.QAlignEnd();

	  if (clipping == none) {
			samFile << "XS:i:" << qAlignStart + 1 << "\t";
			samFile << "XE:i:" << qAlignEnd + 1 << "\t";
    }
		else if (clipping == hard or clipping == soft or clipping == subread) {
			DNALength xs = prefixHardClip;
			DNALength xe = read.length - suffixHardClip;
			if (alignment.tStrand == 1) {
				xs = suffixHardClip;
				xe = read.length - prefixHardClip;
			}
			
			samFile << "XS:i:" << xs + 1 << "\t"; // add 1 for 1-based indexing in sam
			assert(read.length - suffixHardClip == prefixHardClip + alignedSequence.length);
			samFile << "XE:i:" << xe + 1 << "\t";
			samFile << "qs:i:" << read.qs << "\t"; // add 1 for 1-based indexing in sam
			assert(read.length - suffixHardClip == prefixHardClip + alignedSequence.length);
			samFile << "qe:i:" << read.qe << "\t";

		}
		samFile << "zm:i:" << read.holeNumber << "\t";
    samFile << "XL:i:" << alignment.qAlignedSeq.length << "\t";
    samFile << "XT:i:1\t"; // reads are allways continuous reads, not
                        // referenced based circular consensus when
                        // output by blasr.
    samFile << "NM:i:" << context.editDist << "\t";
    samFile << "FI:i:" << alignment.qAlignedSeqPos + 1;
    // Add query sequence length
    samFile << "\t" << "XQ:i:" << alignment.qLength << "\t";

		//
		// Write out optional quality values.  If qvlist does not 
		// have any qv's signaled to print, this is a no-op.
		//
		// First transform characters that are too large to printable ones.
		samFile.precision(4);
		samFile << "rq:f:" << read.readQuality << "\t";
		samFile << "np:i:" << read.numPasses << "\t";
		samFile << "cx:i:" << read.subreadContext << "\t";
		samFile << "sn:B:f,"
						<< read.snr[0] << ","
						<< read.snr[1] << ","
						<< read.snr[2] << ","
						<< read.snr[3];
		
		qvlist.FormatQVOptionalFields(alignedSequence);
		qvlist.PrintQVOptionalFields(alignedSequence, samFile);

    samFile << endl;

    alignedSequence.FreeIfControlled();
  }
  template<typename T_Sequence>
  void PrintUnalignedRead(T_Sequence &read,
                      ostream &samFile,
                      AlignmentContext &context,
                      SupplementalQVList &qvlist,
                      Clipping clipping = none,
                      int subreadIndex = 0,
                      int nSubreads = 0) {
    uint16_t flag = SEGMENT_UNMAPPED;
    T_Sequence alignedSequence;
    DNALength clippedReadLength = read.subreadEnd - read.subreadStart;
    DNALength clippedStartPos = read.subreadStart;

    // Return if subread sequence would be out of bounds.
    if (!(clippedStartPos >= 0 && clippedStartPos <= read.length && clippedReadLength >= 0 && clippedReadLength <= read.length)) {
      return;
    }

    alignedSequence.ReferenceSubstring(read, clippedStartPos, clippedReadLength);

    // Return if the subread has no sequence.
    if (alignedSequence.length == 0) {
      return;
    }

    // Unmapped reads have fixed defaults for most fields.
    samFile << read.GetName() << "\t" << flag << "\t*\t0\t0\t*\t*\t0\t0\t";   // RNAME, POS, MAP, CIGAR, RNEXT, PNEXT, TLEN

    ((DNASequence)alignedSequence).PrintSeq(samFile, 0);  // SEQ
    samFile << "\t";
		//
		// The quality is largely uninformative, so it is skipped and the null flag is output.
		//
		samFile <<"*";
    samFile << "\t";

    //
    // Add optional fields
    //
    samFile << "RG:Z:" << context.readGroupId << "\t";

		// must recompute the
    //
    // "RG" read group Id
    // "XQ" query sequence length
    // "XT" # of continues reads, always 1 for blasr
    //
	  if (clipping == subread) {
			DNALength xs = read.subreadStart;
			DNALength xe = read.subreadEnd;
			samFile << "XS:i:" << xs + 1 << "\t"; // add 1 for 1-based indexing in sam
			samFile << "XE:i:" << xe + 1 << "\t";
			samFile << "qs:i:" << read.qs << "\t"; // add 1 for 1-based indexing in sam
			samFile << "qe:i:" << read.qe << "\t";
		}
		else {
			samFile << "XS:i:1\t"
							<< "XE:i:"<<alignedSequence.length <<"\t"
							<< "qs:i:1\t"
							<< "qe:i:"<<alignedSequence.length << "\t";
		}
		samFile << "zm:i:" << read.holeNumber << "\t";
    samFile << "XL:i:" << 0 << "\t";
		samFile.precision(4);
		samFile << "rq:f:" << read.readQuality << "\t";
		samFile << "np:i:" << read.numPasses << "\t";
		samFile << "cx:i:" << read.subreadContext << "\t";
		samFile << "sn:B:f,"
						<< read.snr[0] << ","
						<< read.snr[1] << ","
						<< read.snr[2] << ","
						<< read.snr[3] << "\t";
    samFile << "XT:i:1"; // reads are allways continuous reads, not
                        // referenced based circular consensus when
                        // output by blasr.
    // Add query sequence length
    samFile << "\t" << "XQ:i:" << alignedSequence.length;

		//
		// Write out optional quality values.  If qvlist does not
		// have any qv's signaled to print, this is a no-op.
		//
		// First transform characters that are too large to printable ones.
		qvlist.FormatQVOptionalFields(read);
		qvlist.PrintQVOptionalFields(read, samFile);

    samFile << endl;

    alignedSequence.FreeIfControlled();
  }
};




#endif
