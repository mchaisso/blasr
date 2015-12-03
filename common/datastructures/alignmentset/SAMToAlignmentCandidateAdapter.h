#ifndef SAM_TO_ALIGNMENT_CANDIDATE_ADAPTER_H_
#define SAM_TO_ALIGNMENT_CANDIDATE_ADAPTER_H_

#include "SAMAlignment.h"
#include "SAMFlag.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include <map>

using namespace std;

void InitializeCandidateFromSAM(SAMAlignment &sam,
                                AlignmentCandidate<> &candidate) {
  candidate.qName = sam.qName;
  candidate.tName = sam.rName;
  candidate.tPos  = sam.pos;
  candidate.mapQV = sam.mapQV;
}

//
// Change a list of gap operations into a gap structure between
// blocks. 
//
int ProcessGap(vector<int> &lengths, vector<char> &ops, 
               int &opIndex, int opEnd, GapList &gaps, 
               int &qAdvance, int &tAdvance) {
  //
  // Default is no gap.
  //
  qAdvance = tAdvance = 0;
  gaps.clear();
  if (opIndex >= lengths.size()) {
    return opIndex;
  }

  //
  // Look to see if the suffix of operations stars with a gap.
  //
  if (ops[opIndex] != 'I' and
      ops[opIndex] != 'D') {
    return opIndex;
  }

  while (opIndex < opEnd and (ops[opIndex] == 'I' or ops[opIndex] == 'D')) {
    if (ops[opIndex] == 'I') {
      gaps.push_back(Gap(Gap::Target, lengths[opIndex]));
      qAdvance += lengths[opIndex];
    }
    else {
      assert(ops[opIndex] == 'D');
      gaps.push_back(Gap(Gap::Query, lengths[opIndex]));
      tAdvance += lengths[opIndex];
    }
    opIndex++;
  }
  return opIndex;
}

bool IsClipping(char c) {
  return (c == 'S' or c == 'H');
}

int AdvancePastClipping(vector<int> &lengths, vector<char> &ops, int &opIndex, int &numSoftClipped) {
  int numClipped = 0;
  numSoftClipped = 0;
  while (opIndex < lengths.size() and IsClipping(ops[opIndex])) {
    numClipped += lengths[opIndex];
    if (ops[opIndex] == 'S') {
      numSoftClipped+= lengths[opIndex];
    }

    ++opIndex;
  }
  return numClipped;
}

int IsSkipped(char c) {
  return c == 'N';
}

int AdvancePastSkipped(vector<int> &lengths, vector<char> &ops, int &opIndex) {
  int numSkipped = 0;
  while (opIndex < lengths.size() and IsSkipped(ops[opIndex])) {
    numSkipped += lengths[opIndex];
    opIndex++;
  }
  return numSkipped;
}


//
// Returns true if the corresponding characters are assigned to each
// other in the alignment.  This is either a match or a mismatch.
//

bool IsAssignChar(char c) {
  return (c == 'M' or c == 'X' or c == '=');
}

int ProcessMatch(vector<int> &lengths, vector<char> &ops,
                  int &opIndex, int opEnd) {

  //
  // Make sure this starts on some sort of blockn
  //

  int blockLength = 0;
  while (opIndex < opEnd and (IsAssignChar(ops[opIndex]))) {
    blockLength += lengths[opIndex];
    opIndex++;
  }
  return blockLength;
}
                 

void CIGAROpsToBlocks(vector<int> &lengths, vector<char> &ops,
                      int cigarStart,
                      int cigarEnd,
                      AlignmentCandidate<> &aln) {
  cigarStart = 0;
  cigarEnd   = lengths.size();
  CIGAROpsToBlocks(lengths, ops, cigarStart, cigarEnd, aln);
}
  
int AdvancePosToAlignmentEnd(vector<char> &ops, int &pos) {
  int start = pos;
  while (pos < ops.size() and ops[pos] != 'N' and !IsClipping(ops[pos])) {
    pos++;
  }
  return pos - start;
}

int GetAlignedQueryLengthByCIGARSum(vector<char> &ops, vector<int> &lengths) {
  int i;
  for (i = 0; i < ops.size(); i++) {
    if (ops[i] != 'S' && ops[i] != 'H') {
      break;
    }
  }
  int queryLength = 0;
  for (; i < ops.size() && ops[i] != 'S' && ops[i] != 'H'; i++) {
    if (IsAssignChar(ops[i]) or ops[i] == 'I' or ops[i] == 'N') {
      queryLength += lengths[i];
    }
  }
  return queryLength;
}





int GetAlignedReferenceLengthByCIGARSum(vector<char> &ops, vector<int> &lengths) {
  int i;
  for (i = 0; i < ops.size(); i++) {
    if (ops[i] != 'S' && ops[i] != 'H') {
      break;
    }
  }
  int refLength = 0;
  for (; i < ops.size() && ops[i] != 'S' && ops[i] != 'H'; i++) {
    if (IsAssignChar(ops[i]) or ops[i] == 'D' or ops[i] == 'N') {
      refLength += lengths[i];
    }
  }
  return refLength;
}
      

void CIGAROpsToBlocks(vector<int> &lengths, vector<char> &ops,
                      int &cigarPos,
                      int &cigarEnd,
                      int &qPos, int &tPos,
                      AlignmentCandidate<> &aln) {

  int gapIndex = 0;
  DNALength qStart = qPos, tStart = tPos;
  assert(cigarPos >= cigarEnd or !IsClipping(ops[cigarPos]));

  //
  // Advance past any skipped portion.
  //
  int numSkipped = AdvancePastSkipped(lengths, ops, cigarPos);
  tPos += numSkipped;

  //
  // Process the gaps before the first match.
  //

  //
  // If there is nothing left, just bail.
  GapList gap;

  cigarEnd = cigarPos;
  AdvancePosToAlignmentEnd(ops, cigarEnd);
  if (cigarPos >= cigarEnd) {
    return;
  }
  
  // 
  // Process any gap that the aligner produces before the first match
  // begins.
  //
  int qAdvance, tAdvance;
  ProcessGap(lengths, ops, cigarPos, cigarEnd, gap, qAdvance, tAdvance);
  aln.gaps.push_back(gap);
  qPos += qAdvance;
  tPos += tAdvance;
  //
  // Now add gaps.
  //
  while(cigarPos < cigarEnd) {
    //
    // The next operation must be a match.
    //
    int matchLength = ProcessMatch(lengths, ops, cigarPos, cigarEnd);
    Block b;
    b.qPos = qPos - qStart;
    b.tPos = tPos - tStart;
    b.length = matchLength;
    aln.blocks.push_back(b);
    qPos += b.length;
    tPos += b.length;

    ProcessGap(lengths, ops, cigarPos, cigarEnd, gap, qAdvance, tAdvance);
    aln.gaps.push_back(gap);
    tPos += tAdvance;
    qPos += qAdvance;
  }
  
}

void ReverseAlignmentOperations(vector<int> &lengths, vector<char> &ops) {
  reverse(lengths.begin(), lengths.end());
  reverse(ops.begin(), ops.end());
}

void SAMAlignmentsToCandidates(SAMAlignment &sam,
                               vector<FASTASequence> &referenceSequences,
                               map<string,int> &refIDToListIndex,
                               vector<AlignmentCandidate<> > &candidates, 
                               bool parseSmrtTitle = false,
                               bool keepRefAsForward = true) {
  //
  // First determine how many alignments there are from CIGAR string.
  //
  vector<int> lengths;
  vector<char> ops;
  sam.cigar.Vectorize(lengths, ops);

  DNASequence querySeq;
  // For now just reference the query sequence.
  querySeq.deleteOnExit = false;
  querySeq.seq = (Nucleotide*) sam.seq.c_str();
  querySeq.length = sam.seq.size();

  DNALength samTEnd = 0;
  DNALength samTStart = sam.pos - 1;
  if (keepRefAsForward == false and IsReverseComplement(sam.flag)) {
    ReverseAlignmentOperations(lengths, ops);
    DNASequence rcQuerySeq;
    querySeq.CopyAsRC(rcQuerySeq);
    //
    // Zero out the query seq so that the string memory is not
    // deleted.
    //
    querySeq.seq = NULL;
    querySeq.length = 0;
    querySeq = rcQuerySeq;
    rcQuerySeq.Free();
    samTEnd = GetAlignedReferenceLengthByCIGARSum(ops, lengths);
  }


  int i;
  int offset = 0;
  if (ops.size() == 0) {
    return;
  }
  bool alignmentStarted = false;
  bool onFirstMatch = true;
  int  curAlignment;
  
  //
  // Advance past any clipping.  This advances in both query and
  // reference position.
  //
  int cigarPos = 0;
  int qPos = 0; 
  int tPos = 0;

  DNALength queryPosOffset = 0;
  if (parseSmrtTitle) {
    //
    // The aligned sequence is really a subread of a full
    // sequence. The position of the aligments start at 0, the
    // beginning of the query sequence, but in the sam file, they
    // may appear as subreads, and are offset from the start of the
    // subread.  By convention, the subread coordinates are embedded
    // in the title of the query, if it is a smrtTitle. 
    // Two types of smrtTitle are supported:
    // movie/zmw/start_end
    // movie/zmw/start_end/start2_end2
    vector<string> values;
    ParseSeparatedList(sam.qName, values, '/');
    DNALength qStart = 0, qEnd = 0;
    bool worked = false;

    if (values.size() >= 3) {
      vector<string> offsets;
      ParseSeparatedList(values[2], offsets, '_');
      if (offsets.size() == 2) {
        qStart = atoi(offsets[0].c_str());
        qEnd   = atoi(offsets[1].c_str());
        if (values.size() == 3) {
          worked = true;
        } else if (values.size() == 4) {
          offsets.clear();
          ParseSeparatedList(values[3], offsets, '_');
          if (offsets.size() == 2) {
            qEnd   = qStart + atoi(offsets[1].c_str());
            qStart = qStart + atoi(offsets[0].c_str());
            worked = true;
          }
        }
      }
    }
    if (worked == false) {
      cout<<values.size() << endl;
      cout << "ERROR. Could not parse title " << sam.qName << endl;
      exit(1);
    }
    queryPosOffset = qStart;
  }
  else if (sam.xs) {
    queryPosOffset += sam.xs - 1;
  }


  while (cigarPos < lengths.size()) {
    int numClipped;
    //
    // Sequence clipping becomes offsets into the q/t alignedSeqPos
    //


    int numSoftClipped;
    numClipped = AdvancePastClipping(lengths, ops, cigarPos, numSoftClipped);

    //
    // End loop now.
    //
    if (cigarPos >= lengths.size()) {
      break;
    }
    qPos += numSoftClipped;

    //
    // Skipped sequences are just advances in the tPos.
    //
    int numSkipped = AdvancePastSkipped(lengths, ops, cigarPos);
    tPos += numSkipped;

    if (cigarPos >= lengths.size()) {
      break;
    }


    AlignmentCandidate<> alignment;
    //
    // The aligned sequence must start at a match therefore the tpos
    // and qpos are 0.
    //
    alignment.qPos = 0;
    alignment.tPos = 0;

    DNALength qAlignStart = qPos;  // qAlignStart is the start of the alignment relative to the sequence in the SAM file.
    DNALength tAlignStart = tPos;  // tAlignStart is the start of the alignment in the genome.
    
    int cigarEnd = cigarPos;
    AdvancePosToAlignmentEnd(ops, cigarEnd);

    CIGAROpsToBlocks(lengths, ops,          
                     cigarPos, cigarEnd,
                     qPos, tPos,
                     alignment);


    DNALength queryLengthSum = GetAlignedQueryLengthByCIGARSum(ops, lengths);
    DNALength refLengthSum   = GetAlignedReferenceLengthByCIGARSum(ops, lengths);
    alignment.qAlignedSeqLength = qPos - qAlignStart;
    alignment.tAlignedSeqLength = tPos - tAlignStart;

    //
    // Assign candidate sequences.
    //
    // First, the query sequence is straight from the SAM line.
		if (qAlignStart + alignment.qAlignedSeqLength > querySeq.length) {
			cerr << "Could not create query sequence." << endl;
			candidates.clear();
			return;
		}

    ((DNASequence*)&alignment.qAlignedSeq)->Copy(querySeq, qAlignStart, alignment.qAlignedSeqLength);
    
    // The SAM Alignments a
    alignment.qStrand = IsReverseComplement(sam.flag);
    alignment.tStrand = 0;
    alignment.mapQV   = sam.mapQV;

    //
    // Assign the offsets into the original sequence where the
    // subsequence starts.
    //

    alignment.qAlignedSeqPos = queryPosOffset + qAlignStart;    
    alignment.tAlignedSeqPos = samTStart + tAlignStart;
    
    if (sam.rName == "*") {
      //
      // No reference, do not add the alignment to the list of
      // candidates.
      //
      continue;
    }
    else {
      int refIndex;
      int s = refIDToListIndex.size();
      if (refIDToListIndex.find(sam.rName) == refIDToListIndex.end()) {
        cout <<" ERROR.  SAM Reference " << sam.rName << " is not found in the list of reference contigs." << endl;
        exit(1);
      }
      
      refIndex = refIDToListIndex[sam.rName];
     
      alignment.tLength = referenceSequences[refIndex].length;
      alignment.qLength = sam.seq.size(); 
      alignment.qName = sam.qName;
      alignment.tName = sam.rName;


      if (keepRefAsForward == false and alignment.qStrand == 1) {

        //
        // Now that the reference sequence has been copied, if it is
        // on the reverse strand, make the reverse complement for
        // proper printing.
        //
        alignment.tAlignedSeqPos = samTStart + (samTEnd - tAlignStart - alignment.tAlignedSeqLength);
				if (alignment.tAlignedSeqLength > referenceSequences[refIndex].length ||
						alignment.tAlignedSeqPos    > referenceSequences[refIndex].length ||
						alignment.tAlignedSeqLength + alignment.tAlignedSeqPos > referenceSequences[refIndex].length + 2) {
					//alignment.tAlignedSeqPos is 1 based and unsigned.
					cout << "WARNING. The mapping of read " << alignment.qName  
							 << " to reference "      << alignment.tName 
							 << " is out of bounds."  << endl
							 << "         StartPos (" << alignment.tAlignedSeqPos  
							 << ") + AlnLength (" << alignment.tAlignedSeqLength 
							 << ") > RefLength (" << referenceSequences[refIndex].length
							 << ") + 2 "          << endl;
					continue;
				}
        ((DNASequence*)&alignment.tAlignedSeq)->Copy(referenceSequences[refIndex], alignment.tAlignedSeqPos, alignment.tAlignedSeqLength);             
        alignment.tAlignedSeq.ReverseComplementSelf();
        // either ref or read is defined as being in the forward
        // orientation.  Here, since refAsForward is false, the read
        // is forward.  Since the read is forward, the aligned
        // sequences are stored as the reverse complement of the read
        // and the references.
        //
        alignment.tStrand = 1;
        alignment.qStrand = 0;
      }
      else {
        if (alignment.tAlignedSeqLength > referenceSequences[refIndex].length ||
						alignment.tAlignedSeqPos    > referenceSequences[refIndex].length ||
						alignment.tAlignedSeqLength + alignment.tAlignedSeqPos > referenceSequences[refIndex].length + 2) {
					//alignment.tAlignedSeqPos is 1 based and unsigned. 
					cout << "WARNING. The mapping of read " << alignment.qName  
							 << " to reference "      << alignment.tName 
							 << " is out of bounds."  << endl
							 << "         StartPos (" << alignment.tAlignedSeqPos  
							 << ") + AlnLength (" << alignment.tAlignedSeqLength 
							 << ") > RefLength (" << referenceSequences[refIndex].length
							 << ") + 2 "          << endl;
					continue;
				}
        ((DNASequence*)&alignment.tAlignedSeq)->Copy(referenceSequences[refIndex], 
                                                     alignment.tAlignedSeqPos, 
                                                     alignment.tAlignedSeqLength);
      }
    }

    if (alignment.blocks.size() > 0) {
      candidates.push_back(alignment);
    }
  }
  if (candidates.size() > 0 and keepRefAsForward == false and candidates[0].tStrand == 1) {
    reverse(candidates.begin(), candidates.end());
  }
  querySeq.Free();
}


#endif
