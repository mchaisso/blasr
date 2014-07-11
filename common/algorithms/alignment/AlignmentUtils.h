#ifndef ALIGNMENT_UTILS_H_
#define ALIGNMENT_UTILS_H_

#include <vector>
#include "defs.h"
#include "NucConversion.h"
#include "utils.h"
#include "datastructures/matrix/FlatMatrix.h"
#include "datastructures/alignment/Alignment.h"

using namespace std;


enum AlignmentType { Local,     // Standard Smith-Waterman
										 Global,    // Standard Needleman-Wuncsh
										 QueryFit,  // No gap penalties from query.  
										 TargetFit, // For when the query engulfs the target.  
										 Overlap,    // No gap penalty at the beginning of
																 // query, nor at the end of text **
																 // not yet implemented for k-band **, 
										 FrontAnchored, // Require the alignment to align
																		// pos 0,0 in the matrix 
										 EndAnchored,   // Require the alignment to align
																		// the pos m,n in the matrix 
										 // Score alignment types solely compute the score
										 // of an alignment and do not store the actual
										 // alignment itself.  This is fast for filtering
										 // potential alignments that may be re-aligned
										 // later to store the actual alignment.
										 Fit,  // No gap penalties at the beginning nor ends of alignments.
										 TSuffixQPrefix, // Same as overlap 
										 TPrefixQSuffix, // so that the order of the alignment does not have to be reversed
										 ScoreGlobal,
										 ScoreLocal,
										 ScoreQueryFit,
										 ScoreTargetFit,
										 ScoreOverlap,
										 ScoreFrontAnchored,
										 ScoreEndAnchored,
										 ScoreTSuffixQPrefix,
										 ScoreTPrefixQSuffix,
										 //
										 // A LocalBoundaries alignment is in-between a
										 // score-only and full-alignment. The full path
										 // matrix is computed, but rather 
										 // than computing an alignment, simply the  
										 // (qStart, qLength), (tStart, tLength)
										 // coordinates of the alignment are returned.  
										 //
										 LocalBoundaries,
										 //
										 // A SignificanceLimited alignment is a banded
										 // alignment that continues alignment until the
										 // score drops below a certain threshold below
										 // the maximum score.
										 //
										 SignificanceLimited
};

template<typename T_ScoreFn>
int ComputeAlignmentScore(string &queryStr, string &textStr, T_ScoreFn &scoreFn, 
                          bool useAffineScore = false) { 
	if (queryStr.size() != textStr.size()) {
		cout << "Computing alignment score using invalid alignment string." << endl;
		cout << "Bailing out."<<endl;
		exit(1);
	}
	VectorIndex i;
	int score = 0;
  int alignStrLen = queryStr.size();
	for(i = 0; i < alignStrLen; i++) {
		if (queryStr[i] != '-' and
				textStr[i] != '-') {
			score += scoreFn.scoreMatrix[ThreeBit[(int)queryStr[i]]][ThreeBit[(int)textStr[i]]];
		}
		else {
      if (useAffineScore) {
        //
        // Compute affine gap scoring.  For now this uses symmetric insertion/deletion penalties. 
        //
        int gapEnd = i;
        while (gapEnd < queryStr.size() and gapEnd < textStr.size() and 
               (queryStr[gapEnd] == '-' or textStr[gapEnd] == '-')) {
          ++gapEnd;
        }
        int gapLength = gapEnd - i;
				int affineScore = gapLength * scoreFn.affineExtend + scoreFn.affineOpen;
				int gapScore;
				if (queryStr[gapEnd-1] == '-') {
					gapScore = scoreFn.del * gapLength;
				}
				else {
					gapScore = scoreFn.ins * gapLength;
				}
				if (gapScore < affineScore) {
					score += gapScore;
				}
				else {
					score += affineScore;
				}
        //
        // Advance past gap -1, so that at the top of the for loop i
        // will be at the end of the gap.
        //
        i = gapEnd - 1;
      }
      else {
        //
        // Use non-affine gap scoring.
        //
        if (queryStr[i] == '-' and textStr[i] != '-') {
          score += scoreFn.del;
        }
        else if (queryStr[i] != '-' and textStr[i] == '-') {
          score += scoreFn.ins;
        }
        else {
          score += scoreFn.scoreMatrix[4][4];
        }
      }
		}
	}
	return score;
}
	

template<typename T_QuerySequence, typename T_TargetSequence, typename T_ScoreFn>
  int ComputeAlignmentScore(Alignment &alignment,
                            T_QuerySequence &query,
                            T_TargetSequence &text,
                            T_ScoreFn &scoreFn, 
                            bool useAffinePenalty = false) {
  VectorIndex b, q, t, l, bi, gi;
  int alignmentScore = 0;

	for (b = 0; b < alignment.blocks.size(); b++ ) {
		for ((q = alignment.qPos + alignment.blocks[b].qPos,
					t = alignment.tPos + alignment.blocks[b].tPos,
					l = 0) ; 
				 l < alignment.blocks[b].length; q++, t++, l++) {
			alignmentScore += scoreFn.Match(text, t, query, q);
		}
    if (alignment.gaps.size() == alignment.blocks.size() + 1) {
      for (gi = 0; gi < alignment.gaps[b+1].size(); gi++) {
        if (alignment.gaps[b+1][gi].seq == Gap::Target) {
          if (useAffinePenalty and 
							(scoreFn.affineOpen + alignment.gaps[b+1][gi].length * scoreFn.affineExtend < 
							 alignment.gaps[b+1][gi].length * scoreFn.ins ) ) {
            alignmentScore += scoreFn.affineOpen + alignment.gaps[b+1][gi].length * scoreFn.affineExtend;
          }
          else {
            alignmentScore += alignment.gaps[b+1][gi].length * scoreFn.ins;
          }
        }
        else {
          if (useAffinePenalty and (scoreFn.affineOpen + alignment.gaps[b+1][gi].length * scoreFn.affineExtend < 
																		alignment.gaps[b+1][gi].length * scoreFn.del)) {
            alignmentScore += scoreFn.affineOpen + alignment.gaps[b+1][gi].length * scoreFn.affineExtend;
          }
          else {
            alignmentScore += alignment.gaps[b+1][gi].length * scoreFn.del;
          }
        }
      }
    }
  }
  return alignmentScore;

}
                            
/*
template<typename T_QuerySequence, typename T_TargetSequence, typename T_ScoreFn>
int ComputeAlignmentScore(Alignment &alignment,
													T_QuerySequence &query,
													T_TargetSequence &text,
                          T_ScoreFn &scoreFn) {
	VectorIndex b;
	int alignmentScore = 0;
	VectorIndex q, t;
	DNALength l;

	//
	// Compute the average match score for the case that
	// there are unassigned alignments.
	//
	int totalMismatchPenalty = 0, averageMismatchPenalty;
	for (q = 0; q < 4; q++ ) {
		for (t = 0; t < 4; t++ ){ 
			if (q != t) totalMismatchPenalty += scoreFn.scoreMatrix[q][t];
		}
	}

	averageMismatchPenalty = (int) (totalMismatchPenalty/ 12.0);

	for (b = 0; b < alignment.blocks.size(); b++ ) {
		for ((q = alignment.qPos + alignment.blocks[b].qPos,
					t = alignment.tPos + alignment.blocks[b].tPos,
					l = 0) ; 
				 l < alignment.blocks[b].length; q++, t++, l++) {
			alignmentScore += scoreFn.matchScores[ThreeBit[query.seq[q]]][ThreeBit[text.seq[t]]];
		}
		if (alignment.blocks.size() > 0 and b < alignment.blocks.size() - 1) {
			// 
			// Add the gap penalties.
			//
			int queryGap = (alignment.blocks[b+1].qPos - 
											(alignment.blocks[b].qPos +
											 alignment.blocks[b].length));
			int textGap  = (alignment.blocks[b+1].tPos - 
											(alignment.blocks[b].tPos + 
											 alignment.blocks[b].length));

			int commonGap = queryGap;
			if (textGap < queryGap) {
				commonGap = textGap;
			}
			
			queryGap -= commonGap;
			textGap  -= commonGap;
			alignmentScore += averageMismatchPenalty * commonGap;
			
			alignmentScore += (queryGap * scoreFn.ins);
			alignmentScore += (textGap * scoreFn.del);
		}
	}
	return alignmentScore;
}
											 
*/

inline int GetNumberWidth(unsigned int value) {
	// 0 has a width of 1.
	int width = 1;
	while (value / 10 > 0)  {
		value = value / 10;
		++width;
	}
	return width;
}

template<typename T_Alignment>
inline void PrintCompareSequencesAlignmentStats(T_Alignment &alignment, ostream &out) {
  int lastBlock;

  lastBlock = alignment.blocks.size() -1;
	
  int qLength, tLength;
  if (lastBlock >= 0) {
    qLength = (alignment.blocks[lastBlock].qPos 
               + alignment.blocks[lastBlock].length) ;
    tLength =  (alignment.blocks[lastBlock].tPos 
                + alignment.blocks[lastBlock].length);
  }
  else {
    qLength = tLength = 0;
  }

  //
  // First print the query 
  //
  int alignmentQStart;
  if (lastBlock >= 0) {
    alignmentQStart = alignment.qPos + alignment.qAlignedSeqPos;
  }
  else {
    alignmentQStart = alignment.qAlignedSeqPos;
  }
  out << alignment.qName
      << " " << alignment.qLength
      << " " << alignmentQStart
      << " " << alignmentQStart + qLength;

  if (alignment.qStrand == 0) {
    out << " +";
  }
  else {
    out << " -";
  }
	
  int alignmentTStart;
  if (lastBlock >= 0) {
    alignmentTStart = alignment.tPos + alignment.tAlignedSeqPos;
  }
  else {
    alignmentTStart = 0;
  }
  out << " " << alignment.tName 
      << " " << alignment.tLength;
  if (alignment.tStrand == 0) {
    out << " " << alignmentTStart
        << " " << alignmentTStart + tLength;
		out << " + ";
  }
  else {
    out	<< " " << alignment.tLength - (alignmentTStart + tLength)
        << " " << alignment.tLength - (alignmentTStart);
		out << " - ";
  }
		
	out << alignment.score 
			<< " " << alignment.nMatch 
			<< " " << alignment.nMismatch
			<< " " << alignment.nIns
			<< " " << alignment.nDel
			<< " " << (int) alignment.mapQV 
			<< " ";
}

template<typename T_Alignment>
inline int ReadCompareSequencesAlignmentStats(istream &in, T_Alignment &alignment) {
	int qEnd, tEnd;
	int qLength, tLength;
	char qStrand, tStrand;
	if (!(in >> alignment.qName )) return 0;
	if (!(in >> qLength)) return 0;
	if (!(in >> alignment.qPos)) return 0;
	if (!(in >> qEnd)) return 0;
	if (!(in >> qStrand)) return 0;
	if (!(in >> alignment.tName)) return 0;
	if (!(in >> tLength)) return 0;
	if (!(in >> alignment.tPos)) return 0;
	if (!(in >> tEnd)) return 0;
	if (!(in >> tStrand)) return 0;
	if (!(in >> alignment.score)) return 0;
	if (!(in >> alignment.nMatch)) return 0;
	if (!(in >> alignment.nMismatch)) return 0;
	if (!(in >> alignment.nIns)) return 0;
	if (!(in >> alignment.nDel)) return 0;
	return 1;
}


template<typename T_Alignment>
inline int ReadCompSeqAlignment(istream &in, T_Alignment &alignment) {
	if (!ReadCompareSequencesAlignmentStats(in, alignment)) return 0;
	string alignStr;
	if (!(in >> alignment.qString)) return 0;
	if (!(in >> alignStr)) return 0;
	if (!(in >> alignment.tString)) return 0;
	string eol;
	getline(in, eol);
	return 1;
}

/*
 * This should be changed to read any type of alignment, since templates are being
 * used.
 */
inline void ReadCompSeqAlignments(string &compSeqAlignmentFileName, vector<CompSeqAlignment> &alignments) {
	ifstream in;
	CrucialOpen(compSeqAlignmentFileName, in);
	CompSeqAlignment alignment;
	while (ReadCompSeqAlignment(in, alignment)) {
		alignments.push_back(alignment);
	}
}

inline void PrintAlignmentStats(Alignment &alignment, ostream &out) {
	out << "    nMatch: " << alignment.nMatch << endl;
	out << " nMisMatch: " << alignment.nMismatch << endl;
	out << "      nIns: " << alignment.nIns << endl;
	out << "      nDel: " << alignment.nDel << endl;
	out << "      %sim: " << alignment.pctSimilarity << endl;
	out << "     Score: " << alignment.score << endl; 
}

template<typename T_QuerySequence, typename T_TargetSequence>
void  AppendGapCharacters(Gap &gap, 
                          T_QuerySequence &query, T_TargetSequence &text, 
                          DNALength &q, DNALength &t,
                          char mismatchChar, char gapChar,
                      string &textStr, string &alignStr, string &queryStr) {
  int gp;
  for (gp = 0; gp < gap.length; gp++) {
    if (gap.seq == Gap::Query) {
      textStr.push_back(text[t]);
      alignStr.push_back(mismatchChar);
      queryStr.push_back(gapChar);
      t++;
    }
    else if (gap.seq == Gap::Target) {
      textStr.push_back(gapChar);
      alignStr.push_back(mismatchChar);
      queryStr.push_back(query[q]);
      q++;
    }
  }
}

template<typename T_Alignment, typename T_QuerySequence, typename T_TargetSequence>
void CreateAlignmentStrings(T_Alignment &alignment, 
														T_QuerySequence &query, T_TargetSequence &text, 
														string &textStr, string &alignStr, string &queryStr, DNALength queryLength=0, DNALength textLength=0) {
	DNALength q = alignment.qPos;
	DNALength t = alignment.tPos;
	DNALength qPos, tPos;
	DNALength  g;
	char mismatchChar = '*';
	char matchChar = '|';
	char gapChar = '-';
  char gapSeparationChar = ' ';
	if (alignment.blocks.size() == 0) {
		textStr = "";
		alignStr = "";
		queryStr = "";
		return;
	}

  if (alignment.gaps.size() == 0) {
    //
    // If there is no gap list, add the gaps as an offset here.
    //
    if (alignment.blocks[0].qPos > 0 or
        alignment.blocks[0].tPos > 0) {
      // commonGapLen should be the shorter gap.
      qPos = alignment.blocks[0].qPos;
      tPos = alignment.blocks[0].tPos;
      int commonGapLen = qPos;
      if (commonGapLen > tPos) {
        commonGapLen = tPos;
      }
      for (g = 0; g < commonGapLen; g++ ) {
        textStr.push_back(text[t]);
        alignStr.push_back(mismatchChar);
        queryStr.push_back(query[q]);
        t++;
        q++;
      }
      tPos -= commonGapLen;
      qPos -= commonGapLen;

      //
      // one of tPos or qPos is now 0.
      // The other represents extra sequence
      // that should be output before starting the alignment
      DNALength p;
      for (p = 0; p < tPos; p++) {
        textStr.push_back(text[t]);
        alignStr.push_back(gapSeparationChar);
        queryStr.push_back(gapChar);
        t++;
      }
      for (p = 0; p < qPos; p++) {
        textStr.push_back(gapChar);
        alignStr.push_back(gapSeparationChar);
        queryStr.push_back(query[q]);
        q++;
      }
    }
  }

  //
  // Add gap characters if they are before the beginning of the alignment. 
  // This shouldn't happen, but for some local alignments, it can.
  //
	DNALength b, bl, gi;
  if (alignment.gaps.size() > 0) {
    // The first gap is before the first block of characters.
    int gi;
    for (gi = 0; gi < alignment.gaps[0].size(); gi++) {
      AppendGapCharacters(alignment.gaps[0][gi], query, text, q, t, mismatchChar, gapChar, textStr, alignStr, queryStr);
    }
  }

	for (b = 0; b < alignment.size() ; b++) {
		for (bl = 0; bl < alignment.blocks[b].length; bl++ ) {
			queryStr.push_back(query[q]);
			textStr.push_back(text[t]);
			assert(queryLength == 0 or q < queryLength);
			assert(textLength == 0 or t < textLength);
			if (TwoBit[query[q]] != 
					TwoBit[text[t]])
				alignStr.push_back(mismatchChar);
			else
				alignStr.push_back(matchChar);
			q++;
			t++;
		}
		//
		//  There are no gaps to count after the last block, so 
		//  don't add the gapped characters for this.
		//
		if (alignment.blocks.size() == 0)
			continue;
		if (b == alignment.blocks.size() - 1) {
			continue;
		}
    if (alignment.gaps.size() > 0) {
      for (gi = 0; gi < alignment.gaps[b+1].size(); gi++) {
        int gp;
        AppendGapCharacters(alignment.gaps[b+1][gi], query, text, q, t, gapSeparationChar, gapChar, textStr, alignStr, queryStr);
      }
    }
    else {


      int queryGapLen = (alignment.blocks[b+1].qPos - 
                         alignment.blocks[b].qPos - alignment.blocks[b].length);
      int textGapLen  = (alignment.blocks[b+1].tPos - 
                         alignment.blocks[b].tPos - alignment.blocks[b].length);

      if (queryGapLen > 0 or textGapLen > 0) {
        // commonGapLen should be the shorter gap.
        int commonGapLen = queryGapLen; 
        if (queryGapLen > textGapLen) {
          commonGapLen = textGapLen;
        }
        textGapLen -= commonGapLen;
        queryGapLen -= commonGapLen;

        for (g = 0; g < queryGapLen; g++, q++ ){
          textStr.push_back(gapChar);
          alignStr.push_back(gapSeparationChar);
          queryStr.push_back(query[q]);
        }
        for (g = 0; g < textGapLen; g++, t++ ){
          textStr.push_back(text[t]);
          alignStr.push_back(gapSeparationChar);
          queryStr.push_back(gapChar);

        }

        for (g = 0; g < commonGapLen; g++ ) {
          textStr.push_back(text[t]);
          alignStr.push_back(gapSeparationChar);
          queryStr.push_back(query[q]);
          t++;
          q++;
        }
      }
    }
	}
}

template<typename T_Alignment, typename T_ScoreFn>
  void ComputeAlignmentStats(T_Alignment &alignment, Nucleotide* qSeq, Nucleotide *tSeq, T_ScoreFn &scoreFn, bool useAffineScore = false) {
  //int matchMatrix[5][5], int ins, int del, bool useAffineScore = false) {
	int qp = 0, tp = 0;
	int nMatch = 0, nMismatch = 0, nIns =0, nDel = 0;
	float pctSimilarity;
	string textStr, alignStr, queryStr;
	CreateAlignmentStrings(alignment, qSeq, tSeq, textStr, alignStr, queryStr);
	int i;
	int alignLength = textStr.size();

	for (i = 0; i < alignLength; i++ ) {
		if ((textStr[i] != '-') and (queryStr[i] != '-')) {
			if (ThreeBit[textStr[i]] == ThreeBit[queryStr[i]]) {
				nMatch++;
			}
			else {
				nMismatch++;
			}
			tp++;
			qp++;
		}
		else if (textStr[i] == '-' and queryStr[i] != '-') {
			nIns++;
			qp++;
		}
		else if (queryStr[i] == '-' and textStr[i] != '-') {
			nDel++;
			tp++;
		}
	}
	if (tp + qp > 0) {
		if (textStr.size() + queryStr.size() > 0) {
			pctSimilarity = (nMatch*2.0) / (textStr.size() + queryStr.size()) * 100;
		}
		else {
			pctSimilarity = 0;
		}
	}
	else {
		pctSimilarity = 0;
	}

	alignment.score = ComputeAlignmentScore(queryStr, textStr, scoreFn, useAffineScore);
	alignment.nMatch = nMatch;
	alignment.nMismatch = nMismatch;
	alignment.nDel = nDel;
	alignment.nIns = nIns;
	alignment.pctSimilarity = pctSimilarity;
}

int ComputeDrift(Block &cur,  Block &next ) {

	int tGap  = (next.tPos - cur.TEnd());
	int qGap = (next.qPos - cur.QEnd());

	int commonGap = 0;

	if (tGap > 0 and qGap > 0) {
		commonGap = abs(tGap - qGap);
	}
	tGap -= commonGap;
	qGap -= commonGap;

	return tGap - qGap;
}

template<typename T_Alignment>
int ComputeDrift(T_Alignment &alignment) {
	VectorIndex b;
	int qGap = 0, tGap = 0, commonGap = 0;
	int drift = 0;
	int maxDrift = 0;
	int driftBetweenBlocks;
	if (alignment.blocks.size() == 0)
		return 0;
	for (b = 0; b < alignment.blocks.size() - 1; b++) { 
		driftBetweenBlocks = ComputeDrift(alignment.blocks[b], alignment.blocks[b+1]);
		drift += driftBetweenBlocks;
		if (abs(drift) > maxDrift) 
			maxDrift = abs(drift);
	}
	return maxDrift;
}

template<typename T_Alignment>
void RemoveAlignmentPrefixGaps(T_Alignment &alignment) {
  if (alignment.gaps.size()  == 0) {
    return;
  }
  int g;
  int tStart = 0, qStart = 0;
  for (g = 0; g < alignment.gaps[0].size(); g++) {
    if (alignment.gaps[0][g].seq == Gap::Target) {
      qStart += alignment.gaps[0][g].length;
    }
    else if (alignment.gaps[0][g].seq == Gap::Query) {
      tStart += alignment.gaps[0][g].length;
    }
  }
  int b;
  int nBlocks;
  for (b = 0, nBlocks = alignment.blocks.size(); b < nBlocks; b++) {
    alignment.blocks[b].qPos -= qStart;
    alignment.blocks[b].tPos -= tStart;
  }
  alignment.gaps[0].clear();
  alignment.tPos += tStart;
  alignment.qPos += qStart;
}
#endif
