#ifndef DATASTRUCTURES_ALIGNMENT_BYTE_ALIGNMENT_H_
#define DATASTRUCTURES_ALIGNMENT_BYTE_ALIGNMENT_H_
#include <vector>
#include "algorithms/alignment/AlignmentUtils.h"
#include "DNASequence.h"

using namespace std;
/*
 * These arrays are for going from the HDF byte alignment format to characters, or
 * from characters to the upper or lower nybble for translating from character 
 * alignments to byte alignments.
 */

typedef vector<unsigned char> ByteAlignment;
static char QueryChar[256] = {
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', //0	 
' ', ' ', ' ', ' ', ' ', ' ', 'A', 'A', 'A', 'A', //10	 
'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', // 20 
'A', 'A', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', // 30 
'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', ' ', ' ', // 40 
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', // 50 
' ', ' ', ' ', ' ', 'G', 'G', 'G', 'G', 'G', 'G', // 60 
'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', 'G', // 70 
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', // 80 
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', // 90 
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', // 100
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', // 110
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'T', 'T', // 120
'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'T', // 130
'T', 'T', 'T', 'T', ' ', ' ', ' ', ' ', ' ', ' ', // 140
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', // 150
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', // 160
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', // 170
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', // 180
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', // 190
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', // 200
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', // 210
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', // 220
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', // 230
' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', // 240
' ', ' ', ' ', ' ', ' ', ' ' 										// 250
};

static char RefChar[256] = {
' ', 'A', 'C', ' ', 'G', ' ', ' ', ' ', 'T', ' ', 
' ', ' ', ' ', ' ', ' ', ' ', ' ', 'A', 'C', ' ', 
'G', ' ', ' ', ' ', 'T', ' ', ' ', ' ', ' ', ' ', 
' ', ' ', ' ', 'A', 'C', ' ', 'G', ' ', ' ', ' ', 
'T', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'A', 
'C', ' ', 'G', ' ', ' ', ' ', 'T', ' ', ' ', ' ', 
' ', ' ', ' ', ' ', ' ', 'A', 'C', ' ', 'G', ' ', 
' ', ' ', 'T', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 
' ', 'A', 'C', ' ', 'G', ' ', ' ', ' ', 'T', ' ', 
' ', ' ', ' ', ' ', ' ', ' ', ' ', 'A', 'C', ' ', 
'G', ' ', ' ', ' ', 'T', ' ', ' ', ' ', ' ', ' ', 
' ', ' ', ' ', 'A', 'C', ' ', 'G', ' ', ' ', ' ', 
'T', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'A', 
'C', ' ', 'G', ' ', ' ', ' ', 'T', ' ', ' ', ' ', 
' ', ' ', ' ', ' ', ' ', 'A', 'C', ' ', 'G', ' ', 
' ', ' ', 'T', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 
' ', 'A', 'C', ' ', 'G', ' ', ' ', ' ', 'T', ' ', 
' ', ' ', ' ', ' ', ' ', ' ', ' ', 'A', 'C', ' ', 
'G', ' ', ' ', ' ', 'T', ' ', ' ', ' ', ' ', ' ', 
' ', ' ', ' ', 'A', 'C', ' ', 'G', ' ', ' ', ' ', 
'T', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'A', 
'C', ' ', 'G', ' ', ' ', ' ', 'T', ' ', ' ', ' ', 
' ', ' ', ' ', ' ', ' ', 'A', 'C', ' ', 'G', ' ', 
' ', ' ', 'T', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 
' ', 'A', 'C', ' ', 'G', ' ', ' ', ' ', 'T', ' ', 
' ', ' ', ' ', ' ', ' ', ' '};



static unsigned char RefAlignmentByte[256] = {0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  1,  0,  2,  0,  0, 
																							0,  4,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  8,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  1,  0,  2, 
																							0,  0,  0,  4,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  8,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																							0,  0,  0,  0,  0,  0};
	

static unsigned char QueryAlignmentByte[256] = { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  16,  0,  32,  0,  0, 
																								 0,  64,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  128,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  16,  0,  32, 
																								 0,  0,  0,  64,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  128,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
																								 0,  0,  0,  0,  0,  0};


void MakeReverseComplementByteAlignment(const unsigned char *byteAlignment,
                                        UInt length,
                                        unsigned char *byteAlignmentRC) {
  unsigned char q,t;
  int i;
  for (i = 0; i < length; i++) {
    if (QueryChar[byteAlignment[i]] == ' ') { q = ' '; }
    else { q = ReverseComplementNuc[QueryChar[byteAlignment[i]]]; }

    if (RefChar[byteAlignment[i]] == ' ') { t = ' '; }
    else { t = ReverseComplementNuc[RefChar[byteAlignment[i]]]; }
    
    byteAlignmentRC[length - i - 1] = QueryAlignmentByte[q] + RefAlignmentByte[t];
  }
}


void ByteAlignmentToQueryString(const unsigned char* byteAlignment,
																UInt length,
																char* charAlignment) {
	int i;
	for (i = 0; i < length; i++) {
		charAlignment[i] = QueryChar[byteAlignment[i]];
	}
}


void ByteAlignmentToRefString(const unsigned char* byteAlignment,
															UInt length,
															char* charAlignment) {
	int i;
	for (i = 0; i < length; i++) {
		charAlignment[i] = RefChar[byteAlignment[i]];
	}
}

void RemoveGaps(string &gappedStr, string &ungappedStr) {
	ungappedStr = gappedStr;
	int i, i2;
	i = i2 = 0;
	for (i = 0;  i < ungappedStr.size(); i++ ){
		if (ungappedStr[i] != ' ') {
			ungappedStr[i2] = ungappedStr[i];
			i2++;
		}
	}
	ungappedStr.resize(i2);
}


void GappedStringsToAlignment(string &gappedQuery, string &gappedRef, Alignment &alignment) {
int qPos = 0, rPos = 0;
	int i = 0; // position in alignment string
	while (i < gappedQuery.size()) {
		while (i < gappedQuery.size() and (gappedQuery[i] == ' ' or gappedRef[i] == ' ')) {
			if (gappedQuery[i] != ' ') {
				qPos++;
			}
			if (gappedRef[i] != ' ') {
				rPos++;
			}
			i++;
		}
		int queryBlockStart, queryBlockEnd, refBlockStart, refBlockEndl;
		Block b;
		b.qPos = qPos;
		b.tPos = rPos;
		while (i < gappedQuery.size() and gappedQuery[i] != ' ' and gappedRef[i] != ' ') {
			i++;
			qPos++;
			rPos++;
		}
		b.length = qPos - b.qPos;
		alignment.blocks.push_back(b);
	}
}


void ByteAlignmentToAlignment(vector<unsigned char> &byteAlignment, Alignment &alignment) {
	string readSequence, refSequence;
	readSequence.resize(byteAlignment.size());
	refSequence.resize(byteAlignment.size());
	
	ByteAlignmentToQueryString(&byteAlignment[0], byteAlignment.size(), &readSequence[0]);
	ByteAlignmentToRefString(&byteAlignment[0], byteAlignment.size(), &refSequence[0]);				
	GappedStringsToAlignment(readSequence, refSequence, alignment);
}

void AlignmentToByteAlignment(Alignment &alignment, 
                              DNASequence &querySeq, DNASequence &refSeq,
                              vector<unsigned char> &byteAlignment) {
  string refStr, alignStr, queryStr;
  CreateAlignmentStrings(alignment, querySeq, refSeq, 
                         refStr, alignStr, queryStr);
  byteAlignment.resize(refStr.size());
  int i;
  for (i = 0; i < refStr.size(); i++) {
    byteAlignment[i] = RefAlignmentByte[refStr[i]] + QueryAlignmentByte[queryStr[i]];
  }
}

bool IsMatch(vector<unsigned char> &byteAlignment, int i) {
	if (QueryChar[byteAlignment[i]] != ' ' and 
			RefChar[byteAlignment[i]] != ' ' and 
			(QueryChar[byteAlignment[i]] == RefChar[byteAlignment[i]])) {
		return true;
	}
	else {
		return false;
	}
}

void CountStats(vector<unsigned char> &byteAlignment, 
								int &nMatch, int &nMismatch, int &nIns, int &nDel, 
								int start=0, int end=-1) {
	int i;
	if (end == -1) {
		end = byteAlignment.size();
	}

	nMatch = nMismatch = nIns = nDel = 0;

	for (i = start; i < end; i++) {
		if (QueryChar[byteAlignment[i]] == ' ') {
			nDel++;
		}
		else if (RefChar[byteAlignment[i]] == ' ') {
			nIns++;
		}
		else if (RefChar[byteAlignment[i]] != QueryChar[byteAlignment[i]]) {
			nMismatch++;
		}
		else {
			nMatch++;
		}
	}
}

int CountBasesInReference(vector<unsigned char> &byteAlignment) {
	int i;
	int nBases = 0;
	for (i = 0; i < byteAlignment.size(); i++) {
		if (RefChar[byteAlignment[i]] != ' ') {
			nBases++;
		}
	}
	return nBases;
}

int CountBasesInQuery(vector<unsigned char> &byteAlignment) {
	int i;
	int nBases = 0;
	for (i = 0; i < byteAlignment.size(); i++) {
		if (QueryChar[byteAlignment[i]] != ' ') {
			nBases++;
		}
	}
	return nBases;
}

int CountNMatches(vector<unsigned char> &byteAlignment) {
	int nMatches = 0;
	int i;
	for (i = 0; i < byteAlignment.size(); i++) {
		if (IsMatch(byteAlignment, i)) {
			nMatches++;
		}
	}
	return nMatches;
}

float ComputePacBioAccuracy(vector<unsigned char> &byteAlignment) {
	int m, mm, i, d;
	CountStats(byteAlignment, m, mm, d, i);
	int readLength = CountBasesInQuery(byteAlignment);
	return 1 - (1.0*mm + d + i)/readLength;
}


float ComputePercentIdentity(vector<unsigned char> &byteAlignment) {
	int i;
	int nMatch = CountNMatches(byteAlignment);
	return (1.0*nMatch) / byteAlignment.size();
}

void CreateSequenceToAlignmentMap(vector<unsigned char> &byteAlignment, 
                                  vector<int> &baseToAlignmentMap) {
  int alignPos, ungappedAlignPos;
  int alignmentLength = byteAlignment.size();
  baseToAlignmentMap.resize(alignmentLength);
  for (ungappedAlignPos = 0, alignPos = 0; 
       alignPos < alignmentLength; 
       alignPos++) {

    if (QueryChar[byteAlignment[alignPos]] != ' ') {
      baseToAlignmentMap[ungappedAlignPos] = alignPos;
      ++ungappedAlignPos;
    }
  }
  baseToAlignmentMap.resize(ungappedAlignPos);
}			

void CreateAlignmentToSequenceMap(vector<unsigned char> &byteAlignment, 
                                  vector<int> &alignmentToBaseMap) {
  int alignPos, ungappedAlignPos;
  int alignmentLength = byteAlignment.size();
  alignmentToBaseMap.resize(alignmentLength);
  for (ungappedAlignPos = 0, alignPos = 0; 
       alignPos < alignmentLength; 
       alignPos++) {
    if (QueryChar[byteAlignment[alignPos]] != ' ') {
      alignmentToBaseMap[alignPos] = ungappedAlignPos;
      ++ungappedAlignPos;
    }
  }
}			

#endif
