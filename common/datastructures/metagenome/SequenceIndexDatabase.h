#ifndef SEQUENCE_INDEX_DATABASE_H_
#define SEQUENCE_INDEX_DATABASE_H_

#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>
#include <stdlib.h>
#include <sstream>
#include <algorithm>
#include "Types.h"
#include "DNASequence.h"
#include "utils/StringUtils.h"
#include <map>
using namespace std;

#define SEQUENCE_INDEX_DATABASE_MAGIC 1233211233
template<typename TSeq>
class SequenceIndexDatabase {
 public:
	vector<DNALength> growableSeqStartPos;
	vector<string> growableName;
  
	DNALength *seqStartPos;
	bool deleteSeqStartPos;
	char **names;
	bool	deleteNames;
	int  *nameLengths;
	bool deleteNameLengths;
	int nSeqPos;
	bool deleteStructures;
	map<string, DNALength> startPos;
	map<string, DNALength> endPos;
  //
  // This is stored after reading in the sequence.
  //
  vector<string> md5;

	SequenceIndexDatabase(int final=0) {
    nSeqPos = 0;
		if (!final) {
			growableSeqStartPos.push_back(0);
		}
		names = NULL; deleteNames = false;
		nameLengths = NULL; deleteNameLengths = false;
		seqStartPos = NULL; deleteSeqStartPos = false;
		deleteStructures = false;
	}

	void BuildNameMaps() {
	  int i;
	  for (i = 0; i < nSeqPos-1; i++) {
		startPos[names[i]] = seqStartPos[i];
		endPos[names[i]] = seqStartPos[i+1]-1;
	  }
	}
	
	void GetBoundaries(string name, DNALength &queryStartPos, DNALength &queryEndPos) {
	  if (startPos.find(name) != startPos.end()) {
		queryStartPos = startPos[name];
		queryEndPos = endPos[name];
	  }
	  else {
		queryStartPos = 0;
		queryEndPos = 0;
	  }
	}

  DNALength GetLengthOfSeq(int seqIndex) {
    assert(seqIndex < nSeqPos-1);
    return seqStartPos[seqIndex+1] - seqStartPos[seqIndex] - 1;
  }

  // Return index of a reference sequence with name "seqName".
  int GetIndexOfSeqName(string seqName) {
    for(int i = 0; i < nSeqPos - 1; i++) {
      if (seqName == string(names[i])) {
        return i;
      }
    }
    return -1;
  }
  
  void GetName(int seqIndex, string &name) {
    assert(seqIndex < nSeqPos-1);
    name = names[seqIndex];
  }

  void MakeSAMSQString(string &sqString) {
    stringstream st;
    int i;
    for (i = 0; i < nSeqPos-1; i++) {
      st << "@SQ\tSN:" << names[i] << "\tLN:" << GetLengthOfSeq(i);
      if (md5.size() == nSeqPos-1) {
        st << "\tM5:" << md5[i];
      }
      st << endl;
    }
    sqString = st.str();
  }
	
  DNALength ChromosomePositionToGenome(int chrom, DNALength chromPos) {
    assert(chrom < nSeqPos);
    return seqStartPos[chrom] + chromPos;

  }

	int SearchForIndex(DNALength pos) {
		// The default behavior for the case
		// that there is just one genome.
		if (nSeqPos == 1) {
			return 0;
		}
    
    DNALength* seqPosIt = upper_bound(seqStartPos+1, seqStartPos + nSeqPos, pos);

    return seqPosIt - seqStartPos - 1;
	}

	string GetSpaceDelimitedName(unsigned int index) {
		int pos;
		assert(index < nSeqPos);
		string name;
		for (pos = 0; pos < nameLengths[index]; pos++) {
			if (names[index][pos] == ' ' or 
					names[index][pos] == '\t' or 
          names[index][pos] == '\0') {
				break;
			}
		}
		name.assign(names[index], pos);
		return name;
	}
		
	int SearchForStartBoundary(DNALength pos) {

		int index = SearchForIndex(pos);
		if (index != -1) {
			return seqStartPos[index];
		}
		else {
			return -1;
		}
	}

	int SearchForEndBoundary(DNALength pos) {
		
		int index = SearchForIndex(pos);
		if (index != -1) {
			return seqStartPos[index + 1];
		}
		else {
			return -1;
		}
	}

	DNALength SearchForStartAndEnd(DNALength pos, DNALength &start, DNALength &end) {
		int index = SearchForIndex(pos);
		if (index != -1) {
			start = seqStartPos[index];
			end   = seqStartPos[index+1];
			return 1;
		}
		else {
			start = end = -1;
			return 0;
		}
	}

	void WriteDatabase(ofstream &out) {
	  int mn = SEQUENCE_INDEX_DATABASE_MAGIC;
		out.write((char*) &mn, sizeof(int));
		out.write((char*) &nSeqPos, sizeof(int));
		out.write((char*) seqStartPos, sizeof(DNALength) * nSeqPos);
		int nSeq = nSeqPos - 1;
		out.write((char*) nameLengths, sizeof(int) * nSeq);
		int i;
		//
		// The number of sequences is 1 less than the number of positions
		// since the positions include 0 as a boundary.
		//
		char nullchar = '\0';
		for (i = 0; i < nSeq; i++) { 
			//
			// nameLengths has space for the null char, so the length of the
			// name = nameLengths[i]-1. Write a nullchar to disk so that it
			// may be read in later with no work.
			// 
			out.write((char*) names[i], sizeof(char) * (nameLengths[i]-1));
      out.write((char*) &nullchar, sizeof(char));
		}
	}

	void ReadDatabase(ifstream &in) {
		int mn;
		// Make sure this is a read database, since the binary input
		// is not syntax checked.
		in.read((char*) &mn, sizeof(int));
		if (mn != SEQUENCE_INDEX_DATABASE_MAGIC) {
			cout << "ERROR: Sequence index database is corrupt!" << endl;
			exit(1);
		}

		//
		// Read in the boundaries of each sequence.
		//
		deleteStructures = true;

		in.read((char*) &nSeqPos, sizeof(int));
		seqStartPos = new DNALength[nSeqPos];
		deleteSeqStartPos = true;
		in.read((char*) seqStartPos, sizeof(DNALength) * nSeqPos);
		int nSeq = nSeqPos - 1;

		// Get the lengths of the strings to read.
		nameLengths = new int[nSeq];
		deleteNameLengths = true;
		in.read((char*)nameLengths, sizeof(int) * nSeq);

		// Get the titles of the sequences.
		names = new char*[nSeq];
		deleteNames = true;
		char *namePtr;
		int i;
		for (i = 0; i < nSeq; i++) { 
			namePtr = new char[nameLengths[i]];
			if (nameLengths[i] > 0) {
				in.read(namePtr, nameLengths[i]);
			}
			namePtr[nameLengths[i]-1] = '\0';
			names[i] = namePtr;
		}
	}

	void SequenceTitleLinesToNames() {
		int seqIndex;
        vector<string> tmpNameArray;
		for (seqIndex = 0; seqIndex < nSeqPos-1; seqIndex++) {
			string tmpName;
			AssignUntilFirstSpace(names[seqIndex], nameLengths[seqIndex], tmpName);
			delete[] names[seqIndex];
			names[seqIndex] = new char[tmpName.size()+1];
			strcpy(names[seqIndex], tmpName.c_str());
			names[seqIndex][tmpName.size()] = '\0';
			nameLengths[seqIndex] = tmpName.size();

            tmpNameArray.push_back(tmpName);
		}
        // Make sure that reference names are unique.
        sort(tmpNameArray.begin(), tmpNameArray.end());
        for(int j = 0; j < tmpNameArray.size() - 1; j++) {
            if (tmpNameArray[j] == tmpNameArray[j+1]) {
                cout << "Error, reference with name \"" 
                     << tmpNameArray[j] 
                     << "\" in the reference genome is not unique"<<endl;
                exit(1);
            }
        }
	}

	VectorIndex AddSequence(TSeq &sequence) {
		int endPos = growableSeqStartPos[growableSeqStartPos.size() - 1];
	  int growableSize = growableSeqStartPos.size();
		growableSeqStartPos.push_back(endPos + sequence.length + 1);
		string fastaTitle;
		sequence.GetFASTATitle(fastaTitle);
		growableName.push_back(fastaTitle);
		return growableName.size();
	}

	void Finalize() {
		deleteStructures  = true;
		seqStartPos = &growableSeqStartPos[0];
		nSeqPos = growableSeqStartPos.size();
		int nSeq = nSeqPos - 1;
		names = new char*[nSeq];
		deleteNames = true;
		unsigned int i;
		nameLengths = new int[nSeq];
		deleteNameLengths = true;
		for (i = 0; i < nSeq; i++) {
			names[i] = new char[growableName[i].size() + 1];
			memcpy((char*) names[i], (char*) growableName[i].c_str(), growableName[i].size());
			names[i][growableName[i].size()] = '\0';
			nameLengths[i] = growableName[i].size() + 1;
		}
	}
	

	void FreeDatabase() {
		int i;
		if (deleteStructures == false) {
			return;
		}
		if (names != NULL and deleteNames) {
			int nSeq = nSeqPos - 1;
			for (i = 0; i < nSeq; i++ ){
				delete[] names[i];
			}
			delete[] names;
		}
		if (nameLengths != NULL) {
			delete[] nameLengths;
		}
		if (seqStartPos != NULL and deleteSeqStartPos) {
			delete[] seqStartPos;
		}
	}
};


template< typename TSeq >
class SeqBoundaryFtr {
 public:
	SequenceIndexDatabase<TSeq> *seqDB;

	SeqBoundaryFtr(SequenceIndexDatabase<TSeq> *_seqDB) {
		seqDB = _seqDB;
	}

  int GetIndex(DNALength pos) {
    return seqDB->SearchForIndex(pos);
  }

  int GetStartPos(int index) {
    assert(index < seqDB->nSeqPos);
    return seqDB->seqStartPos[index];
  }

  void GetBoundaries(const char * name, DNALength &queryStart, DNALength &queryEnd) {
	seqDB->GetBoundaries(name, queryStart, queryEnd);
  }

	DNALength operator()(DNALength pos) {
		return seqDB->SearchForStartBoundary(pos);
	}

  //
  // This is misuse of a functor, but easier interface coding for now.
  DNALength Length(DNALength pos) {
    DNALength start, end;
    seqDB->SearchForStartAndEnd(pos, start, end);
    return end - start - 1;
  }
};

#endif
