#ifndef FASTA_READER_H_
#define FASTA_READER_H_


#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
#include <string.h>
#include <limits.h>
#include <sys/types.h>
#include <unistd.h>
#include <ext/rope>
#include "sys/mman.h"
#include "sys/fcntl.h"
#include "datastructures/metagenome/SequenceIndexDatabase.h"
#include "FASTASequence.h"

using namespace std;

class FASTAReader {
 protected:
	ifstream file;

	int padding;
	char endOfReadDelim;
	char readStartDelim;
  bool doToUpper;

  unsigned char *convMat;
	//
	// Quick check to see how much to read.
	//

 public:
  bool computeMD5;
  string curReadMD5;
	bool storeName;
  void Init() {
		padding = 0;
		endOfReadDelim = '>';
		readStartDelim = '>';
    doToUpper = false;
    convMat = PreserveCase;
    computeMD5 = false;
		storeName = false;
  }
  FASTAReader() {
    Init();
	};

 
FASTAReader(string &fileName) {
  Init(); // initialze defaults.
  Init(fileName); // open file.
  padding = 0;
  endOfReadDelim = '>';
  readStartDelim = '>';
  }

	void SetSpacePadding(int _padding) {
		assert(_padding >= 0);
		padding = _padding;
	}
		
  void SetToUpper() {
    doToUpper = true;
    convMat   = AllToUpper;
  }

	//
	// Synonym for Init() for consistency.
	//
	int Initialize(string &seqInName) {
		return Init(seqInName);
	}

  int Init(string &seqInName, int passive=0) {
		// close file handle just in case
		file.close();
		file.clear();
		file.open(seqInName.c_str());

		if (file.good() == false) {
				cout << "Could not open FASTA file " << seqInName << endl;
				exit(1);
		}
		return 1;
	}
	
  bool AdvanceToTitleStart(char delim='>') {
		while (!file.eof() && file.peek() != delim) {
			char c;
			c = file.get();
		}
		if (file.eof()) {
			return false;
		}
		else {
			return true;
		}
	}


	void CheckValidTitleStart(char delim='>') {
		char c = file.peek();
		if (file.eof() == true or c != delim) {
			cout << "ERROR, sequence entry must begin with \"" << delim << "\"" << endl;
			exit(1);
		}
	}

	long ReadAllSequencesIntoOne(FASTASequence &seq, SequenceIndexDatabase<FASTASequence> *seqDBPtr=NULL) {
		vector<FASTASequence> sequences;
		ReadAllSequences(sequences);
		int i;
		long seqLength=0;
		for (i=0;i<sequences.size();i++) {
			seqLength += sequences[i].length;
		}
		
		if (seqLength > UINT_MAX) {
			cout << "ERROR! Reading fasta files greater than 4Gbytes is not supported." << endl;
			exit(1);
		}
		seq.seq = new Nucleotide[seqLength+sequences.size()];
		unsigned int seqPos = 0;
		for (i=0;i<sequences.size();i++) {
			memcpy(&seq.seq[seqPos], sequences[i].seq, sequences[i].length);
			seq.seq[seqPos+sequences[i].length] = 'N';
			string md5str;
			seqPos+= sequences[i].length+1;
			if (seqDBPtr != NULL) {
				seqDBPtr->growableName.push_back(sequences[i].title);
				MakeMD5((const char*) sequences[i].seq, sequences[i].length, md5str);
				seqDBPtr->md5.push_back(md5str);
				seqDBPtr->growableSeqStartPos.push_back(seqPos);
			}
			sequences[i].Free();
		}
		seq.length = seqPos;
		if (seqDBPtr != NULL) {
			seqDBPtr->Finalize();
		}
		return seq.length;
	}

	void ReadTitle(char *&title, int &titleLength) {
		// 
		// Extract the title.  The length of the title does not include the newline.
		//
		string line;
		getline(file, line);
		title = new char[line.size()];
		if (line.size() > 1) {
			memcpy(title, &line.c_str()[1], line.size()-1);
			title[line.size()-1] = '\0';
			titleLength = line.size()-1;
		}
		else {
			title = NULL;
			titleLength = 0;
		}
	}
	
  int GetNext(FASTASequence &seq, char delim='>') {
		if (file.eof()) {
			return 0;
		}
		

		int foundStart = AdvanceToTitleStart(delim);
		if (foundStart == false) {
			return false;
		}
  
		// 
		// Make sure there is a '>'
		//
		CheckValidTitleStart(delim);
		
		ReadTitle(seq.title, seq.titleLength);

		if (storeName) {
			string name = seq.GetName();
			seq.CopyTitle(name);
		}
		//
		// Read in the next sequence.
		//

		// Count the length of the sequence.
		__gnu_cxx::crope inputRead;
		int nLines =0;
		
		while(file.good() and file.eof() == false and file.peek() != endOfReadDelim) {
			string line;
			getline(file, line);
			stringstream lineStrm(line);
			while (lineStrm) {
				string word;
				lineStrm >> word;
				if (word.size() > 0) {
					inputRead.append(word.c_str());
					nLines+=1;
				}
			}
		}
		seq.length = inputRead.size();
		seq.seq = new Nucleotide[seq.length];
		inputRead.copy((char*) seq.seq);

    if (computeMD5) {
      MakeMD5((const char*) seq.seq, seq.length, curReadMD5);
    }
		return 1;
	}

   int Advance(int nSeq) {

		int nAdvanced = 0;
		// base case -- it's always ok to advance 0
		if (nSeq == 0) { return 1; }
		// Make sure this starts on a sequence.
		AdvanceToTitleStart();
		while (file.eof() == false and nSeq > 0) {
			string line;
			getline(file, line);
			if (file.peek() == endOfReadDelim) {
				nSeq-=1;
			}
		}
		return true;
	}

	int CriticalGetNext(FASTASequence &seq) {
		if (!GetNext(seq)) {
			cout << "Could not read a sequence." << endl;
			exit(1);
		}
	}
	int ConcatenateNext(FASTASequence &cur) {
		FASTASequence next;
		int retVal;	
		if ((retVal = GetNext(next))) {
			next.CleanupASCII();
			cur.Concatenate((Nucleotide*) "N");
			cur.Concatenate(next);	
			delete[] next.seq;
		}
		return retVal;
	}

	void Close() {
		file.close();
		file.clear();
	}
	
	void ReadAllSequences(vector<FASTASequence> &sequences) {
		//
		// Step 1, compute the number of reads in the file.
		// 
		
		long p;
		p = 0;
		int nSeq = 0;
		FASTASequence seq;
		while (GetNext(seq)) {
			sequences.push_back(seq);
		}
	}

}; 


#endif
