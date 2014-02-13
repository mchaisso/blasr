#ifndef CMPSEQ_REVERSE_COMPRESS_INDEX_H_
#define CMPSEQ_REVERSE_COMPRESS_INDEX_H_
#include <iostream>
#include <fstream>

class ReverseCompressIndex {
 public:
	int *index;
	int indexLength;
	int binSize;
	int maxRun;
	int size() { return indexLength;}
  
  ReverseCompressIndex() {
    index = NULL;
    indexLength = binSize = maxRun = 0;
  }

	void Write(ofstream &out) {
		out.write((char*) &indexLength, sizeof(int));
		out.write((char*) &binSize, sizeof(int));
		out.write((char*) &maxRun, sizeof(int));
		out.write((char*) index, sizeof(int) * indexLength);
	}

  void Read(ifstream &in) {
		in.read((char*) &indexLength, sizeof(int));
		in.read((char*) &binSize, sizeof(int));
		in.read((char*) &maxRun, sizeof(int));
		index = new int[indexLength];
		in.read((char*) index, sizeof(int) *indexLength);
	} 
	
	void ShallowCopy(ReverseCompressIndex &rhs) {
		index = rhs.index;
		indexLength = rhs.indexLength;
		binSize = rhs.binSize;
		maxRun  = rhs.maxRun;
	}
};


#endif
