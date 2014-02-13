#ifndef ALIGNMENT_BLOCK_H_
#define ALIGNMENT_BLOCK_H_

#include <iostream>
#include <fstream>
#include "DNASequence.h"

using namespace std;
class Block {
 public:
	//
	// An alignment is a collection of blocks. The qPos and tPos in a block
	// is relative to the beginning of the alignment rather than the
	// target or query.
	//
	 
	DNALength qPos, tPos, length;
	friend ostream &operator<<(ostream &out, const Block &b) {
		out << " q: " << b.qPos << " t: " << b.tPos << " len: " << b.length;
		return out;
	}

	Block& Assign(Block &rhs) {
		qPos = rhs.qPos;
		tPos = rhs.tPos;
		length = rhs.length;
		return *this;
	}

	DNALength QEnd() {
		return qPos + length;
	}

	DNALength TEnd() {
		return tPos + length;
	}

  void Clear() {
    qPos = tPos =  length = 0;
  }
};



#endif
