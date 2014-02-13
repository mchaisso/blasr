#ifndef  DNA_SEQUENCE_H_
#define  DNA_SEQUENCE_H_
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include "NucConversion.h"
#include "Types.h"
using namespace std;

typedef uint32_t DNALength;
typedef unsigned char Nucleotide;

class DNASequence {
 public:
  DNALength length;
  Nucleotide *seq;
	bool deleteOnExit;
	DNALength size() {
		return length;
	}

  void TakeOwnership(DNASequence &rhs) {
    if (deleteOnExit) {
      if (seq != NULL) {
        delete[] seq;
      }
    }
    seq = rhs.seq;
    length = rhs.length;
    deleteOnExit = rhs.deleteOnExit;
  }

  
  void Append(const DNASequence &rhs, DNALength appendPos=0) {
    //
    // Simply append rhs to this seuqence, unless appendPos is nonzero
    // in which case rhs is inserted at attendPos, overwriting this
    // sequence from appendPos to the end.
    //
    Nucleotide *newSeq;
    //
    // Handle the two cases (appendPos == 0 and appendPos > 0)
    // separately in order to handle memory deallocation correctly.
    //
    if (appendPos == 0) {
      DNALength  newSeqLength = length + rhs.length;
      newSeq = new Nucleotide[newSeqLength];
      memcpy(newSeq, seq, length);
      memcpy(&newSeq[length], rhs.seq, rhs.length);

      if (length != 0) {
        delete[] seq;
      }
      seq = newSeq;
      length = newSeqLength;
      deleteOnExit = true;
    }
    else {
      if (appendPos + rhs.length < length) {
        memcpy(&seq[appendPos], rhs.seq, rhs.length);
        length = appendPos + rhs.length;
      }
      else {
        DNALength lengthCopy = length;
        length = appendPos;
        DNALength newSeqLength;
        newSeqLength = length + rhs.length;
        newSeq = new Nucleotide[newSeqLength];
        memcpy(newSeq, seq, length);
        memcpy(&newSeq[length], rhs.seq, rhs.length);
        if (deleteOnExit and lengthCopy != 0) {
          delete[] seq;
        }
        seq = newSeq;
        length = newSeqLength;
        deleteOnExit = true;
      }
    }
  }

  // Copie FROM rhs to this DNASequence. 
	typedef Nucleotide T_Block;
	DNASequence &Copy(const DNASequence &rhs, DNALength rhsPos=0, DNALength rhsLength=0) {
		if (length != 0) {
			if (seq != NULL)
			    delete[] seq;
            seq = NULL;
            length = 0;
		}

    //
    // When initializing a vector of DNASequence's, the copy
    // constructor will initialze a list and call this
    // function with a zero-length DNASequence as the rhs to
    // initialize every element in the vector   The check
    // below will fail on zero-length sequences, so add a boundary
    // condition check before that to allow the copy-constructor to
    // work.
    //
    if (rhs.length == 0) {
      return *this;
    }
      
    //
    // Silently ignoring this case could lead to problems later on,
    // catastrophically assert here if the input is not valid.
	  // In case rhsLength + rhsPos > ULONG_MAX (4294967295), check 
    // both rhsLength and rhsPos, fix bug 21794
	  //

    if (not (rhsLength <= rhs.length     && 
             rhsPos    <= rhs.length + 1 &&
             rhsLength + rhsPos <= rhs.length + 2 )) {
      cout << "ERROR.  The subsequence to copy is out of bounds." << endl
           << "        Failed to copy a subsequence starting at " << rhsPos << endl
           << "        with length "<< rhsLength 
           << " from a sequence of length " << rhs.length << "." << endl;
      exit(1);
	  }
    
    if (rhsLength == 0) {
      rhsLength = rhs.length - rhsPos;
    }
    if (rhsLength == 0) {
      seq = NULL;
    }
    else {
      seq = new Nucleotide [rhsLength];
      memcpy(seq, &rhs.seq[rhsPos], rhsLength);
    }
		length = rhsLength;
		deleteOnExit = true;
		return *this;
	}

	void ShallowCopy(const DNASequence &rhs) {
		seq = rhs.seq;
		length = rhs.length;
		deleteOnExit = false;
	}
	int GetStorageSize() {
		return (length * sizeof(Nucleotide));
	}

	DNASequence &operator=(const DNASequence &rhs){ 
		Copy(rhs);
		return *this;
	}

	int bitsPerNuc;
	DNASequence() {
		seq = NULL;
		length = 0;
		bitsPerNuc = 8;
		deleteOnExit = false;
	}
	//
	// synonym for printseq
	//
	void Print(ostream &out, int lineLength = 50) {
		PrintSeq(out, lineLength);
	}

	void PrintSeq(ostream &out, int lineLength = 50) {
    if (lineLength == 0) {
			string line;
			line.assign((char*)seq, length);
			out << line;
    }
    else {
      //
      // Make sure this isn't 
      assert(lineLength > 0);
      DNALength curPos = 0;
      int curLineLength = lineLength;
      while (curPos < length) {
        if (curPos + curLineLength > length) {
          curLineLength = length - curPos;
        }
        string line;
        line.assign((char*) &seq[curPos], curLineLength);
        out << line << endl;
        curPos += curLineLength;
      }
    }
	}

	void Allocate(DNALength plength) {
		if (seq != NULL) {
			delete[] seq;
		}
		seq = new Nucleotide [plength];
		length = plength;
		deleteOnExit = true;
	}

	void ReferenceSubstring(const DNASequence &rhs, UInt pos=0, int substrLength=0) {
		//
		// This makes a reference therefore it should not be deleted.
		//
        assert(pos >= 0 && pos <= rhs.length &&
               substrLength >= 0 && substrLength <= rhs.length);
		if (substrLength == 0) {
			substrLength = rhs.length - pos;
		}
		assert(pos + substrLength <= rhs.length);
		seq = &rhs.seq[pos];
		length = substrLength;
		deleteOnExit = false;
	}

  DNALength MakeRCCoordinate(DNALength forPos ) {
    return length - forPos - 1;
  }
  
  void CopyAsRC(DNASequence &rc, DNALength pos=0, DNALength rcLength =0) {
    //
    // Different way of acocunting for position. The position is on
    // the rc strand, not the forward strand.
    //
    if (rcLength == 0) {
      rcLength = length - pos;
    }
    DNALength rcStart = length - (pos + rcLength);
    rc.Resize(rcLength);
    DNALength i;
    for (i = 0; i < rcLength; i++) {
      rc.seq[i] = ReverseComplementNuc[seq[rcStart - 1 + (rcLength - i)]];
    }

    // The reverse complement controls its own memory now.
    rc.deleteOnExit = true;
  }

	void MakeRC(DNASequence &rc, DNALength pos=0, DNALength rcLength=0) {
    if (rcLength == 0) {
      rcLength = length - pos;
    }
    
		rc.Allocate(rcLength);
		DNALength i;
		for (i = 0; i < rcLength; i++) {
			rc.seq[rcLength - i - 1] = ReverseComplementNuc[seq[i+pos]];
		}
    rc.length = rcLength;
    rc.deleteOnExit = true;
	}

	void ToTwoBit() {
		DNALength i;
		for (i = 0; i < length; i++) {
			seq[i] = TwoBit[seq[i]];
		}
		bitsPerNuc = 2;
	}

	void ToThreeBit() {
		DNALength i;
		if (bitsPerNuc != 3) 
			for (i = 0; i < length; i++) { seq[i] = ThreeBit[seq[i]]; }
		bitsPerNuc = 3;
	}

	void ToFourBit() {
		DNALength i;
		if (bitsPerNuc != 4) 
			for (i = 0; i < length; i++) { seq[i] = FourBit[seq[i]]; }
		bitsPerNuc = 4;
	}
	
	void ConvertThreeBitToAscii() {
		DNALength i;
		for (i = 0; i < length; i++ ){
			seq[i] = ThreeBitToAscii[seq[i]];
		}
	}

	void ToAscii() {
		DNALength i;
		if (bitsPerNuc != 8) {
			for (i = 0; i < length; i++ ){ 
				seq[i] = FourBitToAscii[seq[i]];
			}
			bitsPerNuc = 8;
		}
	}
		
	void Assign(DNASequence &ref, DNALength start=0, DNALength plength=0) {
		if (seq != NULL) {
			delete[] seq;
            seq = NULL;
            length = 0;
		}
		if (plength) {
			length = plength;
			seq = new Nucleotide[length];
			memcpy(seq, &ref.seq[start], length);
		}
		else if (start) {
			length = ref.length - start;
			seq = new Nucleotide[length];
			memcpy(seq, &ref.seq[start], length);
		}
		else {
			this->Copy(ref);
		}
		deleteOnExit = true;
	}

	void ToLower() {
		DNALength i;
		for (i = 0; i < length; i++) {
			seq[i] = AllToLower[seq[i]];
		}
	}
		
	void ToUpper() {
		DNALength i;
		for (i = 0; i < length; i++) {
			seq[i] = AllToUpper[seq[i]];
		}
	}

	void Concatenate(const Nucleotide *moreSeq, DNALength moreSeqLength) {
		DNALength prevLength = length;
		length += moreSeqLength;
		Nucleotide *prev = seq;
		seq = new Nucleotide[length];
		if (prev != NULL) {
			memcpy(seq, prev, prevLength);
			delete[] prev;
		}
		memcpy((Nucleotide*) &seq[prevLength], moreSeq, moreSeqLength);

	}

	string GetTitle() const {
		return string("");
	}

	void Concatenate(const Nucleotide* moreSeq) {
		DNALength moreSeqLength = strlen((char*) moreSeq);
		Concatenate(moreSeq, moreSeqLength);
	}
	
	void Concatenate(DNASequence &seq) {
		Concatenate(seq.seq, seq.length);
	}

	int Compare(DNALength pos, DNASequence &rhs, DNALength rhsPos, DNALength length) {
		return memcmp(&seq[pos], &rhs.seq[rhsPos], length);
	}

	int LessThanEqual(DNALength pos, DNASequence &rhs, DNALength rhsPos, DNALength length) {
		int res = Compare(pos, rhs, rhsPos, length);
		if (res <= 0) 
			return 1;
		else
			return 0;
	}

	int Equals(DNASequence &rhs, DNALength rhsPos, DNALength length, DNALength pos=0 ) {
		int res = Compare(pos, rhs, rhsPos, length);
		return res == 0;
	}

	int LessThan(DNALength pos,  DNASequence &rhs, DNALength rhsPos, DNALength length) {
		int res=  Compare(pos, rhs, rhsPos, length);
		return (res < 0);
	}

	void CleanupASCII() {
		DNALength i;
		for (i = 0; i < length; i++ ){
			if (ThreeBit[seq[i]] == 255) {
				seq[i] = 'N';
			}
		}
	}
	Nucleotide operator[](int i) {
		return GetNuc(i);
	}

	Nucleotide GetNuc(DNALength i) {
		return seq[i];
	}

  DNALength GetRepeatContent() {
    DNALength i;
    DNALength nRepeat = 0;
    for (i =0 ; i < length;i++) {
      if (tolower(seq[i]) == seq[i]) { nRepeat++;}
    }
    return nRepeat;
  }

  void CleanupOnFree() {
    deleteOnExit = true;
  }

  void FreeIfControlled() {
    if (deleteOnExit) {
      Free();
    }
  }

	virtual void Free() {
		if (deleteOnExit == false) { return; }
		if (seq != NULL) {
			delete[] seq;
			seq = NULL;
			length = 0;
		}
	}
	void Resize(DNALength newLength) {
		if (seq != NULL) {
			delete[] seq;
		}
		seq = new  Nucleotide[newLength];
		length = newLength;
		deleteOnExit = true;
	}
	DNALength GetSeqStorage() {
		return length;
	}
};

template<typename T>
DNALength ResizeSequence(T &dnaseq, DNALength newLength) {
	assert(newLength > 0);
	if (dnaseq.seq != NULL) {
		delete[] dnaseq.seq;
	}
	dnaseq.seq = new Nucleotide[newLength];
	dnaseq.length = newLength;
	dnaseq.deleteOnExit = true;
	return newLength;
}

#endif
