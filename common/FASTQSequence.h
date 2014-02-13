#ifndef FASTQ_SEQUENCE_H_
#define FASTQ_SEQUENCE_H_

#include "FASTASequence.h"
#include "NucConversion.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "qvs/QualityValue.h"
#include "qvs/QualityValueVector.h"
#include "datastructures/matrix/Matrix.h"
#include "cmpseq/CompressedSequence.h"
using namespace std;


class FASTQSequence : public FASTASequence {
 public:
  static int charToQuality;
  QualityValueVector<QualityValue> qual;
	QualityValueVector<QualityValue> deletionQV;
	QualityValueVector<QualityValue> preBaseDeletionQV;
	QualityValueVector<QualityValue> insertionQV;
	QualityValueVector<QualityValue> substitutionQV;
	QualityValueVector<QualityValue> mergeQV;
	Nucleotide *deletionTag;
	Nucleotide *substitutionTag;
	int subreadStart, subreadEnd;
	QualityValue deletionQVPrior, insertionQVPrior, substitutionQVPrior, preBaseDeletionQVPrior;
  
  QVScale qvScale;

  QVScale GetQVScale() {
    return qvScale;
  }

  void SetQVScale(QVScale qvScaleP) {
    qvScale                   = qvScaleP;
    qual.qvScale              = qvScale;
    deletionQV.qvScale        = qvScale;
    preBaseDeletionQV.qvScale = qvScale;
    insertionQV.qvScale       = qvScale;
    substitutionQV.qvScale    = qvScale;
    mergeQV.qvScale           = qvScale;
  }

	QualityValueVector<QualityValue>* GetQVPointerByIndex(int index) {
		if (index == 0) { return &qual; }
		if (index == 1) { return &insertionQV; }
		if (index == 2) { return &deletionQV; }
		if (index == 3) { return &substitutionQV; }
		if (index == 4) { return &mergeQV; }
		return NULL;
	}

	int GetStorageSize() {
		int total = 0;
		int nQV = 0;
		int nTag =0;
		if (!qual.Empty()) {
			nQV++;
		}
		if (!deletionQV.Empty()) { 
			nQV++;
		}
		if (!preBaseDeletionQV.Empty()) {
			nQV+=4;
		}
		if (!insertionQV.Empty()) {
			nQV++;
		}
		if (!substitutionQV.Empty()) {
			nQV++;
		}
    if (!mergeQV.Empty()) {
      nQV++;
    }
		if (deletionTag != NULL) {
			nTag++;
		}
		if (substitutionTag !=NULL) {
			nTag++;
		}
		total = nQV*sizeof(QualityValue)*length + nTag*sizeof(Nucleotide)*length;
		return total + FASTASequence::GetStorageSize();
	}
	
  FASTQSequence() : FASTASequence() {
		deletionTag       = NULL;
		substitutionTag   = NULL;

		//
		// For now assume a prior distribution to be the variation of the human genome.
		//
		deletionQVPrior = 0.001; 
		insertionQVPrior = 0.001; 
		substitutionQVPrior = 0.001;
		preBaseDeletionQVPrior = 0.001;

		subreadStart = subreadEnd = 0;
    qvScale = PHRED;
	}

	QualityValue GetDeletionQV(DNALength pos) {
		assert(pos < ((unsigned int)-1));
		assert(pos < length);
		if (deletionQV.Empty()) {
			return deletionQVPrior;
		}
		else {
			return deletionQV[pos];
		}
	}

	QualityValue GetMergeQV(DNALength pos) {
		assert(pos < ((unsigned int)-1));
		assert(pos < length);
		if (mergeQV.Empty()) {
			return 0;
		}
		else {
			return mergeQV[pos];
		}
	}

	Nucleotide GetSubstitutionTag(DNALength pos) {
		if (substitutionTag == NULL) {
			return 'N';
		}
		assert(pos < ((unsigned int)-1));
		assert(pos < length);
		return substitutionTag[pos];
	}

	Nucleotide GetDeletionTag(DNALength pos) {
		if (deletionTag == NULL) {
			return 'N';
		}
		assert(pos < ((unsigned int)-1));
		assert(pos < length);
		return deletionTag[pos];
	}

	QualityValue GetInsertionQV(DNALength pos) {
		if (insertionQV.Empty()) {
			return insertionQVPrior;
		}
		assert(pos < ((unsigned int)-1));
		assert(pos < length);
		return insertionQV[pos];
	}

	QualityValue GetSubstitutionQV(DNALength pos) {
		if (substitutionQV.Empty()) {
			return substitutionQVPrior;
		}
		assert(pos < ((unsigned int)-1));
		assert(pos < length);
		return substitutionQV[pos];
	}

	QualityValue GetPreBaseDeletionQV(DNALength pos, Nucleotide nuc) {
		if (preBaseDeletionQV.Empty()) {
			return preBaseDeletionQVPrior;
		}
		assert(pos < ((unsigned int)-1));
		assert(pos < length);
		return preBaseDeletionQV[pos*4 + TwoBit[nuc]];
	}

	void ShallowCopy(const FASTQSequence &rhs) {
		qual.ShallowCopy(rhs.qual);
		FASTASequence::ShallowCopy(rhs);
	}
  
	void ReferenceSubstring(const FASTQSequence &rhs) {
    ReferenceSubstring(rhs, 0, rhs.length);
  }

	void ReferenceSubstring(const FASTQSequence &rhs, DNALength pos) {
    ReferenceSubstring(rhs, pos, rhs.length - pos);
  }

	void ReferenceSubstring(const FASTQSequence &rhs, DNALength pos, DNALength substrLength) {
    SetQVScale(rhs.qvScale);
		if (substrLength == 0) {
			substrLength = rhs.length - pos;
		}
    FASTASequence::ReferenceSubstring(rhs,pos,substrLength);
		if (rhs.qual.Empty() == false) {
			qual.ShallowCopy(rhs.qual, pos);
		}
		if (rhs.deletionQV.Empty() == false) {
			deletionQV.ShallowCopy(rhs.deletionQV, pos);
		}
		if (rhs.mergeQV.Empty() == false) {
			mergeQV.ShallowCopy(rhs.mergeQV, pos);
		}
		if (rhs.insertionQV.Empty() == false) {
			insertionQV.ShallowCopy(rhs.insertionQV, pos);
		}
		if (rhs.preBaseDeletionQV.Empty() == false ){
			preBaseDeletionQV.ShallowCopy(rhs.preBaseDeletionQV, pos);
		}
		if (rhs.deletionTag != NULL) {
			deletionTag = &rhs.deletionTag[pos];
		}
		if (rhs.substitutionTag != NULL) {
			substitutionTag = &rhs.substitutionTag[pos];
		}
		if (rhs.substitutionQV.Empty() == false) {
			substitutionQV.ShallowCopy(rhs.substitutionQV, pos);
		}
		deletionQVPrior = rhs.deletionQVPrior;
		insertionQVPrior = rhs.insertionQVPrior;
		substitutionQVPrior = rhs.substitutionQVPrior;
		preBaseDeletionQVPrior = rhs.preBaseDeletionQVPrior;
	}

	void ClearAndNull(QualityValue *value) {
		if (value != NULL) {
			delete[] value;
		}
		value = NULL;
	}
	void CopyQualityValues(const FASTQSequence &rhs) {
    SetQVScale(rhs.qvScale);
    qual.Copy(rhs.qual, rhs.length);
    deletionQV.Copy(rhs.deletionQV, rhs.length);
    insertionQV.Copy(rhs.insertionQV, rhs.length);
    substitutionQV.Copy(rhs.substitutionQV, rhs.length);
    mergeQV.Copy(rhs.mergeQV, rhs.length);
    //
    // Handle the tags separtely (and verbosely)
    //
		if (rhs.deletionTag) {
			AllocateDeletionTagSpace(rhs.length);
			memcpy(deletionTag, rhs.deletionTag, sizeof(Nucleotide)*rhs.length);
		}
		else { 
			ClearAndNull(deletionTag);
		}

		if (rhs.substitutionTag) {
			AllocateSubstitutionTagSpace(rhs.length);
			memcpy(substitutionTag, rhs.substitutionTag, sizeof(Nucleotide)*rhs.length);
		}
		else {
			ClearAndNull(substitutionTag);
		}
	}

	void AllocateQualitySpace(DNALength qualLength) {
    qual.Allocate(qualLength);
	}

	void AllocateDeletionQVSpace(DNALength qualLength) {
    deletionQV.Allocate(qualLength);
	}
	
  void AllocateMergeQVSpace(DNALength len) {
    mergeQV.Allocate(len);
  }

	void AllocateDeletionTagSpace(DNALength qualLength) {
		if (deletionTag != NULL) delete[] deletionTag;
		deletionTag = new Nucleotide[qualLength];
	}

	void AllocatePreBaseDeletionQVSpace(DNALength qualLength) {
    preBaseDeletionQV.Allocate(qualLength);
	}
	
	void AllocateInsertionQVSpace(DNALength qualLength) {
    insertionQV.Allocate(qualLength);
	}

	void AllocateSubstitutionQVSpace(DNALength qualLength ){ 
    substitutionQV.Allocate(qualLength);
	}
	
	void AllocateSubstitutionTagSpace(DNALength qualLength ){ 
		if (substitutionTag != NULL) delete[] substitutionTag;
		substitutionTag = new Nucleotide[qualLength];
	}
	
	void AllocateRichQualityValues(DNALength qualLength) {
		AllocateDeletionQVSpace(qualLength);
		AllocateDeletionTagSpace(qualLength);
		AllocatePreBaseDeletionQVSpace(qualLength);
		AllocateInsertionQVSpace(qualLength);
		AllocateSubstitutionQVSpace(qualLength);
		AllocateSubstitutionTagSpace(qualLength);
    AllocateMergeQVSpace(qualLength);
	}
		
	void Copy(const FASTQSequence &rhs) {
		FASTASequence::Copy(rhs);
		CopyQualityValues(rhs);
	}

	FASTQSequence& operator=(const FASTQSequence &rhs) {
		this->Copy(rhs);
		return *this;
	}

	FASTQSequence(const FASTQSequence &rhs) {
		substitutionTag = NULL;
		deletionTag = NULL;
		this->Copy(rhs);
	}

	void Assign(FASTQSequence &rhs) {
		// copy the nucleotide part
		FASTASequence::Assign(rhs);
		// copy the qual part
		CopyQualityValues(rhs);
    SetQVScale(rhs.qvScale);
	}
  
  void PrintFastq(ostream &out, int lineLength=50) {
    PrintSeq(out, lineLength, '@');
    if (lineLength == 0) { 
      out << endl;
    }
    PrintFastqQuality(out, lineLength);
    if (lineLength == 0) {
      out << endl;
    }
  }

  void PrintFastqQuality(ostream &out, int lineLength=50) {
    out << "+" << endl;
    PrintAsciiQual(out, lineLength);
  }

	void PrintAsciiRichQuality(ostream &out, int whichQuality, int lineLength=50) {
		unsigned char* qualPtr;
		int charOffset = charToQuality;
		if (whichQuality == 0) {
			qualPtr = qual.data;
		}
		else if (whichQuality == 1) {
			qualPtr = insertionQV.data;
		}
		else if (whichQuality == 2) {
			qualPtr = deletionQV.data;
		}
		else if (whichQuality == 3) {
			qualPtr = substitutionQV.data; 
		}
		else if (whichQuality == 4) {
			qualPtr = mergeQV.data;
		}
		else if (whichQuality == 5) {
			qualPtr = (unsigned char*) substitutionTag;
			charOffset = 0;
		}
		else if (whichQuality == 6) {
			qualPtr = (unsigned char*) deletionTag;
			charOffset = 0;
		}
		int i;
    if (lineLength == 0) {
      for (i = 0; i < length; i++) {
				if (qualPtr != NULL) {
					out << (char)(qualPtr[i]+charOffset);
				}
				else {
					// Fake bad quality
					out << "5";
				}
      }
    }
    else {
      for (i = 0; i < length; i++) {
				assert(((unsigned char) (qualPtr[i] + charOffset) > 32) &&
							 ((unsigned char) (qualPtr[i] + charOffset) < 127));
				if (qualPtr != NULL) {
					out << (char)(qualPtr[i]+charOffset);
				}
				else {
					// Fake pretty bad quality.
					out << "5";
				}
        assert(lineLength != 0);
        if (i > 0 and (i+1) % lineLength==0) {
          out << endl;
        }
      }
      if (i == 0 or i % lineLength != 0) {
        out << endl;
      }
    }

	}
  
  void PrintAsciiQual(ostream &out, int lineLength=50) {
		PrintAsciiRichQuality(out, 0, lineLength);
  }


  void PrintQual(ostream &out, int lineLength = 50) {
		out << ">" << this->title << endl;
    DNALength i;
		for (i = 0; i < length; i++ ){
			out << (int) qual[i];
			if (i > 0 and (i+1) % lineLength == 0)
				out << endl;
			else 
				out << " ";
		}
		if (i == 0 or i % lineLength != 0) {
			out << endl;
		}
  }

	void PrintQualSeq(ostream &out, int lineLength = 50) {
		FASTASequence::PrintSeq(out, lineLength);
		int i;
		lineLength /= 4;
    PrintQual(out, lineLength);
	}

	void MakeRC(FASTQSequence &rc) {
    rc.SetQVScale(qvScale);
    FASTASequence::MakeRC(rc);
		if (qual.Empty() == true) {
			// there is no quality values, don't make an rc.
			return;
		}

		if (rc.qual.Empty() == false) {
      rc.qual.Free();
		}
		rc.AllocateQualitySpace(length);
		int i;
		for (i = 0; i < length; i++ ){
			rc.qual.data[length - i - 1] = qual[i];
		}

		if (deletionQV.Empty() == false) {
			//
			// The read contains rich quality values. Reverse them here.
			//
			rc.AllocateRichQualityValues(length);
			DNALength pos;
      
			for (pos = 0; pos < length; pos++) {
				if (insertionQV.Empty() == false) {
					rc.insertionQV[length - pos - 1] = insertionQV[pos];
				}
				if (substitutionQV.Empty() == false) {
					rc.substitutionQV[length - pos - 1]           = substitutionQV[pos];
				}
				if (deletionQV.Empty() == false) {
					rc.deletionQV[length - pos - 1]               = deletionQV[pos];
				}

				if (mergeQV.Empty() == false) {
					rc.mergeQV[length - pos - 1]               = mergeQV[pos];
				}
        

				if (substitutionTag != NULL) {
					rc.substitutionTag[length - pos - 1] = ReverseComplementNuc[substitutionTag[pos]];
				}
				if (deletionTag != NULL) {
					rc.deletionTag[length - pos - 1]     = ReverseComplementNuc[deletionTag[pos]];
				}
			}
		}
		deletionQVPrior = rc.deletionQVPrior;
		insertionQVPrior = rc.insertionQVPrior;
		substitutionQVPrior = rc.substitutionQVPrior;
		preBaseDeletionQVPrior = rc.preBaseDeletionQVPrior;
	}

	void Free() {
		FASTASequence::Free();
    if (deleteOnExit == true) {
      qual.Free();
      deletionQV.Free();
      preBaseDeletionQV.Free();
      insertionQV.Free();
      substitutionQV.Free();
      mergeQV.Free();
      if (deletionTag != NULL) {
        delete[] deletionTag;
        deletionTag = NULL;
      }
      if (substitutionTag != NULL) {
        delete[] substitutionTag;
        substitutionTag = NULL;
      }
    }
	}

	void LowerCaseMask(int qThreshold) {
		int i;
		if (qual.Empty() == true) return;

		for (i = 0; i < length; i++ ){
			if (qual[i] < qThreshold) {
				seq[i] = tolower(seq[i]);
			}
		}
	}

	float GetAverageQuality() {
		DNALength p;
		float totalQ;
		if (qual.Empty() == true) { return 0.0; }
		assert(qual.Empty() == false);
		assert(length > 0);
		for (p = 0, totalQ = 0.0; p < length; p++) {
			totalQ += qual[p];
		}
		return totalQ / length;
	}

};

//
// Initialize a read with quality probabilities from one with quality values.
//
int FASTQSequence::charToQuality = 33;


#endif
