#ifndef  SMRT_SEQUENCE_H_
#define  SMRT_SEQUENCE_H_
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>

#include "NucConversion.h"
#include "FASTQSequence.h"
#include "Enumerations.h"
#include "datastructures/reads/RegionTable.h"
#include "datastructures/reads/ZMWGroupEntry.h"
using namespace std;
#include <iostream>
#include <sstream>
#include "Types.h"

typedef unsigned char Nucleotide;

class SMRTSequence : public FASTQSequence {
public:
  int16_t xy[2];
  int holeNumber;
  ZMWGroupEntry zmwData;
  PlatformType platform;
  HalfWord *preBaseFrames;
  HalfWord *widthInFrames;
  //
  // The following are fields that are read in from the pulse file.
  // Because they are not standard in bas.h5 files, these fields
  // should not be preallocated when resizing a SMRTSequence, and
  // memory should be managed separately.  For now, these fields all
  // have the same length as the number of bases, but this could
  // change so that all pulse values are stored in a SMRTSequence.
  //
  HalfWord *meanSignal, *maxSignal, *midSignal;
  float *classifierQV;
  unsigned int *startFrame;
  int *pulseIndex;
  DNALength lowQualityPrefix, lowQualitySuffix;
	unsigned int *IPD;
	float snr[4];
	float readQuality;
	int numPasses;
	int cx;
	int subreadContext;
	int qs, qe;
  void SetNull() {
    pulseIndex    = NULL;
    preBaseFrames = NULL;
    widthInFrames = NULL;
    xy[0] = 0; xy[1] = 0;
    // These are not allocted by default.
    meanSignal = maxSignal = midSignal = NULL;
    classifierQV = NULL;
    startFrame   = NULL;
    platform     = NoPlatformType;
    // By default, allow the entire read.
    lowQualityPrefix = lowQualitySuffix = 0;
		IPD          = NULL;
		numPasses = 0;
		snr[0] = snr[1] = snr[2] = snr[3] = 0;
		subreadContext = 0;
		readQuality = 0;
		subreadContext = 0;
		qs = qe = 0;
  }
  
  SMRTSequence() : FASTQSequence() {
    holeNumber = -1;
    SetNull();
  }

  void Allocate(DNALength length) {
    FASTQSequence::AllocateRichQualityValues(length);
    seq           = new Nucleotide[length];
    qual.Allocate(length);
    preBaseFrames = new HalfWord[length];
    widthInFrames = new HalfWord[length];
    pulseIndex    = new int[length];
    subreadEnd    = length;
    deleteOnExit  = true;
  }

  void SetSubreadTitle(SMRTSequence &subread, DNALength subreadStart, DNALength  subreadEnd) {
    stringstream titleStream;
    titleStream << title << "/"<< subreadStart << "_" << subreadEnd;
    subread.CopyTitle(titleStream.str());
  }    

  void SetSubreadBoundaries(SMRTSequence &subread, DNALength &subreadStart, int &subreadEnd) {
    if (subreadEnd == -1) {
      subreadEnd = length;
    }
    assert(subreadEnd - subreadStart <= length);
    subread.subreadStart= subreadStart;
    subread.subreadEnd  = subreadEnd;
		//
		// Set qe & qs for printing purposes.  Only set qe and qs if the main read is not considered a subread.
		//
		if (qe == 0) {
			subread.qe = subreadEnd+1;
			subread.qs = subreadStart+1;
		}
    SetSubreadTitle(subread, subreadStart, subreadEnd);
  }

  void MakeSubreadAsMasked(SMRTSequence &subread, DNALength subreadStart = 0, int subreadEnd = -1) {
    //
    // This creates the entire subread, but masks out the portions
    // that do not correspond to this insert.
    //
    subread.Copy(*this);
    SetSubreadBoundaries(subread, subreadStart, subreadEnd);
    DNALength pos;
    for (pos = 0; pos < subreadStart; pos++) { subread.seq[pos] = 'N'; }
    for (pos = subreadEnd; pos < length; pos++) { subread.seq[pos] = 'N'; }
    // This is newly allocated memory, free it on exit.
    subread.deleteOnExit = true;
  }

	void ReferenceSubstring(const SMRTSequence &rhs) {
    ReferenceSubstring(rhs, 0, rhs.length);
	}
	void ReferenceSubstring(const SMRTSequence &rhs, DNALength pos) {
		ReferenceSubstring(rhs, pos, rhs.length - pos);
	}

	void ReferenceSubstring(const SMRTSequence &rhs, DNALength pos, DNALength substrLength) {
		((FASTQSequence*)this)->ReferenceSubstring(rhs, pos, substrLength);
		if (rhs.IPD != NULL) {
			IPD = &rhs.IPD[pos];
		}
	}

  void MakeSubreadAsReference(SMRTSequence &subread, DNALength subreadStart = 0, int subreadEnd = -1) {
    //
    // Just create a reference to a substring of this read.  
    //
    SetSubreadBoundaries(subread, subreadStart, subreadEnd);
    subread.ReferenceSubstring(*this, subreadStart, subreadEnd - subreadStart);
    // The subread references this read, protect the memory.
    subread.deleteOnExit = false;
  }

  void Copy(const SMRTSequence &rhs) {
    Copy(rhs, 0, rhs.length);
  }
    
    
  void Copy(const SMRTSequence &rhs, int rhsPos, int rhsLength) {
    //
    // Make sure not attempting to copy into self.
    //
    SMRTSequence subseq;
    subseq.ReferenceSubstring(rhs, rhsPos, rhsLength);
    subseq.title = rhs.title;
    subseq.titleLength = strlen(rhs.title);
    if (rhs.length == 0) {
      if (preBaseFrames != NULL) { 
        delete[] preBaseFrames;
        preBaseFrames = NULL;
      }
      if (widthInFrames != NULL) {
        delete[] widthInFrames;
        widthInFrames = NULL;
      }
      if (pulseIndex != NULL) {
        delete[] pulseIndex;
        pulseIndex = NULL;
      }
			if (IPD != NULL) {
				delete[] IPD;
				IPD = NULL;
			}

      ((FASTQSequence*)this)->Copy(subseq);
      //
      // Make sure that no values of length 0 are allocated by returning here.
      //
    }
    else {
      
      assert(rhs.seq != seq);
      assert(rhsLength <= rhs.length);
      assert(rhsPos < rhs.length);
      
      ((FASTQSequence*)this)->Copy(subseq);
      if (rhs.preBaseFrames != NULL) {
        preBaseFrames = new HalfWord[length];
        memcpy(preBaseFrames, rhs.preBaseFrames, length*sizeof(HalfWord));
      }
      if (rhs.widthInFrames != NULL) {
        widthInFrames = new HalfWord[length];
        memcpy(widthInFrames, rhs.widthInFrames, length*sizeof(HalfWord));
      }
      if (rhs.pulseIndex != NULL) {
        pulseIndex = new int[length];
        memcpy(pulseIndex, rhs.pulseIndex, sizeof(int) * length);
      }
      if (rhs.IPD != NULL) {
        IPD = new unsigned int[length];
        memcpy(IPD, rhs.IPD, sizeof(int) * length);
      }

    }
		subreadStart = rhs.subreadStart;
		subreadEnd   = rhs.subreadEnd;
		qs           = rhs.qs;
		qe           = rhs.qe;
		memcpy(snr, rhs.snr, sizeof(float)*4);
		lowQualityPrefix = rhs.lowQualityPrefix;
		lowQualitySuffix = rhs.lowQualitySuffix;
    zmwData = rhs.zmwData;
  }
  
  void Print(ostream &out) {
    out << "SMRTSequence for zmw " << zmwData.holeNumber
        << ", [" << subreadStart << ", " << subreadEnd << ")" << endl;
    DNASequence::Print(out);
  }

  SMRTSequence& operator=(const SMRTSequence &rhs) {
    Copy(rhs);
    return *this;
  }
  
  void Free() {
    FASTQSequence::Free();
    if (deleteOnExit == true) {
      if (preBaseFrames)  {
        delete[] preBaseFrames;
        preBaseFrames = NULL;
      }
      if (widthInFrames) {
        delete[] widthInFrames;
        widthInFrames = NULL;
      }
      if (pulseIndex) {
        delete[] pulseIndex;
        pulseIndex = NULL;
      }
      if (startFrame) {
        delete[] startFrame;
        startFrame = NULL;
      }
			if (IPD) {
				delete [] IPD;
				IPD = NULL;
			}
      // meanSignal, maxSignal, midSignal and classifierQV
      // need to be handled separatedly.
    }
    xy[0] = 0; xy[1] = 0;
    lowQualityPrefix = lowQualitySuffix = 0;
    holeNumber = -1;
  }
      
	bool StoreIPD() {
		if (length == 0) {
			return false;
		}
		if (preBaseFrames == NULL or widthInFrames == NULL) {
			return false;
		}
		IPD = new unsigned int[length];
		int i;
		assert(preBaseFrames);
		assert(widthInFrames);
		for (i = 0; i < length; i++) { 
			IPD[i] = preBaseFrames[i] + widthInFrames[i];
		}
		return true;
	}
  bool StoreXY(int16_t xyP[]) {
    xy[0] = xyP[0];
    xy[1] = xyP[1];
    return true;
  }

	void MakeRC(SMRTSequence &rc) {
		((FASTQSequence*)this)->MakeRC(rc);

		if (rc.IPD != NULL and rc.deleteOnExit) {
			delete[] rc.IPD;
			rc.IPD = NULL;
		}
		if (IPD == NULL) {
			return;
		}
		rc.IPD = new unsigned int[length];
		for (int i = 0; i < length; i++) {
			rc.IPD[length-i-1] = IPD[i];
		}
	}

  bool StorePlatformType(PlatformId pid ){
    if (pid == AstroPlatform) {
      platform = Astro;
    }
    if (pid == SpringfieldPlatform) {
      platform = Springfield;
    }
  }

  bool StorePlatformType(PlatformType ptype) {
    platform = ptype;
    return true;
  }

  bool StoreHoleNumber(int holeNumberP){ 
		zmwData.holeNumber = holeNumber = holeNumberP;
    return true;
  }
  
  bool StoreHoleStatus(unsigned int s) {
    zmwData.holeStatus = s;
    return true;
  }
  
  bool StoreZMWData(ZMWGroupEntry &data) {
    zmwData = data;
    return true;
  }

  bool GetXY(int xyP[]) {
    xyP[0] = xy[0];
    xyP[1] = xy[1];
    return true;
  }

  bool GetHoleNumber(int& holeNumberP) {
    holeNumberP = holeNumber;
    return true;
  }

	void PrintCSVPulseMetric(ostream &out, int whichQuality, int lineLength=50) {
		unsigned int* dataPtr = NULL;
		int charOffset = charToQuality;
		int i;
		if (whichQuality == 7) {
			dataPtr = IPD;
		}
		else {
			cout << "ERROR. Invalid quality " << whichQuality << endl;
			exit(0);
		}
    if (lineLength == 0) {
      for (i = 0; i < length; i++) {
				if (dataPtr != NULL) {
					out << (dataPtr[i]);
				}
				else {
					// Fake bad quality
					out << "5";
				}
				if (i < length-1) {
					out << ",";
				}
      }
    }
    else {
      for (i = 0; i < length; i++) {
				if (dataPtr != NULL) {
					out << (dataPtr[i]);
				}
				else {
					out << "0";
				}
				if (i < length - 1) {
					out << ",";
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

};

#endif
