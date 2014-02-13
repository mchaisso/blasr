#ifndef DATASTRUCTURES_READS_PULSE_FILE_H_
#define DATASTRUCTURES_READS_PULSE_FILE_H_

#include <vector>
#include "Types.h"
#include "Enumerations.h"
#include "data/hdf/PlatformId.h"
#include "PulseBaseCommon.h"
#include "ScanData.h"

class PulseFile : public PulseBaseCommon {
 public:
    unsigned int numFrames;
    PlatformId platformId;
    vector<unsigned int> startFrame;
    vector<uint16_t> plsWidthInFrames;
    int midSignalNDims, maxSignalNDims, meanSignalNDims;
    vector<uint16_t> midSignal;
    vector<uint16_t> maxSignal;
    vector<uint16_t> meanSignal; 
    vector<int>      numEvent;
    vector<int>      pulseStartPositions;
    vector<float>    classifierQV;

    void CopySignal(HalfWord *signalData, // either a vector or matrix
                    int signalNDims,
                    int pulseStartPos,    // 0 if baseToPulseIndex maps to abs position
                    int *baseToPulseIndex,
                    Nucleotide *readSeq,
                    int readLength,
                    HalfWord *readData) {
        // if baseToPulseIndex maps bases to absolute pulse postions
        // pulseStartPos must be 0; 
        // otherwise, pulseStartPos is pulseStartPositions[holeIndex]

      map<char, int> baseMap = GetBaseMap();
      int i;
      if (signalNDims == 1) {
        for (i = 0; i < readLength; i++) {
          readData[i] = signalData[pulseStartPos + baseToPulseIndex[i] ];
        }
      }
      else {
        for (i = 0; i < readLength; i++) {
          readData[i] = signalData[baseToPulseIndex[i]*4 + baseMap[readSeq[i]]];
        }
      }
    }

    template<typename T_FieldType>
      void StoreField(vector<T_FieldType> &source, int *basToPlsIndex, T_FieldType *dest, int destLength) {
      int i;
      for (i = 0 ; i < destLength; i++) {
        dest[i] = source[basToPlsIndex[i]];
      }
    }
    
    template <typename T>
      bool Realloc(T *&ptr, int length) {
      if (ptr != NULL) {
        delete[] ptr;
      }
      ptr = new T[length];
      return ptr != NULL;
    }


    // plsReadIndex: index of this hole number in /PulseCalls/ZMW/HoleNumber.
    // baseToPulseIndex: index from pulse to base from the beginning of the read.
    // read: a SMRTSequence.
    void CopyReadAt(uint32_t plsReadIndex, int *baseToPulseIndex, SMRTSequence &read) {
      int pulseStartPos = pulseStartPositions[plsReadIndex];
      bool allocResult;
      if (midSignal.size() > 0) {
        assert(midSignal.size() > pulseStartPos);
        allocResult = Realloc(read.midSignal, read.length);
        CopySignal(&midSignal[0], 
                   midSignalNDims, 
                   pulseStartPos,
                   baseToPulseIndex,
                   read.seq, read.length,
                   read.midSignal);
      }

      if (maxSignal.size() > 0) {
        assert(maxSignal.size() > pulseStartPos); 
        Realloc(read.maxSignal, read.length);
        CopySignal(&maxSignal[0], 
                   maxSignalNDims, 
                   pulseStartPos,
                   baseToPulseIndex,
                   read.seq, read.length,
                   read.maxSignal);
        }

      if (meanSignal.size() > 0) {
        assert(meanSignal.size() > pulseStartPos); 
        Realloc(read.meanSignal, read.length);
        CopySignal(&meanSignal[0], 
                   meanSignalNDims, 
                   pulseStartPos,
                   baseToPulseIndex,
                   read.seq, read.length,
                   read.meanSignal);
      }
      if (plsWidthInFrames.size() > 0) {
        Realloc(read.widthInFrames, read.length);
        StoreField(plsWidthInFrames, baseToPulseIndex, read.widthInFrames, read.length);
      }
      if (classifierQV.size() > 0) {
        Realloc(read.classifierQV, read.length);
        StoreField(classifierQV, baseToPulseIndex, read.classifierQV, read.length);
      }
      if (startFrame.size() > 0) {
        Realloc(read.startFrame, read.length);
        StoreField(startFrame, baseToPulseIndex, read.startFrame, read.length);
      }
    }
};

#endif
