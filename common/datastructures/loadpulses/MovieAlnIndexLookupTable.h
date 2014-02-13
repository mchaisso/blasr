/*
 * =====================================================================================
 *
 *       Filename:  MovieAlnIndexLookupTable.h
 *
 *    Description:  For cpp/pbihdfutils/LoadPulses.cpp use only.
 *                  Cache all the information required for loading metric from pls/bas
 *                  to cmp file.
 *                  mapping movieAlignmentIndex in movieIndexSets[movieIndex] to:
 *                  alignmentIndex: index of this alignment in cmpFile.alnInfo.alignments
 *                  refGroupIndex : index of this alignment in cmpReader.refAlignGroups
 *                  readGroupIndex: index of this alignment in cmpReader.refAlignGroups\
 *                                  [refGroupIndex]->readGroupds[readGroupIndex]
 *                  offsetBegin   : offset begin for this alignment in cmpFile 
 *                                  dataset /ref/movie/AlnArray
 *                  offsetEnd     : offset enda
 *                  queryStart    : start position of query of this alignment 
 *                  queryEnd      : end position of query of this alignment 
 *                  readIndex     : index of this alignment in baseFile.readStartPositions
 *                  readStart     : start position of this alignment in baseFile
 *                  readLength    : read length of this alignment in baseFile
 *                  plsReadIndex  : index of this alignment in pulseFile.pulseStartPositions
 *                  alignedSequence
 *                                : aligned sequence of this alignment
 *                  
 *        Version:  1.0
 *        Created:  12/19/2012 03:50:21 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#ifndef DATASTRUCTURES_LOADPULSES_MovieAlnIndexLookupTable_H_
#define DATASTRUCTURES_LOADPULSES_MovieAlnIndexLookupTable_H_
#include "Types.h"

class MovieAlnIndexLookupTable {
public: 
    bool  skip;                   
    // as movies may split into multiple parts, skip=true if 
    // this alignment is not found in the movie
     
    int   movieAlignmentIndex;    
    // movieIndexSets[movieIndex][toFrom[movieAlignmentIndex]] 
    
    UInt  alignmentIndex;         
    // cmpFile.alnInfo.alignments[alignmentIndex]
    
    UInt  holeNumber;
    // holeNumber corresponding to this alignment in baseFile
    
    int   refGroupIndex;          
    // cmpReader.refAlignGroups[refGroupIndex]
    
    int   readGroupIndex;         
    // cmpReader.refAlignGroups[refGroupIndex]->readGroups[readGroupIndex]
    
    UInt  offsetBegin, offsetEnd; 
    // offset begin and end for this alignment in /ref/movie/AlnArray
    // = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetBegin/End()
    
    int   queryStart,  queryEnd;  
    // start/end for this read = cmpFile.alnInfo.
    // alignments[alignmentIndex].GetQueryStart/End()

    int   readIndex;              
    // index of this alignment in baseFile.readStartPositions
    // = index of this hole number in BaseCalls/ZMW/HoleNumber
    // baseFile.LookupReadIndexByHoleNumber(holeNumber, out=readIndex)
    
    int   readStart;              
    // start pos of this alignment in baseFile
    // = baseFile.readStartPositions[readIndex]

    int   readLength;             
    // read length of this alignment in baseFile
    // = baseFile.readStartPositions[readIndex+1] - readStart 

    int   plsReadIndex;
    // index of this alignment in pulseFile.pulseStartPositions
    // = index of this hole number in PulseCalls/ZMW/HoleNumbers
    // = pulseFile.LookupReadIndexByHoleNumber(holeNumber, out=plsReadIndex)

    // vector<int> baseToAlignmentMap; 
    // keep all the baseToAlignmentMap in memory for now
    // Note that baseToAlignmentMap is not initialized when 
    // BuildLookupTable is called.
   
    string alignedSequence;

    MovieAlnIndexLookupTable() {
        skip = true;
        //skip this alignment unless we can find all the information
    }

    void SetValue(const bool & skipP,            
                  const int  & movieAlignmentIndexP,
                  const UInt & alignmentIndexP,  
                  const int  & refGroupIndexP,   
                  const int  & readGroupIndexP,  
                  const UInt & holeNumberP,
                  const UInt & offsetBeginP,     
                  const UInt & offsetEndP,
                  const int  & queryStartP,      
                  const int  & queryEndP,
                  const int  & readIndexP,
                  const int  & readStartP,       
                  const int  & readLengthP,
                  const int  & plsReadIndexP) {
        skip = skipP; 
        movieAlignmentIndex = movieAlignmentIndexP;
        alignmentIndex      = alignmentIndexP;
        refGroupIndex       = refGroupIndexP;
        readGroupIndex      = readGroupIndexP;
        holeNumber          = holeNumberP;
        offsetBegin         = offsetBeginP;
        offsetEnd           = offsetEndP;
        queryStart          = queryStartP;
        queryEnd            = queryEndP;
        readIndex           = readIndexP;
        readStart           = readStartP;
        readLength          = readLengthP;
        plsReadIndex        = plsReadIndexP;
    }

    void print() {
        // Print this lookup table for debugging . 
        if (skip) 
            cout << "skip = True, ";
        else 
            cout << "skip = False, ";
        cout << "movieAlnIndex    = " << alignmentIndex << ", refGroupIndex = " << refGroupIndex
             << ", readGroupIndex = " << readGroupIndex << ", holeNumber    = " << holeNumber
             << ", offsetBegin    = " << offsetBegin    << ", offsetEnd     = " << offsetEnd
             << ", queryStart     = " << queryStart     << ", queryEnd      = " << queryEnd
             << ", readIndex      = " << readIndex      << ", readStart     = " << readStart
             << ", readLength     = " << readLength     << ", plsReadIndex  = " << plsReadIndex
             << endl;
    }
};
#endif
