#ifndef RGRAPH_H_
#define RGRAPH_H_

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <map>
#include "Read.h"
#include "Types.h"
#include "Enumerations.h"
#include "RM4Line.h"
#include "PBReadNameParser.h"
using namespace std;


class RGraph {
 public:
  typedef map<string, SameMolecule*> MoleculeMap;
  MoleculeMap moleculeMap;
  typedef map<string, Read*> ReadMap;
  ReadMap readMap;
  int lineNumber;
  string line;

  void WriteStringToBinary(const string &str, ostream &outFile) {
    UInt strLen = str.size();
    outFile.write((const char*) &strLen, sizeof(strLen));
    outFile.write(str.c_str(), strLen);
  }

  void ReadStringToBinary(istream &inFile, string &str) {
    UInt strLen;
    inFile.read((char*)&strLen, sizeof(strLen));
    if (strLen > 0) {
      char *tmpStr = new char[strLen+1];
      inFile.read(tmpStr, strLen);
      tmpStr[strLen] = '\0';
      str = tmpStr;
    }
    else {
      str.clear();
    }
  }
      
    
  void WriteMoleculeMap(ostream &outFile) {
    unsigned int mmSize = moleculeMap.size();
    outFile.write((const char*) &mmSize, sizeof(mmSize));
    MoleculeMap::iterator mmIt, mmEnd;
    mmEnd = moleculeMap.end();
    for (mmIt = moleculeMap.begin(); mmIt != mmEnd; ++mmIt) {
      WriteStringToBinary((*mmIt).first, outFile);
      //      (*mmIt).second->Write(outFile);
    }
  }

  void Write(ostream &outFile) {
    


  }

  RGraph() {
    lineNumber = 0;
  }
  
  Read* FetchRead(string &readName) {
    ReadMap::iterator it = readMap.find(readName);
    if (it != readMap.end()) {
      return (*it).second;
    }
    else {
      return NULL;
    }
  }
  
  void PrintPossibleChimeras(int minCoverage) {
    ReadMap::iterator mapIt, mapEnd;
    
    for (mapIt = readMap.begin(), mapEnd = readMap.end(); mapIt != mapEnd; ++mapIt) {
      if ((*mapIt).second->IsLikelyChimeric(minCoverage)) {
        cout << (*mapIt).first << endl;
        (*mapIt).second->PrintLowCoverageIntervals(100000);
      }
    }
  }

  void PrintLowCoverageReads(int minCoverage) {
    ReadMap::iterator mapIt, mapEnd;
    
    for (mapIt = readMap.begin(), mapEnd = readMap.end(); mapIt != mapEnd; ++mapIt) {
      cout << (*mapIt).first << endl;
      (*mapIt).second->PrintLowCoverageIntervals(minCoverage);
    }

  }
      
  void StoreParsedReadInformation(string &qName, 
                                  string &moleculeName, 
                                  UInt &subreadBegin, 
                                  UInt &subreadEnd) {

    //
    // Parse the read name for appropriate information: molecule name
    // and subread subreads.
    //
    if (PBReadNameParser::GetMoleculeNameFromReadName(qName, moleculeName) == false) {
      cout << "Malformatted title at " << lineNumber << endl;
      cout << line << endl;
      return;
    }

    if (PBReadNameParser::GetReadCoordinatesFromReadName(qName,
                                                         subreadBegin,
                                                         subreadEnd) == false) {
      cout << "Malformatted subread name at " << lineNumber << endl;
      cout << line<< endl;
      exit(1);
    }
  }

  Read* InitializeRead(string &name, UInt subreadBegin, UInt subreadEnd, UInt length, UInt fullReadLength) {
    Read *readPtr = new Read;
    readMap[name] = readPtr;
    readPtr->subreadBegin = subreadBegin;
    readPtr->subreadEnd   = subreadEnd;
    readPtr->readLength   = length;
    readPtr->fullReadLength=fullReadLength;

    readPtr->InitializeEndNodes();
    return readPtr;
  }

  UInt FullReadPosToSubreadForwardPos(Read &read, UInt x, Strand strand) {
    if (strand == Forward) {
      return x - read.subreadBegin;
    }
    else {
      UInt forwardPos = read.fullReadLength - x - 1;
      return forwardPos - read.subreadBegin;
    }
  }

  void SetForwardPositionsRelativeToFullread(Read &read, UInt fullReadAlnBegin, UInt fullReadAlnEnd, Strand alnStrand,
                                            UInt &alnBegin, UInt &alnEnd) {

    if (alnStrand == Forward) {
      alnBegin = fullReadAlnBegin;
      alnEnd   = fullReadAlnEnd;
    }
    else {
      alnBegin = read.readLength - fullReadAlnEnd;
      alnEnd   = read.readLength - fullReadAlnBegin;
    }
  }
  
  void SetForwardPositionsRelativeToSubread(Read &read, UInt fullReadAlnBegin, UInt fullReadAlnEnd, Strand alnStrand,
                                            UInt &forwardSubreadAlnBegin, UInt &forwardSubreadAlnEnd) {
    
    if (alnStrand == Forward) {
      forwardSubreadAlnBegin = FullReadPosToSubreadForwardPos(read, fullReadAlnBegin, alnStrand);
      forwardSubreadAlnEnd   = FullReadPosToSubreadForwardPos(read, fullReadAlnEnd, alnStrand);
    }
    else {
      forwardSubreadAlnEnd   = FullReadPosToSubreadForwardPos(read, fullReadAlnBegin, alnStrand);
      forwardSubreadAlnBegin = FullReadPosToSubreadForwardPos(read, fullReadAlnEnd, alnStrand);

    }
  }


  void StoreOverlapPair(Read *query, Read *target, RM4Line rm4Line) {
    //
    // Store the start of the overlap.
    //
    UInt qAlnBegin, qAlnEnd, tAlnBegin, tAlnEnd;
    SetForwardPositionsRelativeToSubread(*query, rm4Line.qBegin, rm4Line.qEnd, (Strand) rm4Line.qStrand, qAlnBegin, qAlnEnd);
    SetForwardPositionsRelativeToFullread(*target, rm4Line.tBegin, rm4Line.tEnd, (Strand) rm4Line.tStrand, tAlnBegin, tAlnEnd);
    
    ReadNode *queryBeginNode, *queryEndNode, *targetBeginNode, *targetEndNode;
    query->AddInterval(qAlnBegin, qAlnEnd, queryBeginNode, queryEndNode);
    target->AddInterval(tAlnBegin, tAlnEnd, targetBeginNode, targetEndNode);
    
    queryBeginNode->AddOverlap(queryEndNode, targetBeginNode, targetEndNode);
    targetBeginNode->AddOverlap(targetBeginNode, queryBeginNode, queryEndNode);

    queryBeginNode->AddNodeToSuperset(targetBeginNode);
    queryEndNode->AddNodeToSuperset(targetEndNode);
  }

  void StoreOverlap(RM4Line &rm4Line) {
    Read* readPtr = FetchRead(rm4Line.qName);
    string qMoleculeName, tMoleculeName;
    UInt qSubreadBegin, qSubreadEnd, tSubreadBegin, tSubreadEnd;
    
    StoreParsedReadInformation(rm4Line.qName, qMoleculeName,
                               qSubreadBegin, qSubreadEnd);
    
    if (readPtr == NULL) {
      readPtr = InitializeRead(rm4Line.qName, 
                               qSubreadBegin, qSubreadEnd, 
                               qSubreadEnd - qSubreadBegin, rm4Line.qLength);
    }

    if (rm4Line.qName == rm4Line.tName) {
      //
      // This is a self-alignment.  Don't try and assemble the read.
      //
      return;
    }

    StoreParsedReadInformation(rm4Line.tName, tMoleculeName,
                               tSubreadBegin, tSubreadEnd);
    
    Read* targetPtr = FetchRead(rm4Line.tName);
    if (targetPtr == NULL) {
      targetPtr = InitializeRead(rm4Line.tName,
                                 tSubreadBegin, tSubreadEnd,
                                 tSubreadEnd - tSubreadBegin,
                                 rm4Line.tLength);
    }

    if (readPtr->OverlapsWith(targetPtr)) {
      //
      // No need to do symmetric add (may need to look into this
      // later. 
      //
      return;
    }

    /*    
    if (qMoleculeName == tMoleculeName) {
      //
      // For now, don't process same-molecule alignments.
      //
      return;
    }
    */

    //
    // Process the overlap.
    //
    StoreOverlapPair(readPtr, targetPtr, rm4Line);

  }

  int CountConsistentSupernodes() {
    int numConsistent = 0;
    ReadMap::iterator mapIt, mapEnd;
    for (mapIt = readMap.begin(), mapEnd = readMap.end(); 
         mapIt != mapEnd; ++mapIt) {
      numConsistent+= (*mapIt).second->CountConsistentSupersets();
    }
    return numConsistent;
  }

  void MergeOverlaps() {

    ReadMap::iterator mapIt, mapEnd;
    for (mapIt = readMap.begin(), mapEnd = readMap.end(); mapIt != mapEnd; ++mapIt) {
      mapIt->second->MergeOverlaps();
    }
  }

  void RemoveLowCoverage(int minCoverage) {

    ReadMap::iterator mapIt, mapEnd;
    bool graphIsModified = true;
    int iter = 1;
    while (graphIsModified) { 
      graphIsModified = false;
      int numberRemovedReads = 0;
      for (mapIt = readMap.begin(), mapEnd = readMap.end(); mapIt != mapEnd; ++mapIt) {
        if ((*mapIt).second->RemoveIfLowCoverage(minCoverage)) {
          graphIsModified = true;
          numberRemovedReads++;
        }
      }
      cout << "iter " << iter << " removed " << numberRemovedReads << endl;
      ++iter;
    }
  }
  
  void BuildFromRM4(string &rm4FileName) {
    ifstream in;
    in.open(rm4FileName.c_str());
    if (! in ) {
      cout << "Cannot open " << rm4FileName << endl;
      exit(1);
    }
    lineNumber = 0;
    while (in) {
      RM4Line rm4Line;
      line;
      getline(in, line);
      if (line.size() == 0) {
        break;
      }
      stringstream strm(line);
      strm >> rm4Line;
      string moleculeName;
      UInt subreadBegin, subreadEnd;
      StoreParsedReadInformation(rm4Line.qName, 
                                 moleculeName, 
                                 subreadBegin, subreadEnd);

      //
      // Now process overlap.
      //
      StoreOverlap(rm4Line);
      ++lineNumber;
      if (lineNumber % 10000 == 0) {
        cout << "." << flush;
      }
      if (lineNumber % 500000 == 0) {
        cout << endl;
      }
    }
    cout << endl;
  }
};

#endif
