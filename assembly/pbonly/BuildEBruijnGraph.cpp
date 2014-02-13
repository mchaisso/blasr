//
// Specific to e. Bruijn
//
#include "ebruijn/MatchVertex.h"
#include "ebruijn/ReadWordMatch.h"
#include "ebruijn/ReadDB.h"

//
// General I/O routines.
//
#include "FASTAReader.h"
#include "FASTASequence.h"
#include "algorithms/alignment/readers/sam/SAMReader.h"
#include "datastructures/alignmentset/SAMAlignment.h"
#include "datastructures/alignmentset/SAMToAlignmentCandidateAdapter.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "CommandLineParser.h"
#include "utils/StringUtils.h"

#include <vector>
#include <string>
#include <set>
#include <map>

bool verbose = false;

using namespace std;

bool LookupReadIndex(string name,
                     map<string, int> &readNameMap,
                     int &index) {
  bool result = true;
  map<string,int>::iterator it;
  if ((it = readNameMap.find(name)) == readNameMap.end()) {
    return false;
  }
  else {
    index = (*it).second;
    return true;
  }
}

bool LookupReadIndices(AlignmentCandidate<> &alignment,
                       map<string, int> &readNameMap,
                       int &qIndex, int &tIndex) {
  if (LookupReadIndex(alignment.qName, readNameMap, qIndex) == false) {
    return false;
  }

  if (LookupReadIndex(alignment.tName, readNameMap, tIndex) == false) {
    return false;
  }
  // all good (found both), flag that
  return true;
}

bool ParseReadOffsetFromTitle(string title, int &offset) {
  vector<string> values;
  ParseSeparatedList(title, values, '/');
  if (values.size() < 3) {
    return false;
  }
  stringstream strm(values[2]);
  strm >> offset;
  return true;
}

void MarkMatches(AlignmentCandidate<> &alignment,
                 int vertexSize,
                 ReadWordMatch &qReadWordMatches,
                 ReadWordMatch &tReadWordMatches) {

  int qStrand = alignment.qStrand;
  int tStrand = alignment.tStrand;
  
  if (alignment.blocks.size() > 0) {
    int b;
    int qReadOffset;
    int tReadOffset;
    ParseReadOffsetFromTitle(alignment.qName, qReadOffset);
    ParseReadOffsetFromTitle(alignment.tName, tReadOffset);

    for (b = 0; b < alignment.blocks.size(); b++) {
      if (alignment.blocks[b].length >= vertexSize) {
        int bi;
        int qPos, tPos;
        qPos = alignment.qAlignedSeqPos + alignment.qPos + alignment.blocks[b].qPos;// - qReadOffset;
        tPos = alignment.tAlignedSeqPos + alignment.tPos + alignment.blocks[b].tPos;// - tReadOffset;
        if (qStrand == 1) {
          qPos = qReadWordMatches.size() - qPos - 1;
          for (bi = 0; bi < alignment.blocks[b].length ; bi++, qPos--, tPos++) {
           assert(qPos < qReadWordMatches.pos.size());
           assert(tPos < tReadWordMatches.pos.size());
           qReadWordMatches.SetMatch(qPos);
           tReadWordMatches.SetMatch(tPos);
          }
        }
        else { 
          //
         // Now mark the positions where there are matches.
         for (bi = 0; bi < alignment.blocks[b].length ; bi++, tPos++) {
           assert(qPos < qReadWordMatches.pos.size());
           assert(tPos < tReadWordMatches.pos.size());
           qReadWordMatches.SetMatch(qPos);
           tReadWordMatches.SetMatch(tPos);
         }
        }
      }
    }
  }      
  /*
  cout << alignment.qName << endl;
  qReadWordMatches.PrintPos(cout);
  cout << alignment.tName << endl;
  tReadWordMatches.PrintPos(cout);
  cout << endl;
  */
}

void FindParent(vector<int> &parentIndices, int &cur) {
  int start = cur;
  while (parentIndices[cur] != cur) {
    cur = parentIndices[cur];
  }
  int parent = cur;
  int toPromote = start;
  //
  // Move the pointers up for faster search later.
  //
  while (parentIndices[toPromote] != parent) {
    int temp = toPromote;
    toPromote = parentIndices[toPromote];
    parentIndices[temp] = parent;
  }
}


void Promote(vector<int> &parentIndices, int index) {
  int i = index;
  while (parentIndices[i] != i) {
    i = parentIndices[i];
  }
  int parent = i;
  i = index;
  while (parentIndices[i] != i) {
    int temp = parentIndices[i];
    parentIndices[i] = parent;
    i = temp;
  }
}

void PromoteAll(vector<int> &parentIndices) {
  int i;
  for (i = 1; i < parentIndices.size(); i++) {
    Promote(parentIndices, i);
  }
}


void JoinVertices(AlignmentCandidate<> &alignment,
                  int vertexSize,
                  ReadWordMatch &qReadWordMatches,
                  ReadWordMatch &tReadWordMatches,
                  int &curParentIndex,
                  vector<int> &parentIndices) {

  int qStrand = alignment.qStrand;
  int tStrand = alignment.tStrand;
  /*
  cout << alignment.qName << endl;
  cout << alignment.tName << endl; 
  */
  if (alignment.blocks.size() > 0) {
    int b;
    int qReadOffset;
    int tReadOffset;
    ParseReadOffsetFromTitle(alignment.qName, qReadOffset);
    ParseReadOffsetFromTitle(alignment.tName, tReadOffset);

    for (b = 0; b < alignment.blocks.size(); b++) {
      if (alignment.blocks[b].length >= vertexSize) {
        int bi;
        int qPos, tPos;
        qPos = alignment.qAlignedSeqPos + alignment.qPos + alignment.blocks[b].qPos;// - qReadOffset;
        tPos = alignment.tAlignedSeqPos + alignment.tPos + alignment.blocks[b].tPos;
        if (qStrand == 1) {
          qPos = qReadWordMatches.size() - qPos - 1;
        }
        //
        // Now mark the positions where there are matches.
        for (bi = 0; bi < alignment.blocks[b].length ; bi++, tPos++) {
          //
          // These vertices have not yet been joined to anything.
          // Do this now.
          //
          assert(qPos < qReadWordMatches.parents.size());
          assert(tPos < tReadWordMatches.parents.size());
          assert(curParentIndex < parentIndices.size());
          //          cout << "joining " << qPos << " " << tPos << endl;
          if (qReadWordMatches.parents[qPos] == 0 and tReadWordMatches.parents[tPos] == 0) {
            qReadWordMatches.parents[qPos] = tReadWordMatches.parents[tPos] = curParentIndex;
            parentIndices[curParentIndex] = curParentIndex;
            ++curParentIndex;
          }
          else {
            //
            // The first two cases the parents are not joined
            if (qReadWordMatches.parents[qPos] == 0 and tReadWordMatches.parents[tPos] != 0) {
              qReadWordMatches.parents[qPos] = tReadWordMatches.parents[tPos];
            }
            else if (qReadWordMatches.parents[qPos] != 0 and tReadWordMatches.parents[tPos] == 0) {
              tReadWordMatches.parents[tPos] = qReadWordMatches.parents[qPos];
            }
            else {
              //
              // Both of these vertices have previously been assigned
              // a parent.  Create a new one. 
              // 
              int qParent = qReadWordMatches.parents[qPos];
              FindParent(parentIndices, qParent);

              int tParent = tReadWordMatches.parents[tPos];
              FindParent(parentIndices, tParent);
              if (qParent != tParent) {
                parentIndices[curParentIndex] = curParentIndex;

                parentIndices[qParent] = curParentIndex;
                parentIndices[tParent] = curParentIndex;

                qReadWordMatches.parents[qPos] = curParentIndex;
                tReadWordMatches.parents[tPos] = curParentIndex;
                ++curParentIndex;
              }
            }
          }
          if (qStrand == 0) {
                    qPos++;
          }
          else {
                    qPos--;
          }
       } 
      }
    }
  }
  /*
  cout << alignment.qName << endl;
  qReadWordMatches.PrintParents(cout);
  cout << alignment.tName << endl;
  tReadWordMatches.PrintParents(cout);
  cout << endl;
  */
}


bool MarkMatches(AlignmentCandidate<> &alignment,
                 map<string, int> &readNameMap,
                 int vertexSize,
                 ReadWordMatchVector &readWordMatches) {
  int qIndex, tIndex;
  if (LookupReadIndices(alignment, readNameMap, qIndex, tIndex) == false) {
    return false;
  }
  
  MarkMatches(alignment, vertexSize, readWordMatches[qIndex], readWordMatches[tIndex]);
  return true;
}

bool JoinVertices(AlignmentCandidate<> &alignment,
                  int vertexSize,
                  map<string, int> &readNameMap,
                  ReadWordMatchVector &readWordMatches,                  
                  int &curParentIndex,
                  vector<int> &parentIndices) {
  int qIndex, tIndex;
  if (LookupReadIndices(alignment, readNameMap, qIndex, tIndex) == false) {
    return false;
  }

  JoinVertices(alignment, 
               vertexSize,
               readWordMatches[qIndex], readWordMatches[tIndex],
               curParentIndex, parentIndices);
  return true;
}
                 
//
// This is global to all methods.
//
int vertexSize = 12;

//
// Use this to strip the coordinates from a read to make sure subreads
// are not aligned to each other.
//
bool GetReadBase(string &readName, string &base) {
  int e = readName.rfind('/');
  if (e != readName.npos) {
    base = readName.substr(0,e);
    return true;
  }
  else {
    return false;
  }
}

int main(int argc, char* argv[]) {
  
  CommandLineParser clp;
  string readsFileName;
  string alignmentsFileName;
  string outputFileName;
  float minMergeIdentity = 0.70;
  clp.RegisterStringOption("reads", &readsFileName, "Reads used for alignments.");
  clp.RegisterStringOption("alignments", &alignmentsFileName, "SAM formatted alignments.");
  clp.RegisterIntOption("k", &vertexSize, "Minimum match length", CommandLineParser::PositiveInteger);
  clp.RegisterStringOption("outfile", &outputFileName, "Alignment output.");
  clp.RegisterPreviousFlagsAsHidden();
  clp.RegisterFlagOption("v", &verbose, "");
  clp.RegisterFloatOption("minMergeIdentity", 
                          &minMergeIdentity, 
                          "Minimum identity to merge paths.", CommandLineParser::PositiveFloat);
  
  clp.ParseCommandLine(argc, argv);

  if (minMergeIdentity < 0 or minMergeIdentity > 1) {
    cout << "ERROR. minMergeIdentity must be between 0 and 1" << endl;
    exit(0);
  }
  
  vector<FASTASequence> reads;

  FASTAReader fastaReader;
  fastaReader.Initialize(readsFileName);
  fastaReader.ReadAllSequences(reads);

  //
  // It is necessary to go from read title to index in the list of reads. 
  //
  map<string, int> readNameToIndex;
  BuildReadNameToIndexMap(reads, readNameToIndex);

  ReadWordMatchVector readWordMatches;
  InitializeFromReads(reads, readWordMatches);
  
  //
  // Get ready to read in the alignments.
  //
  SAMReader<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> samReader;
  samReader.Initialize(alignmentsFileName);
  AlignmentSet<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> alignmentSet;
  samReader.ReadHeader(alignmentSet);
  
  SAMAlignment samAlignment;
  AlignmentCandidate<> alignment;
  int numAlignedBases = 0;
  int alignmentIndex = 0;
  while ( samReader.GetNextAlignment( samAlignment ) ) {
    vector<AlignmentCandidate<> > alignments;
    SAMAlignmentsToCandidates(samAlignment,
                              reads,
                              readNameToIndex,
                              alignments, false, true);

    int i;
    ++alignmentIndex;
    int a;
    for (a = 0; a < alignments.size();a++) {
      if (alignments[a].qName != alignments[a].tName) {
        MarkMatches(alignments[a], readNameToIndex, vertexSize, readWordMatches);
      }
    }
    if (alignmentIndex % 1000 == 0) {
      cout << alignmentIndex << endl;
    }
  }


  int numMatches = 0;
  int parentIndex = 1;
  int r;
  for (r = 0; r < readWordMatches.size(); r++) {
    readWordMatches[r].CreateParents();
    numMatches += readWordMatches[r].pos.size();
  }

  vector<int> parentIndices;
  parentIndices.resize(2*numMatches + 1);
  fill(parentIndices.begin(), parentIndices.end(), 0);
  //
  // Start indexing off at 1 so that 0 does not need to be treated in
  // a special case.
  //
  int curParentIndex = 1;
  cout << "There are " << numMatches << " matches." << endl;

  samReader.Close();
  samReader.Initialize(alignmentsFileName);
  AlignmentSet<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> alignmentSet2;
  samReader.ReadHeader(alignmentSet2);
  
  numAlignedBases = 0;
  alignmentIndex = 0;
  while ( samReader.GetNextAlignment( samAlignment ) ) {
    vector<AlignmentCandidate<> > alignments;
    SAMAlignmentsToCandidates(samAlignment,
                              reads,
                              readNameToIndex,
                              alignments, false, true);

    int i;
    ++alignmentIndex;
    int a;
    for (a = 0; a < alignments.size();a++) {
      if (alignments[a].qName != alignments[a].tName) {
        JoinVertices(alignments[a], vertexSize, readNameToIndex, readWordMatches, curParentIndex, parentIndices);
      }
    }
    if (alignmentIndex % 1000 == 0) {
      cout << alignmentIndex << endl;
    }
  }
  vector<int> parentCounts;
  parentCounts.resize(parentIndices.size());
  fill(parentCounts.begin(), parentCounts.end(), 0);
  int p;
  PromoteAll(parentIndices);
  int i;
  for (r = 0; r < readWordMatches.size(); r++) {
    for (i = 0; i < readWordMatches[r].parents.size(); i++) {
      readWordMatches[r].parents[i] = parentIndices[readWordMatches[r].parents[i]];
      parentCounts[readWordMatches[r].parents[i]]++;
    }
  }
  /*
  for (i = 0; i < readWordMatches.size(); i++) {
    readWordMatches[i].PrintPos(cout);
    readWordMatches[i].PrintParents(cout);
  }
  */

  map<int,int> hist;
  int numParents = 0;
  for (i = 1; i < parentCounts.size() && parentIndices[i] != 0; i++) {
    if (parentCounts[i] != 0) {
      ++numParents;
    }
    if (hist.find(parentCounts[i]) == hist.end()) {
      hist[parentCounts[i]] = 1;
    }
    else {
      hist[parentCounts[i]]++;
    }
  }

  map<int,int>::iterator histIt;
  cout << " freq count" << endl;
  for(histIt = hist.begin(); histIt != hist.end(); ++histIt) {
    cout << (*histIt).second << " " << (*histIt).first << endl;
  }

  MatchVertexList vertices;
  vertices.resize(numParents);
  cout << "there are " << numParents << " parents. " << endl;
  
}
