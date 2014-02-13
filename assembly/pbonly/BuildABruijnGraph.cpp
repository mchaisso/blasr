#include "FASTAReader.h"
#include "FASTASequence.h"
#include "algorithms/alignment/readers/sam/SAMReader.h"
#include "datastructures/alignmentset/SAMAlignment.h"
#include "datastructures/alignmentset/SAMToAlignmentCandidateAdapter.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "CommandLineParser.h"

#include <vector>
#include <string>
#include <set>

bool verbose = false;

using namespace std;


class NucleotideVertex {
public:
  vector<int> out;
  void AddVertex(int v) {
    if (binary_search(out.begin(), out.end(), v)) {
      return;
    }
    else {
      out.push_back(v);
      sort(out.begin(), out.end());
    }
  }
};


typedef vector<int> VertexList;

class AGVertex {
public:
  vector<int> in, out;
  int na, nc, ng, nt;

  AGVertex() {
    na = nc = ng = nt = 0;
  }
};

class AGEdge {
public:
  int coverage;
  int ends[2];
  AGEdge() {
    coverage = 0;
  }
};


int TransformAListToVertexIndices(vector<int> &aList, int numBases) {
  assert(numBases <= aList.size());
  
  map<int, int> graphVertices;
  int a;
  int curVertex = 0;
  for (a = numBases; a < aList.size(); a++) {
    if (graphVertices.find(aList[a]) == graphVertices.end()) {
      graphVertices[aList[a]] = curVertex;
      curVertex++;
    }
    //
    // Point the parent list at the condensed graph parent.
    aList[a] = graphVertices[aList[a]];
  }

  return graphVertices.size();
}


int CountEndParents(vector<int> vertices) {
  


}

int CountTotalBases(vector<FASTASequence> &reads) {
  int numBases = 0;
  int i;
  for (i = 0; i < reads.size(); i++) {
    numBases += reads[i].length;
  }
  return numBases;
}

void InitializeVerticesToIdentity(vector<int> &vertices, int maxVertexIndex = 0) {
  int i;
  int verticesEnd = vertices.size();
  if (maxVertexIndex > 0) {
    verticesEnd = maxVertexIndex;
  }
  for (i = 0; i < verticesEnd; i++) {
    vertices[i] = i;
  }
  for (i = verticesEnd; i < vertices.size(); i++) {
    vertices[i] = 0;
  }
}


int GetParentIndex(int v, vector<int> &vertices) {
  int first = v;
  assert(v < vertices.size());
  
  while(vertices[v] != v) {
    assert(v < vertices.size());
    v = vertices[v];
  }
  int parent = v;
  //
  // Promote this path up to the parent.
  //
  v = first;
  while (vertices[v] != v) {
    int next = vertices[v];
    vertices[v] = parent;
    v = next;
  }

  return parent;
}

void JoinVertices(int a, int b, VertexList &vertices, int &curVertexEnd) {
  int aParent = GetParentIndex(a, vertices);
  int bParent = GetParentIndex(b, vertices);
  
  vertices[aParent] = curVertexEnd;
  vertices[bParent] = curVertexEnd;
  vertices[curVertexEnd] = curVertexEnd;
  ++curVertexEnd;
}

void BuildReadToVertexStartTable(vector<FASTASequence> &reads,
                                 vector<int> &vertexStart) {

  int i;
  int v = 0;
  vertexStart.resize(reads.size());
  for (i = 0; i < reads.size(); i++) {
    vertexStart[i] = v;
    v += reads[i].length;
  }
}

bool LookupReadIndex(string readName, map<string, int> &readMap, int &index) {
  map<string, int>::iterator it = readMap.find(readName);
  if (it == readMap.end()) {
    return false;
  }
  else {
    index = it->second;
    return true;
  }
}

int JoinVerticesByAlignment(vector<AlignmentCandidate<> > &alignments, 
                            map<string, int> &readNameToIndex, 
                            vector<FASTASequence> &reads,
                            vector<int> &readToVertexStart, 
                            VertexList &vertices,
                            int &curVerticesEnd) {
  int b, bi;
  int a;
  int qPos, tPos;
  int qReadIndex, tReadIndex;
  int qLength, tLength;
  int qVertexStart, tVertexStart;

  for (a = 0; a < alignments.size(); a++) {
    if (LookupReadIndex(alignments[a].qName, readNameToIndex, qReadIndex) == 0) {
      cout << "ERROR. Could not find read " << alignments[a].qName << " in the input set of reads." << endl;
      exit(1);
    }
    if (LookupReadIndex(alignments[a].tName, readNameToIndex, tReadIndex) == 0) {
      cout << "ERROR. Could not find read " << alignments[a].tName << " in the input set of reads." << endl;
      exit(1);
    }
    //
    // Read some values off of arrays (for less typing)
    //
    qVertexStart = readToVertexStart[qReadIndex];
    tVertexStart = readToVertexStart[tReadIndex];
    qLength = reads[qReadIndex].length;
    tLength = reads[tReadIndex].length;
    int qStrand = alignments[a].qStrand;
    int tStrand = alignments[a].tStrand;
    /*    cout << "joining " << alignments[a].qName << " (" << qReadIndex << ") with " << alignments[a].tName << " ("
          << tReadIndex << ")" << endl;*/
    
    if (alignments[a].blocks.size() > 0) {
      qPos = alignments[a].qAlignedSeqPos + alignments[a].qPos + alignments[a].blocks[0].qPos;
      tPos = alignments[a].tAlignedSeqPos + alignments[a].tPos + alignments[a].blocks[0].tPos;
      int qVertex = qVertexStart + qPos;
      int tVertex = tVertexStart + tPos;
      //
      // If the vertices are not already joined, create a common
      // parent that links them.
      //
      if (qVertex != tVertex) {
        JoinVertices(qVertex, tVertex, vertices, curVerticesEnd);
      }
      int lastBlock = alignments[a].blocks.size() - 1;
      qPos = alignments[a].qAlignedSeqPos + alignments[a].qPos + alignments[a].blocks[lastBlock].qPos;
      tPos = alignments[a].tAlignedSeqPos + alignments[a].tPos + alignments[a].blocks[lastBlock].tPos;

      qVertex = qVertexStart + qPos;
      tVertex = tVertexStart + tPos;
      //
      // If the vertices are not already joined, create a common
      // parent that links them.
      //
      if (qVertex != tVertex) {
        JoinVertices(qVertex, tVertex, vertices, curVerticesEnd);
      }
    }      
  }
}


void Flatten(vector<int> &vertices) {
  int i;
  int r;
  int v = 0;
  /*  for (r = 0; r < reads.size(); r++) {
    int vr;
    set<int> parents;
    for (vr = v; i = 0; i < reads[r].length; i++; vr++) {
      parents[vr] = i;
    }
  }
  */
  for (i = 0; i < vertices.size(); i++) {
    //
    // Two steps. First, find the parent.
    //
    int p = GetParentIndex(i, vertices);
    vertices[i] = p;
  }
}

int main(int argc, char* argv[]) {
  
  CommandLineParser clp;
  string readsFileName;
  string alignmentsFileName;
  string outputFileName;
  clp.RegisterStringOption("reads", &readsFileName, "Reads used for alignments.");
  clp.RegisterStringOption("alignments", &alignmentsFileName, "SAM formatted alignments.");
  clp.RegisterStringOption("outfile", &outputFileName, "Alignment output.");
  clp.RegisterPreviousFlagsAsHidden();
  clp.RegisterFlagOption("v", &verbose, "");
  clp.ParseCommandLine(argc, argv);

  
  vector<FASTASequence> reads;

  FASTAReader fastaReader;
  fastaReader.Initialize(readsFileName);
  fastaReader.ReadAllSequences(reads);
  map<string, int> readNameToIndex;
  vector<int> readToVertexStart;
  BuildReadNameToIndexMap(reads, readNameToIndex);
  BuildReadToVertexStartTable(reads, readToVertexStart);
  int numBases = CountTotalBases(reads);
  
  VertexList vertices;
  vertices.resize(numBases*2-1);
  InitializeVerticesToIdentity(vertices, numBases);

  SAMReader<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> samReader;
  samReader.Initialize(alignmentsFileName);
  AlignmentSet<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> alignmentSet;
  samReader.ReadHeader(alignmentSet);
  
  SAMAlignment samAlignment;
  AlignmentCandidate<> alignment;
  int curVerticesEnd = numBases;
  int numAlignedBases = 0;
  int alignmentIndex = 0;
  while ( samReader.GetNextAlignment( samAlignment ) ) {
    vector<AlignmentCandidate<> > alignments;
    SAMAlignmentsToCandidates(samAlignment,
                              reads,
                              readNameToIndex,
                              alignments);

    //
    // Glue vertices together from the graph.
    //
    //    cout << "joining " << alignments[0].qName << " " << alignments[0].tName << endl;
    
    JoinVerticesByAlignment(alignments, 
                            readNameToIndex, 
                            reads,
                            readToVertexStart, 
                            vertices,
                            curVerticesEnd);
    int i;
    ++alignmentIndex;
    if (alignmentIndex % 1000 == 0) {
      cout << alignmentIndex << " " << curVerticesEnd - numBases << endl;
    }
  }

  Flatten(vertices);
  set<int> newVertices;
  int i;
  for (i = numBases; i < numBases*2-1; i++) {
    newVertices.insert(vertices[i]);
  }
  cout << "There are " << newVertices.size() << " parent vertices." << endl;
  //
  // Will never grow vertices past this.
  //

  /*  
  
  cout << "index size " << endl;
  int i;
  map<int,int> hist;
  for (i = 0; i < vertices.size(); i++) {
    int count = vertices[i].out.size();
    if (hist.find(count) == hist.end()) {
      hist[count] = 1;
    }
    else {
      hist[count] = hist[count] + 1;
    }
  }
  map<int,int>::iterator mapIt;
  cout << "count freq" << endl;
  for (mapIt = hist.begin(); mapIt != hist.end(); ++mapIt) {
    cout << mapIt->first << " " << mapIt->second << endl;
  }
  */
}
  
  
  
