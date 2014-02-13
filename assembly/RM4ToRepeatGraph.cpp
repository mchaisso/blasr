#include "RepeatGraph.h"
#include "files/AlignmentReader.h"
#include "FASTASequence.h"
#include "FASTAReader.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "datastructures/metagenome/FASTATitleDictionary.h"
#include "utils.h"

using namespace std;
#include <string>


int main(int argc, char* argv[]) {

	if (argc < 4) {
		cout << "usage: rm4ToRepeatGraph rm4File vertexSequences repeatGraphName" << endl;
		exit(1);
	}

	string rm4FileName;
	string vertexSequenceFileName;
	string repeatGraphName;

	rm4FileName            = argv[1];
	vertexSequenceFileName = argv[2];
	repeatGraphName        = argv[3];

	FASTAReader vertexSequenceReader;
	AlignmentReader<AlignmentCandidate<FASTASequence, FASTASequence> > alignmentReader;


	//
	// Initialize file I/O
	//
	vertexSequenceReader.Init(vertexSequenceFileName);
  alignmentReader.Init(rm4FileName);
	ofstream graphOut;
  CrucialOpen(repeatGraphName, graphOut, std::ios::out);
	
	//
	// Do actual I/O. 
	//

	vector<FASTASequence> vertexSequences;
	vector<AlignmentCandidate<FASTASequence, FASTASequence> > alignments;

	cout << "Reading sequences." << endl;
	vertexSequenceReader.ReadAllSequences(vertexSequences);
	cout << "Reading alignments." <<endl;
	alignmentReader.ReadAllAlignments(alignments);

	cout << "got " << alignments.size() << " alignments." << endl;
	//
	// Create a map from vertex name to index in the graph.
	//
	FASTATitleDictionary titleDictionary;
	titleDictionary.AddAllSequences(vertexSequences);

	RepeatGraph<string> rg;

	rg.vertices.resize(vertexSequences.size()*2);
	
	//
	// Assign the key values so that the vertices may be tracked as the graph changes.
	//
	VectorIndex vertexIndex;
	for (vertexIndex = 0; vertexIndex < vertexSequences.size(); vertexIndex++) {
		string vertexName = vertexSequences[vertexIndex].GetName();
		rg.vertices[vertexIndex*2].SetKey(vertexName);
		string rcName = vertexName+"rc";
		rg.vertices[vertexIndex*2+1].SetKey(rcName);
	}

	
	VectorIndex alnIndex;
	for (alnIndex = 0; alnIndex < alignments.size(); alnIndex+= 2) {
		//		cout << "Connecting " << alignments[alnIndex].tName << " -> " << alignments[alnIndex+1].tName << endl;
		int srcIndex;
		int destIndex;
		if (titleDictionary.LookupSequence(alignments[alnIndex].tName, srcIndex) == 0) {
			cout << "ERROR! Sequence " << alignments[alnIndex].tName
					 << " is in an alignment but does not correspond to a vertex." << endl;
			exit(1);
		}
		if (titleDictionary.LookupSequence(alignments[alnIndex+1].tName, destIndex) == 0) {
			cout << "ERROR! Sequence " << alignments[alnIndex+1].tName
					 << " is in an alignment but does not correspond to a vertex." << endl;
			exit(1);
		}
		//
		// Add the edge to the graph.
		//
		rg.Link(srcIndex*2+alignments[alnIndex].tStrand, destIndex*2+alignments[alnIndex+1].tStrand);
		//
		// Add the reverse complement edge to the graph.
		//
		rg.Link(destIndex*2+ (!alignments[alnIndex+1].tStrand),srcIndex*2+ (!alignments[alnIndex].tStrand));
	}
	/*
	for (vertexIndex = 0; vertexIndex < rg.vertices.size(); vertexIndex++) {
		cout << "vertex " << vertexIndex << " key: " << rg.vertices[vertexIndex].key << " of length " << vertexSequences[vertexIndex/2].length << " has " << rg.vertices[vertexIndex].outEdges.size() << " edges " << endl;
		int outEdgeIndex;
		for (outEdgeIndex = 0; outEdgeIndex < rg.vertices[vertexIndex].outEdges.size(); outEdgeIndex++) {
			cout << rg.edges[rg.vertices[vertexIndex].outEdges[outEdgeIndex]].count << " " << vertexSequences[rg.edges[rg.vertices[vertexIndex].outEdges[outEdgeIndex]].dest/2].GetName() << ", ";
		}
		cout << endl;
	}
	*/
	rg.WriteGraph(graphOut);
	return 0;
	
}
