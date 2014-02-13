#include "RepeatGraph.h"
#include "files/AlignmentReader.h"
#include "FASTASequence.h"
#include "FASTAReader.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "datastructures/metagenome/FASTATitleDictionary.h"
#include "utils.h"
#include "Path.h"
#include "PathCollection.h"

using namespace std;
#include <string>


int main(int argc, char* argv[]) {

	if (argc < 6) {
		cout << "usage: rm4ToPathList rm4File graphFile vertexSequences pathFileName minPathLength" << endl;
		exit(1);
	}
	
	string rm4FileName;
	string vertexSequenceFileName;
	string repeatGraphName;
	string pathFileName;
	int    minPathLength;

	rm4FileName            = argv[1];
	repeatGraphName        = argv[2];
	vertexSequenceFileName = argv[3];
	pathFileName           = argv[4];
	minPathLength          = atoi(argv[5]);

	FASTAReader vertexSequenceReader;
	AlignmentReader<AlignmentCandidate<FASTASequence, FASTASequence> > alignmentReader;

	//
	// Initialize file I/O
	//
	vertexSequenceReader.Init(vertexSequenceFileName);
  alignmentReader.Init(rm4FileName);
	RepeatGraph<string> rg;
	rg.ReadGraph(repeatGraphName);

	ofstream pathsOut;
	CrucialOpen(pathFileName, pathsOut, std::ios::out);

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

	VectorIndex  ai, aiEnd;
  PathCollection<UInt> paths;
	for (ai = 0; ai < alignments.size();) {
		// 
		// Advance aiEnd to point at next aligned read.
		//
		for (aiEnd = ai+1; aiEnd < alignments.size(); aiEnd++ ){
			if (alignments[aiEnd].qName != alignments[ai].qName) {
				break;
			}
		}

		if (aiEnd - ai >= minPathLength) {
			VectorIndex pi;
			VectorIndex vertexIndex;
			Path<UInt> path;
			Path<UInt> pathRC;
			//pathsOut << ai << " " << alignments[ai].qName << " " << alignments[ai].qLength;
			// advance by two here since the
			pi = ai;
			if (titleDictionary.LookupSequence(alignments[pi].tName, (int&) vertexIndex) == 0) {
				cout << "ERROR, vertex name missing from sequence file: " << alignments[pi].tName << endl;
				exit(1);
			}
//			pathsOut << " " << vertexIndex << " " << alignments[pi].sumQVScore;

			path.AddEdge(vertexIndex*2 + alignments[pi].tStrand);
			for (pi++; pi < aiEnd; pi+=2 ){ 
				if (titleDictionary.LookupSequence(alignments[pi].tName, (int&) vertexIndex) == 0) {
					cout << "ERROR, vertex name missing from sequence file: " << alignments[pi].tName << endl;
					exit(1);
				}
				path.AddEdge(vertexIndex*2 + alignments[pi].tStrand);
			}
			
			for (pi = 0; pi < path.size(); pi++ ){
		  	pathRC.AddEdge(ReverseEdgeIndex(path[path.size() - pi - 1]));
      }

			//
			// Check this path out to see if it is anything interesting.
			//
			
		  if (path.size() > 1) {
				paths.AddPath(path);
				paths.AddPath(pathRC);
      }	
		}
		ai = aiEnd;
	}
  paths.Condense();
	paths.RemoveShortPaths(2);
  paths.Print(pathsOut);
	return 0;
}
