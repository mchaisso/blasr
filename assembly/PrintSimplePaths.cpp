#include "RepeatGraph.h"
#include "FASTAReader.h"
#include "FASTASequence.h"
#include <string>
#include <vector>

using namespace std;

int main(int argc, char* argv[]) {
	
	string rgFileName, vertexSeqFileName, pathFileName;

	ofstream pathFileOut
	CrucialOpen(pathFileName,pathFileOut, std::ios::out);

	FASTAReader vertexSequenceReader;
	vertexSequenceReader.Init(vertexSequenceFileName);

	//
	// Input necessary data
	//
	vector<FASTASequence> vertexSequences;
	cout << "Reading sequences." << endl;
	vertexSequenceReader.ReadAllSequences(vertexSequences);

	vector<FASTASequence> vertexRCSequences;
	VectorIndex vertexIndex;	
	for (vertexIndex = 0; vertexIndex < vertexSequences.size(); vertexIndex++ ){
		vertexSequences.MakeRC(vertexRCSequences[vertexIndex]);
	}
	
	
	RepeatGraph<string> rg;
	rg.ReadGraph(rgInName);


	VectorIndex outEdgeIndex;

	for (vertexIndex = 0; vertexIndex < vertices.size(); vertexIndex++ ){
		if ((rg.vertices[vertexIndex].inEdges.size() > 1 or
				 rg.vertices[vertexIndex].inEdges.size() == 0) and
				rg.vertices[vertexIndex].outEdges.size()) {
			//
			// This is a branching vertex.
			// 
			
			

		}



	}



}
