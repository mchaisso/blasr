#include "RepeatGraph.h"
#include "RGEdgeOperations.h"
#include "FASTAReader.h"
#include "FASTASequence.h"
#include <string>
#include <set>
#include <iostream>
using namespace std;



int main(int argc, char* argv[]) {
	string rgInName, rgOutName;
	int minPathLength;
	string vertexSequenceFileName;
	if (argc < 5) {
		cout << "usage: trimShortEnds in.rg  vertexSequences minPathLength out.rg" << endl;
		exit(1);
	}

	rgInName      = argv[1];
	vertexSequenceFileName = argv[2];
	minPathLength = atoi(argv[3]);
	rgOutName     = argv[4];

	ofstream rgOut;
	CrucialOpen(rgOutName, rgOut, std::ios::out);
	FASTAReader vertexSequenceReader;
	vertexSequenceReader.Init(vertexSequenceFileName);

	RepeatGraph<string> rg;
	vector<FASTASequence> vertexSequences;
	rg.ReadGraph(rgInName);
	vertexSequenceReader.ReadAllSequences(vertexSequences);

	VectorIndex vertexIndex;
	VectorIndex outEdgeIndex;
	VectorIndex edgeIndex;
	
	if (rg.edges.size() == 0) {
		cout << "LIKELY INVALID GRAPH. There are no edges." << endl;
		return 0;
	}
	//
	// At first, any edge that exists is connected to a vertex. This
	// will change as low coverage edges are deleted and replaced by
	// high coverage edges from the end of the array.
	//
	for (edgeIndex = 0; edgeIndex < rg.edges.size(); edgeIndex++) {
		rg.edges[edgeIndex].connected = true;
	}
	set<std::pair<VectorIndex, VectorIndex> > srcDestToRemove;
	
	for (vertexIndex = 0; vertexIndex < rg.vertices.size(); vertexIndex++) {
		if (rg.vertices[vertexIndex].inEdges.size() == 0 and
				rg.vertices[vertexIndex].outEdges.size() == 1) {
			//
			// This is a source.  Traverse this until a branching vertex or the end is found.
			//
			vector<VectorIndex> path;
			path.push_back(vertexIndex);
			int pathLength = 0;
			VectorIndex pathVertex;
			VectorIndex pathEdge;
			pathEdge = rg.vertices[vertexIndex].outEdges[0];
			pathVertex = rg.edges[pathEdge].dest;
			while (rg.vertices[pathVertex].inEdges.size() == 1 and
						 rg.vertices[pathVertex].outEdges.size() == 1) {
				path.push_back(pathVertex);
				pathEdge   =  rg.vertices[pathVertex].outEdges[0];
				pathVertex =  rg.edges[pathEdge].dest;
				pathLength += vertexSequences[pathVertex/2].length;
			}
			pathLength += vertexSequences[pathVertex/2].length;
			path.push_back(pathVertex);
			if (pathLength < minPathLength and path.size() < 3) {
				//
				// Remove this path, it is too short.
				// Also remove the complement.
				//
				cout << "trimming path of " << path.size() << " is of sequence length " << pathLength << endl;

				VectorIndex pathIndex;
				for (pathIndex = 0; pathIndex < path.size() - 1; pathIndex++) {
					srcDestToRemove.insert(pair<VectorIndex, VectorIndex>(path[pathIndex], path[pathIndex+1]));
					srcDestToRemove.insert(pair<VectorIndex, VectorIndex>(2*(path[pathIndex+1]/2) + !(path[pathIndex+1]%2),
																																2*(path[pathIndex]/2) + !(path[pathIndex]%2)));
				}
			}
		}
	}

	MarkEdgePairsForRemoval(srcDestToRemove, rg.vertices, rg.edges);
	RemoveUnconnectedEdges(rg.vertices, rg.edges);

	rg.WriteGraph(rgOut);
	return 0;
}
