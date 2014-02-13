#include "RepeatGraph.h"
#include "RGEdgeOperations.h"
#include <string>
#include <set>
#include <iostream>
using namespace std;



int main(int argc, char* argv[]) {
	string rgInName, rgOutName;
	int minCoverage;

	if (argc < 4) {
		cout << "usage: removeTransitiveOverlaps in.rg minCoverage out.rg" << endl;
		exit(1);
	}

	rgInName    = argv[1];
	minCoverage = atoi(argv[2]);
	rgOutName   = argv[3];

	ofstream rgOut;
	CrucialOpen(rgOutName, rgOut, std::ios::out);

	RepeatGraph<string> rg;
	
	rg.ReadGraph(rgInName);

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
	
	VectorIndex numEdges = rg.edges.size();
	for (edgeIndex = 0; edgeIndex < numEdges; ) {
		if (rg.edges[edgeIndex].count >= minCoverage) {
			// This edge is fine.
			edgeIndex++;
		}
		else {
			// This edge needs to be deleted. Find the first edge at the end that
			// will not be deleted anyway, and move it here.

			while (numEdges > edgeIndex and rg.edges[numEdges-1].count < minCoverage) {
				UnlinkDirected(rg.vertices, rg.edges, numEdges-1);
				numEdges--;
			}

			//
			// If exhausted all edges, just break since all are deleted.
			//
			if (numEdges == edgeIndex) {
				continue;
			}
			VectorIndex src  = rg.edges[edgeIndex].src;
			VectorIndex dest = rg.edges[edgeIndex].dest;
			//
			// Get rid of this low coverage edge.
			UnlinkDirected(rg.vertices, rg.edges, edgeIndex);

			//
			// Pack in one from a higher coverage.
			rg.edges[edgeIndex] = rg.edges[numEdges-1];

			//
			// Update the connecting vertex.
			UpdateOutEdge(rg.vertices, rg.edges, numEdges-1, edgeIndex);
			UpdateInEdge(rg.vertices, rg.edges, numEdges-1, edgeIndex);
			--numEdges;

			++edgeIndex;
		}
	}
	rg.edges.resize(numEdges);
	rg.WriteGraph(rgOut);
	return 0;
}
