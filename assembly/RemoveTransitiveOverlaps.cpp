#include "RepeatGraph.h"
#include "RGEdgeOperations.h"
#include <string>
#include <set>
#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;

int DFSForTargetVertices(RepeatGraph<string> &rg,
												 VectorIndex curVertex,
//												 set<VectorIndex> &targetVertices,
												 VectorIndex &targetVertex,
												 set<VectorIndex> &searchVertexList,
												 set<VectorIndex> &searchEdgeList,
//												 set<VectorIndex> &tEdges,
												 int depth) {
  if (depth == 0) {
		return 0;
	}

	searchVertexList.insert(curVertex);
	VectorIndex outEdgeIndex;
	for (outEdgeIndex = 0; outEdgeIndex < rg.vertices[curVertex].outEdges.size(); outEdgeIndex++) {
		VectorIndex edgeIndex = rg.vertices[curVertex].outEdges[outEdgeIndex];
		VectorIndex destVertex = rg.edges[edgeIndex].dest;
//		cout << depth << " edge " << edgeIndex << " from " << curVertex << " goes to " << destVertex << " target?= " << targetVertex << endl;
		
		if (searchEdgeList.find(edgeIndex) != searchEdgeList.end() or 
				searchVertexList.find(destVertex) != searchVertexList.end()) {
			//
			// Already searched this vertex.
			//
			continue;
		}


//		if (targetVertices.find(destVertex) != targetVertices.end()) {
			if (targetVertex == destVertex) {
//			tEdges.insert(edgeIndex);
			return 1;
		}


		searchEdgeList.insert(edgeIndex);
		// Search forward on active edges.

		if (DFSForTargetVertices(rg, destVertex, targetVertex, searchVertexList, searchEdgeList, depth-1)) {
			return 1;
		}
	}
	return 0;
}


int main(int argc, char* argv[]) {
	string rgInName, rgOutName;
	int radius;

	if (argc < 4) {
		cout << "usage: removeTransitiveOverlaps in.rg radius out.rg" << endl;
		exit(1);
	}

	rgInName   = argv[1];
	radius = atoi(argv[2]);
	rgOutName  = argv[3];

	ofstream rgOut;
	CrucialOpen(rgOutName,rgOut, std::ios::out);
	RepeatGraph<string> rg;
	
	rg.ReadGraph(rgInName);

	VectorIndex vertexIndex;
	VectorIndex outEdgeIndex;

	//
	// Remove transitive overlaps from each vertex.
	//
	VectorIndex edgeIndex;
	for (edgeIndex =0 ; edgeIndex < rg.edges.size(); edgeIndex++) {
		rg.edges[edgeIndex].connected = true;
	}
	set<std::pair<VectorIndex, VectorIndex> > srcDestToRemove;

	for (vertexIndex = 0; vertexIndex < rg.vertices.size(); vertexIndex++) {
		// Ensure that there is at least one edge
//    cout << "removing transitive overlaps for " << vertexIndex << " " << rg.vertices[vertexIndex].key << endl;
		VectorIndex outEdgeIndex, outEdge;
		VectorIndex destVertex;

		//
		// Make a quickly query-able list of adjacent vertices.
		//
/*
		cout << endl;
		for (outEdgeIndex = 0; outEdgeIndex < rg.vertices[vertexIndex].outEdges.size(); outEdgeIndex++) {
			outEdge    = rg.vertices[vertexIndex].outEdges[outEdgeIndex];
			destVertex = rg.edges[outEdge].dest;
			cout << destVertex << " is an adjacency " << endl;
			adjacentVertices.insert(destVertex);
		}
*/
		int directPathExists = 0;

		//
		// Search starting at all adjacent vertices for other adjacent
		// vertices, these are transitively linked with a path of length 1.
		//
		if (rg.vertices[vertexIndex].outEdges.size() > 1) {
			for (outEdgeIndex = 0; outEdgeIndex < rg.vertices[vertexIndex].outEdges.size(); outEdgeIndex++) {
				outEdge    = rg.vertices[vertexIndex].outEdges[outEdgeIndex];
				destVertex = rg.edges[outEdge].dest;
				set<VectorIndex> searchVertexSet;
				set<VectorIndex> searchEdgeSet;
				searchVertexSet.insert(vertexIndex);
				searchEdgeSet.insert(outEdge);
				set<VectorIndex> adjacentVertices;
				adjacentVertices.insert(rg.edges[outEdge].dest);
				directPathExists = DFSForTargetVertices(rg, vertexIndex, 
																								//adjacentVertices, //
																								//targets
																								destVertex,
																								searchVertexSet, searchEdgeSet, radius);
				if (directPathExists ) {
					cout << "edge " << outEdge << " -> " << destVertex << " links transitively linked vertices." << endl;
					//
					// Insert edges to remove as adjacencies so that the
					// forward/reverse overlaps are handled without having to
					// look too many edges up.
					//
					srcDestToRemove.insert(pair<VectorIndex, VectorIndex>(rg.edges[outEdge].src, rg.edges[outEdge].dest));
					srcDestToRemove.insert(pair<VectorIndex, VectorIndex>(2*(rg.edges[outEdge].dest/2) + !(rg.edges[outEdge].dest%2),
																		2*(rg.edges[outEdge].src/2)  + !(rg.edges[outEdge].src%2)));
				}
			}
		}
	}

  if (srcDestToRemove.size() == 0) {
    cout << "Removing transitive edges is not necessary as none were found." << endl;
    return 0;
  }

	MarkEdgePairsForRemoval(srcDestToRemove, rg.vertices, rg.edges);
	
	//
	// Remove the edges that link transitively linked vertices.
	//
  VectorIndex lastEdgeIndex = rg.edges.size();
	
	RemoveUnconnectedEdges(rg.vertices, rg.edges);

	rg.WriteGraph(rgOut);
	return 0;
}
