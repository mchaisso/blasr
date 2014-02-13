#ifndef GRAPH_FLOW_H_
#define GRAPH_FLOW_H_


/* 
	 Define some hacky operations for determining flow through a repeat graph.
*/

#include "RepeatGraph.h"
#include "PathOperations.h"
#include "Types.h"

void AssignMinimumFlowToEdges(RepeatGraph<string> &rg, int minPathLength) {

	//
	// First look for a path that is "long".
	//

	VectorIndex edgeIndex;
	VectorIndex vertexIndex;

	for (vertexIndex = 0; vertexIndex < rg.vertices.size(); vertexIndex++ ) {
		if (rg.vertices[vertexIndex].InDegree() != 1 or rg.vertices[vertexIndex].OutDegree() != 1) {
			// 
			// This is a branching vertex, look to see if the path length is long enough
			// to consider a 'real' sequence.
			//
			VectorIndex destVertex;
			VectorIndex outEdgeIndex;
			int pathLength = 0;
			//			cout << " checking vertex: " << rg.vertices[vertexIndex].key << endl;
			for (outEdgeIndex = 0; outEdgeIndex < rg.vertices[vertexIndex].OutDegree(); outEdgeIndex++) {
				vector<VectorIndex> pathEdgeList;
				destVertex = rg.edges[rg.vertices[vertexIndex].outEdges[outEdgeIndex]].dest;
				pathEdgeList.push_back(rg.vertices[vertexIndex].outEdges[outEdgeIndex]);
				while (rg.vertices[destVertex].InDegree() == 1 and rg.vertices[destVertex].OutDegree() == 1) {
					VectorIndex pathEdge = rg.vertices[destVertex].outEdges[0]; 
					destVertex = rg.edges[pathEdge].dest;
					pathEdgeList.push_back(pathEdge);
					++pathLength;
				}

				if (pathEdgeList.size() > minPathLength) {
					//
					// This path is to be trusted.
					//
					//				cout << "path of length " << pathEdgeList.size() << " is to be trusted." << endl;
					UInt pathIndex;
					for (pathIndex = 0; pathIndex < pathEdgeList.size(); pathIndex++ ) {
						/*		cout << "edge (" << rg.vertices[rg.edges[pathEdgeList[pathIndex]].src].key << ", " 
									<< rg.vertices[rg.edges[pathEdgeList[pathIndex]].dest].key << ") requires one flow." << endl;*/
						rg.edges[pathEdgeList[pathIndex]].minFlow = 1;
						
						//
						// This is a source vertex, so it has flow in of 1 from a
						// virtual super-source.
						//
						if (rg.vertices[rg.edges[pathEdgeList[pathIndex]].src].InDegree() == 0) {
							rg.vertices[rg.edges[pathEdgeList[pathIndex]].src].flowIn = 1;
						}
						if (rg.vertices[rg.edges[pathEdgeList[pathIndex]].dest].OutDegree() == 0) {
							rg.vertices[rg.edges[pathEdgeList[pathIndex]].dest].flowOut = 1;
						}
					}
				}
			}		}
	}
}


void AssignVertexFlowBalance(RepeatGraph<string> &rg) {

	VectorIndex vertexIndex;

	for (vertexIndex = 0; vertexIndex < rg.vertices.size(); vertexIndex++ ) {
		if (rg.vertices[vertexIndex].InDegree() == 0) {
			// This vertex emits flow.
			rg.vertices[vertexIndex].flowBalance = 1;
		}
		else if (rg.vertices[vertexIndex].OutDegree() == 0) {
			// This vertex can absorb some flow.
			rg.vertices[vertexIndex].flowBalance = -1;
		}
		else {
			// require that flow is balanced at this vertex.
			rg.vertices[vertexIndex].flowBalance = 0;
		}
	}
}

int GetNodeFlowRequirement(RepeatGraph<string> &rg, UInt vertexIndex) {
	UInt outEdgeIndex, outEdge;
	UInt netFlow = 0;
	UInt netMinFlow = 0;

	for (outEdgeIndex = 0; outEdgeIndex < rg.vertices[vertexIndex].outEdges.size(); outEdgeIndex++) {
		outEdge = rg.vertices[vertexIndex].outEdges[outEdgeIndex];
		netFlow += rg.edges[outEdge].flow;
		netMinFlow += rg.edges[outEdge].minFlow;
	}
	return netMinFlow - netFlow;
}

bool GetOpenEdge(RepeatGraph<string> &rg, UInt curVertex, UInt &openEdge) {
	UInt outEdgeIndex, outEdge;
	for (outEdgeIndex = 0; outEdgeIndex < rg.vertices[curVertex].outEdges.size(); outEdgeIndex++) {
		outEdge = rg.vertices[curVertex].outEdges[outEdgeIndex];
		if (rg.edges[outEdge].flow < rg.edges[outEdge].minFlow) {
			openEdge = outEdge;
			return true;
		}
	}
	return false;
}
	

bool DFSForCapableVertexAndEdge(RepeatGraph<string> &rg, UInt curVertex, int maxDepth, UInt &capVertex, UInt &capEdge, vector<UInt> &capPath) {
	UInt outEdge, outEdgeIndex;
	if (maxDepth == 0) {
		return false;
	}
	for (outEdgeIndex = 0; outEdgeIndex < rg.vertices[curVertex].outEdges.size(); outEdgeIndex++) {
		UInt curCapPathSize = capPath.size();
		outEdge = rg.vertices[curVertex].outEdges[outEdgeIndex];
		UInt destVertex;
		destVertex = rg.edges[outEdge].dest;
		UInt edgeToDest;
		AdvanceAndStoreSimplePath(rg, curVertex, outEdge, destVertex, edgeToDest, capPath);
		cout << maxDepth << " path of len " << capPath.size() << endl;
		cout << rg.vertices[destVertex].InDegree() << " " << rg.vertices[destVertex].OutDegree() << endl;
		if (GetNodeFlowRequirement(rg, destVertex) > 0) {
			// This node has excess flow capacity going out of it.
			// look up which edge it is.
			cout << "found excess capacity at " << rg.vertices[destVertex].key << endl;
			assert(GetOpenEdge(rg, destVertex, capEdge) == true);
			UInt capDestVertex, capDestEdge;
			// 
			// Store all the way to the vertex that the flow will be pushed
			// to.
			//
			AdvanceAndStoreSimplePath(rg, curVertex, capEdge, capDestVertex, capDestEdge, capPath);
			cout << "capacity pushed to " << rg.vertices[capDestVertex].key << endl;
			capEdge = capDestEdge;
			capVertex = capDestVertex;
			return true;
		}
		else if (DFSForCapableVertexAndEdge(rg, destVertex, maxDepth-1, capVertex, capEdge, capPath)) {
			// early exit.
			return true;
		}
		else {
			// No capable vertex on this path. Get rid of the path that was stored.
			capPath.resize(curCapPathSize);
		}			
	}
	return false;
}

void BalanceKirchhoffFlow(RepeatGraph<string> &rg) {
	int maxIter = 100;
	//
	// First assign excess flow to vertices that are start nodes.
	//
	bool flowIsBalanced = false;
	int iter = 0;
	while (flowIsBalanced == false and iter < maxIter) {
		cout << "iter: " << iter << endl;
		++iter;
		//
		// Check for flow-unbalanced vertices.
		//
		UInt vertexIndex;
		flowIsBalanced = true;
		for (vertexIndex = 0; vertexIndex < rg.vertices.size(); ++vertexIndex) {
			int flowExcess = rg.vertices[vertexIndex].flowIn - rg.vertices[vertexIndex].flowOut;
			//			cout << "fe: " << flowExcess << endl;
			cout << "VERTEX: " << rg.vertices[vertexIndex].key << endl;
			if (flowExcess > 0) {
				cout << "vertex: " << vertexIndex << " has excess flow: " << flowExcess << endl;
				cout << " in:" << rg.vertices[vertexIndex].flowIn 
						 << " out: " << rg.vertices[vertexIndex].flowOut
						 << " degrees, in: " << rg.vertices[vertexIndex].InDegree()
						 << " out: " << rg.vertices[vertexIndex].OutDegree() << " " << flowExcess << endl;
				flowIsBalanced = false;
				//
				// This vertex has excess flow to get rid of.
				//
				UInt edgeIndex;
				// Search through out edges for one that has not yet met its
				// flow capacity.
				UInt e;
				// cout << "vertex has " <<
				// rg.vertices[vertexIndex].outEdges.size() << " out edges. "
				// << endl;
				bool foundFlowCapableEdge = false;
				for (e = 0; e < rg.vertices[vertexIndex].outEdges.size(); e++ ) {
					edgeIndex = rg.vertices[vertexIndex].outEdges[e];
					cout << "edge flow of " << edgeIndex << " " << rg.edges[edgeIndex].flow << " " << rg.edges[edgeIndex].minFlow << endl;
					if (rg.edges[edgeIndex].flow < rg.edges[edgeIndex].minFlow) {
						//
						// Send some flow out this edge.
						//
						foundFlowCapableEdge = true;
						cout << "sending flow to " << rg.vertices[rg.edges[edgeIndex].dest].key << endl;
						rg.vertices[vertexIndex].flowOut++;
						cout << rg.vertices[vertexIndex].flowIn << " " << rg.vertices[vertexIndex].flowOut << endl;
						
						UInt lastEdge, lastVertex;
						vector<UInt> pathEdges;
						AdvanceAndStoreSimplePath(rg, vertexIndex, edgeIndex, lastVertex, lastEdge, pathEdges);
						UInt pathIndex;
						//
						// Update the flow along this path.
						//
						for (pathIndex = 0; pathIndex < pathEdges.size(); pathIndex++ ) {
							rg.edges[pathEdges[pathIndex]].flow++;
						}
						rg.vertices[lastVertex].flowIn++;
						cout << "incremented flow of " << rg.vertices[lastVertex].key << " to " << rg.vertices[lastVertex].flowIn << endl;
						cout << rg.vertices[lastVertex].InDegree() << " " << rg.vertices[lastVertex].OutDegree() << endl;
						break;
					}
				}
				if (foundFlowCapableEdge == false) {
					UInt capVertex, capEdge;
					vector<UInt> capPath;
					cout << "no flow capable edge found." << endl;
					if (DFSForCapableVertexAndEdge(rg, vertexIndex, 5, capVertex, capEdge, capPath)) {
						cout << "HOORRAAY! found a capable path of length " << capPath.size() << endl;
						UInt pathIndex;
						for (pathIndex = 0; pathIndex < capPath.size(); pathIndex++) {
							//
							// Increment the flow along this path.
							//
							rg.edges[capPath[pathIndex]].flow++;
						}
						// the last edge should go to a non-branching vertex.
						assert(rg.vertices[capVertex].InDegree() != 0 or
									 rg.vertices[capVertex].OutDegree() != 0);
						//
						// Now, flow was pushed out on the path to capEdge, now
						// push that along to the next branching vertex.
						//
						rg.vertices[vertexIndex].flowOut++;
						rg.vertices[capVertex].flowIn++; 
					}
				}
			}
		}
	}
}



#endif
