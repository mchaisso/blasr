#include "RepeatGraph.h"
#include "FASTAReader.h"
#include "FASTASequence.h"
#include "GraphFlow.h"
#include <string>
#include <vector>
#include <sstream>

using namespace std;

void PrintUsage() {
	cout << "usage: printScaffolds rgFileName vertexSeqFileName scaffoldDirName [-separate] [-repeats]" << endl;
}

int main(int argc, char* argv[]) {
	
	if (argc < 4) {
		PrintUsage();
		exit(1);
	}

	string rgFileName, vertexSeqFileName, scaffoldDirName;
	
	rgFileName         = argv[1];
	vertexSeqFileName  = argv[2];
	scaffoldDirName    = argv[3];

	string repeatFileName = "";
	bool printRepeatsSeparately = false;
	int argi = 4;
	bool printSeparate=false;
	while (argi < argc) {
		if (strcmp(argv[argi], "-separate") == 0) {
			printSeparate=true;
		}
		else if (strcmp(argv[argi], "-repeats") == 0) {
			printRepeatsSeparately = true;
			repeatFileName = argv[++argi];
		}
		else {
			cout << "bad option: " << argv[argi] << endl;
			PrintUsage();
			exit(1);
		}
		++argi;
	}
	
	FASTAReader vertexSequenceReader;
	vertexSequenceReader.Init(vertexSeqFileName);

	//
	// Input necessary data
	//
	vector<FASTASequence> vertexSequences;
	vertexSequenceReader.ReadAllSequences(vertexSequences);
	RepeatGraph<string> rg;
	rg.ReadGraph(rgFileName);

	vector<FASTASequence> vertexRCSequences;
	VectorIndex vertexIndex;	
	vertexRCSequences.resize(vertexSequences.size());
	for (vertexIndex = 0; vertexIndex < vertexSequences.size(); vertexIndex++ ){
		vertexSequences[vertexIndex].MakeRC(vertexRCSequences[vertexIndex]);
	}
	
	VectorIndex outEdgeIndex;
	int scaffoldIndex = 0;
	ofstream scaffoldOut;

	if (printSeparate==false) {
		// scaffold dir name is really a file name here.
		CrucialOpen(scaffoldDirName, scaffoldOut, std::ios::out);
	}
	for (vertexIndex = 0; vertexIndex < rg.vertices.size(); vertexIndex++ ){
		rg.vertices[vertexIndex].traversed = false;
	}

	//
	// Set up flow for calling multiplicity.
	//
	/*
		Test all this out later.
	AssignMinimumFlowToEdges(rg, 2);
	AssignVertexFlowBalance(rg);
	BalanceKirchhoffFlow(rg);

	UInt edgeIndex;
	for (edgeIndex = 0; edgeIndex < rg.edges.size(); edgeIndex++) {
		if (rg.edges[edgeIndex].flow > 1) {
			cout << edgeIndex << " " << rg.edges[edgeIndex].flow << endl;
		}
	}
	*/

	int numPrintedVertices = 0;
	for (vertexIndex = 0; vertexIndex < rg.vertices.size(); vertexIndex++ ){
		//
		// Look to see if this vertex is a branching vertex.
		//
		if ((rg.vertices[vertexIndex].inEdges.size() != 1 or
				 rg.vertices[vertexIndex].outEdges.size() != 1) and
				rg.vertices[vertexIndex].traversed == false) {

			//
			// This is a branching vertex.  Print all paths from this vertex, but not the vertex
			// itself if it appears repetitive. 
			//
			VectorIndex outEdgeIndex;
			bool printedThisVertex = false;
			for (outEdgeIndex = 0; outEdgeIndex < rg.vertices[vertexIndex].outEdges.size(); outEdgeIndex++ ){
				//
				// This is a branching vertex.
				// 

				VectorIndex pathIndex;
				stringstream scaffoldFileNameStrm;
				cout << " printing scaffold: " << scaffoldIndex << endl;
				if (printSeparate) {
					scaffoldFileNameStrm << scaffoldDirName << "/" << scaffoldIndex << ".fasta";
					string scaffoldFileName = scaffoldFileNameStrm.str();
					CrucialOpen(scaffoldFileName, scaffoldOut, std::ios::out);
				}			
				++scaffoldIndex;

				//
				// Store the nonbranching path in a list so that it may be quickly processed.
				//
				bool pathIsPrinted = false;
				vector<VectorIndex> path;
				if (rg.vertices[vertexIndex].InDegree() == 0 and rg.vertices[vertexIndex].OutDegree() == 1) {
					path.push_back(vertexIndex);
				}
				VectorIndex pathVertex = rg.edges[rg.vertices[vertexIndex].outEdges[outEdgeIndex]].dest;				
				while(rg.vertices[pathVertex].inEdges.size() == 1 and
							rg.vertices[pathVertex].outEdges.size() == 1) {
					if (rg.vertices[pathVertex].traversed == true) {
						pathIsPrinted = true;
						break;
					}
					path.push_back(pathVertex);
					// Mark the forward and reverse complement as traversed.
					pathVertex = rg.edges[rg.vertices[pathVertex].outEdges[0]].dest;
					//
				}
				//
				// Look to see if this is the end of a simple path, if so, add it to the scaffold.
				//
				pathVertex = rg.edges[rg.vertices[vertexIndex].outEdges[outEdgeIndex]].dest;				
				if (rg.vertices[pathVertex].OutDegree() == 0 and rg.vertices[pathVertex].InDegree() == 1) {
					path.push_back(pathVertex);
				}
				//
				// Determine the sequences in the scaffold and the total scaffold length.
				//
				if (pathIsPrinted == false) {
					VectorIndex p;
					DNALength scaffoldLength = 0;
					for (p = 0; p < path.size(); p++ ){
						scaffoldLength += vertexSequences[path[p]/2].length;
						rg.vertices[path[p]].traversed = true;
						//					rg.vertices[2*(path[p]/2)+ !(path[p]%2)].traversed = true;
						++numPrintedVertices;
					}
					cout << "path is of size " << path.size() << " length " << scaffoldLength << endl;
					if (!printSeparate) {
						scaffoldOut << ">" << scaffoldIndex << " " << path.size() << " " << scaffoldLength << endl;
					}
					for (p = 0; p < path.size(); p++) {
						if (printSeparate) {
							scaffoldOut << ">" << p << " " << path[p]/2 << " " << vertexSequences[path[p]/2].length << endl;
						}
						if (path[p]%2 == 0) {
							((DNASequence)vertexSequences[path[p]/2]).PrintSeq(scaffoldOut);
						}
						else {
							((DNASequence)vertexRCSequences[path[p]/2]).PrintSeq(scaffoldOut);
						}
						rg.vertices[path[p]].traversed = true;
						rg.vertices[2*(path[p]/2) + !(path[p]%2)].traversed = true;
					}
					if (printSeparate) {
						scaffoldOut.close();
						scaffoldOut.clear();
					}
				}
			}
		}
	}

	ofstream* outPtr;
	ofstream repeatOut;
	if (printRepeatsSeparately) {
		CrucialOpen(repeatFileName, repeatOut, std::ios::out);
		outPtr = &repeatOut;
	}
	else {
		outPtr = &scaffoldOut;
	}

	for (vertexIndex = 0; vertexIndex < rg.vertices.size(); vertexIndex++ ){
		if (rg.vertices[vertexIndex].traversed == false) {
			//
			// Print this vertex sequence only.  It is repetitive, or isolated.
			//
			*outPtr << ">" << scaffoldIndex << endl;
			++scaffoldIndex;
			if (vertexIndex%2 == 0) {
				((DNASequence)vertexSequences[vertexIndex/2]).PrintSeq(*outPtr);
			}
			else {
				((DNASequence)vertexRCSequences[vertexIndex/2]).PrintSeq(*outPtr);
			}
			rg.vertices[vertexIndex].traversed = true;
			rg.vertices[2*(vertexIndex/2)+ !(vertexIndex%2)].traversed = true;
		}
	}	

	cout << "printed: " << numPrintedVertices << " of " << rg.vertices.size() << endl;
}
