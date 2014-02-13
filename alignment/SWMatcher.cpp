#include <assert.h>
#include <string>
#include <iostream>
#include <vector>

#include "FASTAReader.h"
#include "FASTASequence.h"
#include "algorithms/alignment.h"
#include "algorithms/alignment/DistanceMatrixScoreFunction.h"
#include "algorithms/alignment/IDSScoreFunction.h"

using namespace std;

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: swMatcher query target [-indel i] [-local] [-showalign] " << endl
				 << "       [-type queryfit|overlap|global] [-match m ] [-mismatch m]" << endl
				 << "    or [-local] [-queryfit] [-overlap] [-fixedtarget] [-fixedquery]" << endl
				 << "       [-printmatrix]"<< endl;
    cout << "   Unless -showalign is specified, output is tabular and in the formt:"<<endl
         << "   query_length target_length align_score query_start query_end target_start target_end"<<endl;
		exit(1);
	}

	string queryName, targetName;
	queryName = argv[1];
	targetName = argv[2];
	int argi = 3;
	int indelCost = 3;
	int showAlign = 0;
	AlignmentType alignType = Global;
	int match = 0;
	int mismatch = 0;
	int fixedTarget = 0;
	int fixedQuery  = 0;
	bool printMatrix = false;
  int insertion = 4;
  int deletion  = 5;
	while (argi < argc) {
		if (strcmp(argv[argi], "-insertion") == 0) {
			insertion = atoi(argv[++argi]);
		}
    else if (strcmp(argv[argi], "-deletion") == 0) {
      deletion = atoi(argv[++argi]);
    }
		else if (strcmp(argv[argi], "-local") == 0) {
			alignType = Local;
		}
		else if (strcmp(argv[argi], "-showalign") == 0) {
			showAlign = 1;
		}
		else if (strcmp(argv[argi], "-fixedtarget") == 0) {
			fixedTarget = 1;
		}
		else if (strcmp(argv[argi], "-fixedquery") == 0) {
			fixedQuery = 1;
		}
		else if (strcmp(argv[argi], "-type") == 0) {
			++argi;
			if (strcmp(argv[argi], "queryfit") == 0) {
				alignType = QueryFit;
			}
			else if (strcmp(argv[argi], "targetfit") == 0) {
				alignType = TargetFit;
			}
			else if (strcmp(argv[argi], "overlap") == 0) {
				alignType = Overlap;
			}
			else if (strcmp(argv[argi], "global") == 0) {
				alignType = Global;
			}
			else if (strcmp(argv[argi], "tpqs") == 0) {
				alignType = TPrefixQSuffix;
			}
			else if (strcmp(argv[argi], "tsqp") == 0 ){ 
				alignType = TSuffixQPrefix;
			}
			else {
				cout <<" ERROR, aligntype must be one of queryfit, overlap, or global" << endl;
				exit(1);
			}
		}
		else if(strcmp(argv[argi], "-printmatrix") == 0) {
			printMatrix = true;
		}
		else if (strcmp(argv[argi], "-local") == 0) {
			alignType = Local;
		}
		else if (strcmp(argv[argi], "-queryfit") == 0) {
			alignType = QueryFit;
		}
		else if (strcmp(argv[argi], "-overlap") == 0) {
			alignType = Overlap;
		}
		else if (strcmp(argv[argi], "-match") == 0) {
			match = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-mismatch") == 0) {
			mismatch = atoi(argv[++argi]);
		}
		++argi;
	}
	DistanceMatrixScoreFunction<FASTASequence, FASTASequence> scoreFn;
	scoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);
	scoreFn.ins = insertion;
	scoreFn.del = deletion;

	FASTASequence query, target;
	FASTAReader queryReader, targetReader;
	queryReader.Init(queryName);
	
	targetReader.Init(targetName);

	if (fixedTarget) {
		targetReader.GetNext(target);
	}
	if (fixedQuery) {
		queryReader.GetNext(query);
	}
	//
	// Prepare the target database;
	//

	//
	// Prepare the query match set.
	//

	int seqIndex = 0;

	vector<int> scoreMat;
	vector<Arrow> pathMat;
	int alignScore;
	MatchedAlignment alignment;

	if (match != 0) {
		int i;
		for (i = 0; i < 4; i++ ) {
			LocalAlignLowMutationMatrix[i][i] = match;
		}
	}

	int i,j;
	for (i = 0; i < 5; i++) {
		for (j = 0; j < 5 ; j++) {
			if (i == j) continue;
			SMRTDistanceMatrix[i][j] += 3;
		}
	}

	cout << "qlen tlen score" << endl;
	while ((fixedQuery or queryReader.GetNext(query)) and 
				 (fixedTarget or targetReader.GetNext(target))) {
		alignment.qName.assign(query.title, query.titleLength);
		alignment.tName.assign(target.title, target.titleLength);
		alignment.blocks.clear();
		alignment.qPos = 0;
		alignment.tPos = 0;
		alignment.qStart = 0;
		alignment.tStart = 0;
		if (query.length == 0 or target.length == 0)
			continue;

		alignScore = SWAlign(query, target, scoreMat, pathMat, 
												 alignment, scoreFn, alignType, false, printMatrix);

		cout << query.length << " " << target.length << " " << alignScore << endl;
		cout << alignment.qPos << " " << alignment.QEnd() 
         << " " << alignment.tPos << " " << alignment.TEnd() << endl;

		if (showAlign) {
			ComputeAlignmentStats(alignment, query.seq, target.seq, scoreFn);
			PrintAlignmentStats(alignment, cout);			
			StickPrintAlignment(alignment, query, target, cout);
		}
		++seqIndex;
	}

	return 0;
}


