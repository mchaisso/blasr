#include <assert.h>
#include <string>
#include <iostream>
#include <vector>

#include "../common/FASTAReader.h"
#include "../common/FASTASequence.h"
#include "../common/algorithms/alignment.h"
#include "../common/algorithms/alignment/AffineKBandAlign.h"
#include "../common/algorithms/alignment/DistanceMatrixScoreFunction.h"
#include "../common/datastructures/alignment/AlignmentCandidate.h"
using namespace std;

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: sdpMatcher query target [-k k] [-ins i] [-del d]" << endl;
		exit(1);
	}

	string queryName, targetName;
	queryName = argv[1];
	targetName = argv[2];
	int argi = 3;
	int k;
	int ins = 3;
	int del = 3;
  bool targetfixed = false;
	while (argi < argc ) {
		if (strcmp(argv[argi], "-k") == 0) {
			k = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-del") == 0) {
			del = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-ins") == 0) {
			ins = atoi(argv[++argi]);
		}
    else if (strcmp(argv[argi], "-targetfixed") == 0) {
      targetfixed = true;
    }
    else {
      cout << "ERROR: bad option: " << argv[argi] << endl;
      exit(1);
    }
		++argi;
	}


	FASTASequence query, target;
	FASTAReader queryReader, targetReader;
	queryReader.Init(queryName);
	
	targetReader.Init(targetName);
	//
	// Prepare the target database;
	//

	//
	// Prepare the query match set.
	//

	int seqIndex = 0;
	int numNoFragments = 0;

	vector<int> scoreMat;
	vector<Arrow> pathMat;
	int alignScore;

	int computedAlignScore;
	DistanceMatrixScoreFunction<FASTASequence, FASTASequence> scoreFn;
	scoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);
  if (targetfixed) {
    targetReader.GetNext(target);
  }

	while (queryReader.GetNext(query) and (targetfixed or targetReader.GetNext(target))) {
    cout << "got " << query.title << endl;
		AlignmentCandidate<FASTASequence, FASTASequence> alignment, aalignment;
		alignment.blocks.clear();
		alignment.qPos = 0;
		alignment.tPos = 0;
		if (query.length == 0 or target.length == 0)
			continue;

		vector<int> kscoremat, ascoremat;
		vector<Arrow> kpathmat, apathmat;
		vector<int> hpInsScoreMat, insScoreMat;
		vector<Arrow> hpInsPathMat, insPathMat;
		
		alignScore = AffineKBandAlign(query, target, SMRTDistanceMatrix,
																	ins, ins - 2, // homopolymer affine penalty
																	ins, ins - 1, // regular affine
																	del, k, 
																	ascoremat, apathmat,
																	hpInsScoreMat, hpInsPathMat,
																	insScoreMat, insPathMat,
																	aalignment, Local);

		aalignment.qPos = 0;
		aalignment.tPos = 0;

		PrintCompareSequencesAlignment(aalignment, query.seq, target.seq, cout, false);
		PrintAlignmentStats(aalignment, cout);			
		++seqIndex;
	}

	return 0;
}


 
