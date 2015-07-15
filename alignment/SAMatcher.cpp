#include "FASTASequence.h"
#include "FASTQSequence.h"
#include "FASTAReader.h"
#include "datastructures/suffixarray/SuffixArray.h"
#include "datastructures/suffixarray/SuffixArrayTypes.h"
#include "algorithms/alignment/AlignmentPrinter.h"
#include "algorithms/alignment/sdp/VariableLengthSDPFragment.h"
#include "algorithms/alignment/sdp/NonoverlappingSparseDynamicProgramming.h"
#include "algorithms/anchoring/MapBySuffixArray.h"
#include "algorithms/anchoring/GlobalChain.h"
#include "algorithms/anchoring/BasicEndpoint.h"
#include "datastructures/anchoring/MatchPos.h"
#include "datastructures/anchoring/AnchorParameters.h"

// Declare some static variables.


int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: samatcher queryfile targetfile" <<endl;
		exit(1);
	}

	string queryFileName, targetFileName;
	int minMatchLength = 5;
	int maxExpand      = 0;
	queryFileName  = argv[1];
	targetFileName = argv[2];
	int argi = 3;
	AnchorParameters anchorParams;
	
	while (argi < argc) {
		if (strcmp(argv[argi], "-minmatch") == 0) {
			anchorParams.minMatchLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-maxexpand") == 0) {
			anchorParams.expand = atoi(argv[++argi]);
		}
		else {
			cout << "ERROR! Invalid argument: " << argv[argi]<< endl;
			exit(1);
		}
		++argi;
	}

	FASTQSequence query, target;
	FASTAReader queryReader, targetReader;
	queryReader.Init(queryFileName);
	targetReader.Init(targetFileName);


	while(1) {
		if (!queryReader.GetNext(query)) break;
		if (!targetReader.GetNext(target)) break;

		query.ToUpper();
		target.ToUpper();

		//
		// Build the suffix array on the target.
		//
		DNASuffixArray sarray;
		target.ToThreeBit();		
		vector<int> alphabet;
		sarray.InitThreeBitDNAAlphabet(alphabet);
		sarray.LarssonBuildSuffixArray(target.seq, target.length, alphabet);
		cout <<"done building suffix array." << endl;
		target.ToAscii();

		//
		// Find the list of anchors.
		//
	
		query.PrintSeq(cout);
		cout << "target: " << endl;
		target.PrintSeq(cout);
		MatchPosList matchPosList;
		int numKeysMatched;
		anchorParams.useLookupTable = false;
		numKeysMatched   = 
			MapReadToGenome(target, sarray, query, sarray.lookupPrefixLength,
											matchPosList, anchorParams);
	
		//
		// Now, convert the matchPosList to a set of fragments
		// that can be used in the sdp.
		//
		SortMatchPosList(matchPosList);
		vector<ChainedFragment> fragments;
		fragments.resize(matchPosList.size());
		VectorIndex i;
		for (i = 0; i < matchPosList.size(); i++) {
			fragments[i].x = matchPosList[i].t;
			fragments[i].y = matchPosList[i].q;
			fragments[i].length = fragments[i].weight = matchPosList[i].w;
			//		cout << fragments[i].x << " " << fragments[i].y << " " << fragments[i].weight << endl;
		}
		cout << "stored a total of: " << fragments.size() << " fragments." << endl;

		int maxFragmentChainLength;
		vector<DNALength> maxFragmentChain;
	
		maxFragmentChainLength = GlobalChain<ChainedFragment, BasicEndpoint<ChainedFragment> >(fragments, maxFragmentChain);
		MatchedAlignment alignment;		
		std::reverse(maxFragmentChain.begin(), maxFragmentChain.end());
		alignment.AllocateBlocks(maxFragmentChain.size());
		for (i = 0; i < maxFragmentChain.size(); i++) {
			alignment.blocks[i].qPos = fragments[maxFragmentChain[i]].y;
			alignment.blocks[i].tPos = fragments[maxFragmentChain[i]].x;
			alignment.blocks[i].length = fragments[maxFragmentChain[i]].length;
			cout << "( " << fragments[maxFragmentChain[i]].x << " "
					 << fragments[maxFragmentChain[i]].y << " "
					 << fragments[maxFragmentChain[i]].length << ") ";
		}
		cout << endl;
		alignment.tStart = alignment.qStart = 0;
		alignment.tPos = alignment.qPos = 0;
		StickPrintAlignment(alignment, query, target, cout);
	}
	return 0;
}
