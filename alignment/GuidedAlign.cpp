#include "../common/algorithms/alignment/GuidedAlign.h"
#include "../common/algorithms/alignment/SDPAlign.h"
#include "../common/algorithms/alignment/IDSScoreFunction.h"
#include "../common/algorithms/alignment/DistanceMatrixScoreFunction.h"
#include "../common/algorithms/alignment/AlignmentPrinter.h"
#include "../common/datastructures/alignment/AlignmentCandidate.h"
#include "../common/SMRTSequence.h"
#include "../common/FASTQSequence.h"
#include "../common/defs.h"
#include "../common/files/ReaderAgglomerate.h"
#include "../common/algorithms/alignment/ScoreMatrices.h"

using namespace std;
int main(int argc, char* argv[]) {
	
	string queryFileName, targetFileName;
	if (argc < 3) {
		cout << "Usage: guidedalign query target [sdptuple]" << endl;
		exit(1);
	}
	queryFileName = argv[1];
	targetFileName = argv[2];
	int sdpTupleSize = 4;
	if (argc > 3) {
		sdpTupleSize = atoi(argv[3]);
	}
	
	ReaderAgglomerate reader;
	FASTQSequence query, target;

	reader.Initialize(queryFileName);
	reader.GetNext(query);
	reader.Close();
	reader.Initialize(targetFileName);
	reader.GetNext(target);
	reader.Close();
	
	int alignScore;
	/*
	Alignment sdpAlignment;
	int nSDPHits = 0;
	alignScore = SDPAlign(query, target,
												SMRTDistanceMatrix, 
												4, 4, sdpTupleSize, 4, 0.90,
												sdpAlignment, nSDPHits, Local, false, false);
	int b;
	for (b = 0; b < sdpAlignment.blocks.size(); b++) {
		sdpAlignment.blocks[b].qPos += sdpAlignment.qPos;
		sdpAlignment.blocks[b].tPos += sdpAlignment.tPos;
		}
	Guide guide;
	int bandSize = 16;
	AlignmentToGuide(sdpAlignment, guide, bandSize);
	StoreMatrixOffsets(guide);
	int guideSize = ComputeMatrixNElem(guide);
	int i;
	*/

	vector<int> scoreMat;
	vector<Arrow> pathMat;
	vector<double> probMat, optPathProbMat;
  vector<float> lnSubVect, lnInsVect, lnDelVect, lnMatchVect;
  //	AlignmentCandidate<FASTASequence, FASTASequence> alignment;
  Alignment alignment;
	DistanceMatrixScoreFunction<DNASequence, DNASequence> distScoreFn;
	distScoreFn.del = 3;
	distScoreFn.ins = 3;
	distScoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);

	alignScore = GuidedAlign(query, target, distScoreFn, 10,
                           // in order after edit distance:
                           // pairwise-ins, pairwise-del, k, sdp-ins, sdp-del, sdp-insrate
                           //                           distScoreFn, 

                           5,5,.15,
                           alignment, Local, false, 8);
  //	StickPrintAlignment(alignment, query, target, cout);
}

