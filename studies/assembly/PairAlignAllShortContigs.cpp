#include <string>
#include <vector>
#include <iostream>
#include "FASTAReader.h"
#include "FASTASequence.h"
#include "algorithms/alignment/SWAlign.h"
#include "algorithms/alignment/ScoreMatrices.h"
#include "algorithms/alignment/printers/StickAlignmentPrinter.h"
#include "datastructures/matrix/Matrix.h"
#include "utils.h"

using namespace std;


int main(int argc, char* argv[]) {


	FASTAReader reader;
	FASTASequence read;
	int maxLength = 100;
	if (argc < 3) {
		cout << "usage: pairAlignAllContigs inFile maxLength equivalencies [-minIdent i]" << endl; 
		exit(1);
	}
	string readsFileName, equivalenciesFileName;
	readsFileName = argv[1];
	maxLength = atoi(argv[2]);
	equivalenciesFileName = argv[3];
	int argi = 4;
  float minIdentity = 80;
	while (argi < argc) {
		if (strcmp(argv[argi], "-minIdent") == 0) {
			minIdentity = atoi(argv[++argi]);
		}
		++argi;
	}
	vector<FASTASequence> reads, readsRC;;
	reader.Init(readsFileName); 
  reader.ReadAllSequences(reads);
	readsRC.resize(reads.size());
	int r;
	for (r =0; r < reads.size();r++) {
	  reads[r].MakeRC(readsRC[r]);
  }
	ofstream equivOut;
	CrucialOpen(equivalenciesFileName, equivOut);

	Matrix<int> alignScores;
	Matrix<float> alignIdentities;
	alignScores.Resize(reads.size(), reads.size());
	alignIdentities.Resize(reads.size(), reads.size());
	vector<int> scoreMat;
	vector<Arrow> pathMat;
	int i, j;
	int alignScore;
	FASTASequence readi, readj;
	FASTASequence rcReadi, rcReadj;
	
	for (i = 0; i < reads.size(); i++) {
		float maxFrontIdent, maxEndIdent;
		int   maxFrontIdentIndex, maxEndIdentIndex;
		maxFrontIdent = 0; maxEndIdent = 0;
		maxFrontIdentIndex = 0;
		maxEndIdentIndex   = 0;
		int maxFrontIdentLength = 0;
		int maxEndIdentLength  = 0;
		int maxFrontLength     = 0;
		int	maxEndLength       = 0;
		int nmaxFrontLengthIndex = 0;
		int maxEndLengthIndex  = 0;
		float maxFrontLengthIdent = 0;
		float maxEndLengthIdent = 0;
		int maxFrontLengthIndex = 0;
		equivOut << reads[i].GetName();
		for (j = 0; j < reads.size(); j++ ){
			// 
			// Store the two ends of the alignment.
			//
			alignScore = 0;
			int rcAlignScore;
			Alignment alignment;
			Alignment rcAlignment;
			Alignment *optAlignment;
			if (i != j) {
  	  	if (reads[i].length < maxLength and reads[j].length < maxLength) {
  				alignScore = SWAlign(reads[i], reads[j], SMRTDistanceMatrix, 3, scoreMat, pathMat, alignment, Global);
  			}
				if (reads[i].length < maxLength and reads[j].length < maxLength) {
          rcAlignScore = SWAlign(reads[i], readsRC[j], SMRTDistanceMatrix, 3, scoreMat, pathMat, rcAlignment, Global);
        }	
  			ComputeAlignmentStats(alignment, reads[i].seq, reads[j].seq, SMRTDistanceMatrix, 3,3 );
        ComputeAlignmentStats(rcAlignment, reads[i].seq, readsRC[j].seq, SMRTDistanceMatrix, 3,3 );

  			if (alignment.pctSimilarity > minIdentity or rcAlignment.pctSimilarity > minIdentity) {
  				equivOut << " " << reads[j].GetName();
  			}	
      }
	}
		equivOut << endl;	
	}

	return 0;
}
