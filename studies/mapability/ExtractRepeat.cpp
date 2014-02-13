#include "FASTASequence.h"
#include "FASTAReader.h"
#include "algorithms/alignment/SWAlign.h"
#include "algorithms/alignment/KBandAlign.h"
#include "algorithms/alignment/ScoreMatrices.h"
#include <string>
#include <sstream>
using namespace std;
int main(int argc, char* argv[]) {
	string genomeFileName, subseqFileName;
	if (argc != 3) {
		cout << "usage: extractRepeats genome repeat" << endl;
		exit(1);
	}

	genomeFileName = argv[1];
	subseqFileName = argv[2];

	FASTASequence genome, sub;
	FASTAReader reader;
	reader.Init(genomeFileName);
	reader.GetNext(genome);
	reader.Init(subseqFileName);
	reader.GetNext(sub);

	genome.ToUpper();
	sub.ToUpper();	
	DNALength genomePos;
	FASTASequence genomeSub;
	int kband = (int) (0.15) * sub.length;
	vector<int> scoreMat;
	vector<Arrow> pathMat;
	int readIndex = 0;
	cout << "starting extraction" << endl;
	for (genomePos = 0; genomePos < genome.length - sub.length + 1; genomePos++) {
		genomeSub.seq = &genome.seq[genomePos];
		genomeSub.length = sub.length;
		int alignScore;
		Alignment alignment;
		alignScore = SWAlign(genomeSub, sub,
												 EditDistanceMatrix, 1, //1,kband,
												 scoreMat, pathMat,
												 alignment, QueryFit);
												 
		if (alignScore < 0.25 * sub.length) {
			stringstream titlestrm;
			titlestrm << readIndex << "|" 
								 << genomePos << "|"
								<< genomePos + sub.length << " " << alignScore/ (1.0*sub.length);
			FASTASequence subcopy;
			subcopy.CopyTitle(titlestrm.str());
			subcopy.seq = &genome.seq[genomePos];
			subcopy.length = sub.length;
			subcopy.PrintSeq(std::cout);
			genomePos += sub.length;
		}
	}
}


	
	
