#include "../common/FASTAReader.h"
#include "../common/FASTASequence.h"
#include "../common/utils.h"
#include <string>
#include <sstream>
using namespace std;

int main(int argc, char* argv[]) {
	string seqInName, seqOutName, dotOutName;
	if (argc < 4) {
		cout << "usage: exciseRepeats inName repMaskOutFile outName" << endl;
		exit(1);
	}

	seqInName = argv[1];
	dotOutName = argv[2];
	seqOutName = argv[3];
	FASTAReader reader;
	reader.Initialize(seqInName);
	FASTASequence origSeq;
	reader.GetNext(origSeq);
	
	ifstream dotOutFile;
	CrucialOpen(dotOutName, dotOutFile);
	ofstream seqOutFile;
	ofstream seqOut;
	CrucialOpen(seqOutName, seqOut, std::ios::out);
	string dotOutLine;
	getline(dotOutFile, dotOutLine);
	getline(dotOutFile, dotOutLine);
	getline(dotOutFile, dotOutLine);
	while(getline(dotOutFile, dotOutLine)) {
		stringstream lineStrm(dotOutLine);
		int swScore;
		float pctDiv, pctDel, pctIns;
		string query;
		int qPosBegin, qPosEnd;
		string left;
		char strand;
		string matchingRepeat;
		string repClass;
		string repPos, repEnd, repLeft;
		int id;
		lineStrm >> swScore >> pctDiv >> pctDel >> pctIns >> query >> qPosBegin >> qPosEnd >> left >> strand >> matchingRepeat >> repClass >> repPos >> repEnd >> repLeft >> id;
		DNALength seqPos;
		for (seqPos = qPosBegin; seqPos < qPosEnd; seqPos++) {
			origSeq.seq[seqPos] = 'X';
		}
	}

	DNALength seqPos, unexPos;
	unexPos = 0;
	for (seqPos = 0; seqPos < origSeq.length; seqPos++) {
		if (origSeq.seq[seqPos] != 'X') {
			origSeq.seq[unexPos] = origSeq.seq[seqPos];
			unexPos++;
		}
	}
	origSeq.length = unexPos;

	origSeq.PrintSeq(seqOut);
	return 0;
}

		

	
					 
					 
