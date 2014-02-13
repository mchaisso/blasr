#include "FASTAReader.h"
#include <iostream>
#include <sstream>

#include "utils.h"

int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << "usage: extractGFFSequences gffFile [seqDir]" << endl; 
		cout << " seqDir should point to the directory relative to the current " << endl
				 << "    directory where the reference sequences are" << endl;
		exit(1);
	}
	string gffFileName = argv[1];
  string dirName;
	dirName = ".";
	if (argc > 2) {
		dirName = argv[2];
	}
	
	string curSeqName, curFileName;
	curSeqName = "";
	
	ifstream gffIn;
	CrucialOpen(gffFileName, gffIn, std::ios::in);
	FASTAReader reader;
	FASTASequence seq;
	while(gffIn) {
		string line;
		getline(gffIn, line);
		stringstream linestrm(line);
		string seqName, dupType, sim;
		int startPos, endPos;
		char strand;

		linestrm >> seqName >> dupType >> sim >> startPos >> endPos;
		if (seqName != curSeqName) {
			if (curSeqName != "" ){ 
				reader.Close();
				seq.Free();
			}
			curFileName = dirName;
			curFileName.append("/");
			curFileName.append(seqName);
			curFileName.append(".fa");
			curSeqName = seqName;
			reader.Init(curFileName);
			reader.GetNext(seq);
		}

		FASTASequence dupSeq;
		stringstream dupSeqTitleStrm;
		dupSeqTitleStrm << curSeqName << "_" << startPos << "_" << endPos;
		dupSeq.CopyTitle(dupSeqTitleStrm.str());
		dupSeq.seq = &seq.seq[startPos];
		dupSeq.length = endPos - startPos;
		dupSeq.PrintSeq(cout);
	}
}
