#include <string>
#include <iostream>
#include "../common/utils.h"
#include "../common/FASTASequence.h"
#include "../common/FASTAReader.h"

int main(int argc, char* argv[]) {
	
	if (argc < 3 ){ 
		cout << "usage: catreads readsFile catFile " << endl;
		cout << "       Joins all reads in reads file into one fasta sequence, with an 'N'" <<endl
				 << "       separating each read." <<endl;
		exit(1);
	}
	string readsFileName = argv[1];
	string catFileName   = argv[2];

	FASTAReader reader;
	reader.Init(readsFileName);
	FASTASequence catSeq;
	ofstream catFile;
	CrucialOpen(catFileName, catFile);
	reader.ReadAllSequencesIntoOne(catSeq);
	catSeq.PrintSeq(catFile);
	return 0;
}
	

	
