#include "../common/FASTASequence.h"
#include "../common/FASTAReader.h"
#include "../common/utils.h"
#include "../common/statistics/statutils.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

int main(int argc, char *argv[]) {
	string sequencesInName, sequencesOutName;
	if (argc <3){ 
		cout << "usage: scramble in out" << endl;
		exit(1);
	}
	sequencesInName = argv[1];
	sequencesOutName= argv[2];
	vector<FASTASequence*> sequences;
	vector<int> sequenceIndices;

	FASTAReader reader;
	reader.Init(sequencesInName);
	ofstream out;
	CrucialOpen(sequencesOutName, out, std::ios::out);
	

	FASTASequence read;
	FASTASequence*readPtr;
	while(reader.GetNext(read)) {
		readPtr = new FASTASequence;
		*readPtr = read;
		sequences.push_back(readPtr);
	}

	int i;
	for (i = 0; i < sequences.size(); i++) {
		sequenceIndices.push_back(i);
	}

	for (i = 0; i < 10*sequences.size(); i++ ){
		//
		// shuffle indices.
		//
		int idx1;
		int idx2;
		idx1 = RandomInt(sequences.size());
		idx2 = RandomInt(sequences.size());
		int tmp;
		tmp  = sequenceIndices[idx1];
		sequenceIndices[idx1] = sequenceIndices[idx2];
		sequenceIndices[idx2] = tmp;
	}
	
	for (i = 0; i < sequenceIndices.size(); i++ ){
		sequences[sequenceIndices[i]]->PrintSeq(out);
	}
	return 0;
}
		
	
		
	
