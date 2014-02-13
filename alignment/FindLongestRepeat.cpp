#include "datastructures/suffixarray/SuffixArrayTypes.h"
#include "FASTAReader.h"

#include <iostream>


using namespace std;

int main(int argc, char* argv[]) {

	DNASuffixArray sarray;
	string genomeName = argv[1];
	string sarrayName = genomeName  + ".sa";
	if (sarray.Read(sarrayName)) {
		
	}
	else {
		cout << "Invalid suffix array." << endl;
	}
	
	FASTAReader reader;
	FASTASequence genome;
	reader.Initialize(genomeName);
	reader.GetNext(genome);
	int i;
	int maxRepeat = 0;
	int maxRepeatI = 0;
	genome.ToUpper();
	int maxRepeatJ = 0;
	for (i = 0; i < genome.length-1; i++) {
		int l;
		if (genome.seq[sarray.index[i]] == 'N'  or
				genome.seq[sarray.index[i+1]] == 'N') {
			continue;
		}
		for (l = 0; l < genome.length - max(sarray.index[i], sarray.index[i+1]) - 1 ; l++) {
			if (genome.seq[sarray.index[i]+l] != genome.seq[sarray.index[i+1]+l] or
					genome.seq[sarray.index[i]+l] == 'N' or
					genome.seq[sarray.index[i+1]+l] == 'N') {
				if (l > maxRepeat) {
					maxRepeat = l;
					maxRepeatI = sarray.index[i];
					maxRepeatJ = sarray.index[i+1];
				}
				break;
			}
		}
	}
	cout << "Max repeat: " << maxRepeat << " " << maxRepeatI << endl;
	FASTASequence repeat;
	stringstream titleStrm;
	titleStrm << genome.title << ":" << maxRepeatI << "-" << maxRepeatI + maxRepeat << ";" << genome.title << ":" << maxRepeatJ << "-" << maxRepeatJ + maxRepeat;
	repeat.CopyTitle(titleStrm.str());
	repeat.ReferenceSubstring(genome, maxRepeatI, maxRepeat);
	repeat.PrintSeq(cout);
	return 0;

}
