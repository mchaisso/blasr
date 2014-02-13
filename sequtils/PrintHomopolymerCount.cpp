#include "../common/FASTAReader.h"
#include "../common/FASTASequence.h"


#include <map>
#include <string>
#include <iostream>

int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << "usage: printHPC in.fa" << endl;
		exit(1);
	}
	string seqFileName = argv[1];

	FASTAReader reader;
	reader.Init(seqFileName);
	FASTASequence seq;
	reader.GetNext(seq);
	int p, pn;
	pn = 0;
	map<int,int> hpCount;
	for (p = 0; pn < seq.length; p = pn) {
		pn = p + 1;
		while(pn < seq.length and seq.seq[p] == seq.seq[pn]) {
			++pn;
		}
		hpCount[pn-p]++;
	}
	map<int,int>::iterator mapIt;
	int first = 2;
	int sum = 0;
	for( mapIt = hpCount.begin(); mapIt != hpCount.end(); ++mapIt) {
		if (first == 0) {
			sum += mapIt->second;
		}
		else {
			--first;
		}
		cout << mapIt->first << " " << mapIt->second << " " << sum << endl;
	}
	return 0;
}
		
		
	
