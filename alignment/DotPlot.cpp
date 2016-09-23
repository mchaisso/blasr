#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <sstream>
#include <pthread.h>
#include <stdlib.h>
#include <time.h>
#include <signal.h>
#include <execinfo.h>
#include "FASTASequence.h"
#include "SMRTSequence.h"
#include "FASTAReader.h"
#include "SeqUtils.h"
#include "algorithms/anchoring/MapBySuffixArray.h"
#include "datastructures/suffixarray/SharedSuffixArray.h"
#include "datastructures/suffixarray/SuffixArrayTypes.h"
#include "datastructures/anchoring/AnchorParameters.h"	
#include "tuples/DNATuple.h"

typedef DNASuffixArray T_SuffixArray;

void PrintUsage() {
	cout << "usage: dotplot query target min_k [--maxCount] [--useLcp]" << endl;
	cout << "  --maxCount Limits the number of matches per position." << endl
			 << "  --useLcp   Use the longest match at every position (will likely not show STR)." << endl;
}

void QuickMatch(FASTASequence &target, DNASuffixArray &sa, FASTASequence &query, 
								AnchorParameters &anchorParameters, 
								vector<ChainedMatchPos> &matchPosList, int maxMatchesPerPosition=0) {
	SMRTSequence smrtQuery;
	smrtQuery.seq = query.seq;
	smrtQuery.length = query.length;
	smrtQuery.subreadStart = 0;
	smrtQuery.subreadEnd = query.length;
	vector<DNALength> matchLow, matchHigh, matchLength;
	LocateAnchorBoundsInSuffixArray(target, sa, smrtQuery, 0, matchLow, matchHigh, matchLength, anchorParameters);	

	int i;
	for (i = 0; i < matchLow.size(); i++) {
		DNALength mp;
		if (matchLength[i] >= anchorParameters.minMatchLength) {
			if (maxMatchesPerPosition == 0 or matchHigh[i] - matchLow[i] <= maxMatchesPerPosition) {
				for (mp = matchLow[i]; mp < matchHigh[i]; mp++ ) {
					matchPosList.push_back(ChainedMatchPos(sa.index[mp], i, matchLength[i]));
				}
			}
		}
	}
	vector<bool> toRemove(matchPosList.size(), false);
	i = 0;
	while (i < matchPosList.size()) {
		int l=1;
		while (i + l < matchPosList.size() and
					 matchPosList[i].q + l == matchPosList[i+l].q and 
					 matchPosList[i].t + l == matchPosList[i+l].t and
					 matchPosList[i].l == matchPosList[i+l].l + l) {
			toRemove[i+l] = true;
			l++;
		}
		i+=l;
	}
	int j;
	i = 0;
	for (j = 0; j < matchPosList.size(); j++ ){
		if (toRemove[j] == false) {
			matchPosList[i] = matchPosList[j];
			i++;
		}
	}
	matchPosList.resize(i);
	cerr <<"done with " << query.title << " " << matchPosList.size() << endl;
}

int main(int argc, char* argv[]) {
	string queryName, targetName;

	int minK;

	if (argc < 4) {
		PrintUsage();
		exit(1);
	}

	queryName = argv[1];
	targetName = argv[2];
	minK = atoi(argv[3]);
	string outFileName = "/dev/stdout";
	ofstream outFile;
	int argi = 4;
	int maxPerPosition=0;
	bool useLCP = false;
	while (argi < argc){ 
		if (strcmp(argv[argi], "-o") == 0) {
			++argi;
			cerr << "opening " << argv[argi] << endl;
			outFileName = argv[argi];
		}
		if (strcmp(argv[argi], "-m") == 0 or strcmp(argv[argi], "--maxCount") == 0) {
			++argi;
			maxPerPosition = atoi(argv[argi]);
		}
		if (strcmp(argv[argi], "--useLCP") == 0) {
			++argi;
			useLCP = true;
		}
		
		
		++argi;
	}

	ofstream dotOut(outFileName.c_str());
	FASTAReader reader;
	reader.Initialize(queryName);
	
	FASTASequence query, queryRC, target;
	reader.GetNext(query);
	reader.Initialize(targetName);
	reader.GetNext(target);

  vector<ChainedMatchPos> matchPosList;
  vector<ChainedMatchPos> rcMatchPosList;

	query.MakeRC(queryRC);

	//
	// Build the suffix array.
	//
	DNASuffixArray sarray;

	target.ToThreeBit();		
	vector<int> alphabet;
	sarray.InitThreeBitDNAAlphabet(alphabet);
	sarray.LarssonBuildSuffixArray(target.seq, target.length, alphabet);
	target.ConvertThreeBitToAscii();
	cerr << "done building sa" << endl;
	AnchorParameters anchorParameters;
	anchorParameters.minMatchLength        = minK;
	if (useLCP == false) {
		anchorParameters.maxLCPLength          = minK+1;
	}
	else {
		anchorParameters.maxLCPLength          = 0;
	}
	
	anchorParameters.stopMappingOnceUnique = true;
 
	QuickMatch(target, sarray, query, anchorParameters,matchPosList, maxPerPosition);
	QuickMatch(target, sarray, queryRC, anchorParameters, rcMatchPosList, maxPerPosition  );
	
	int i;
	for (i = 0; i < matchPosList.size(); i++ ){
		dotOut << matchPosList[i].q << "\t" << matchPosList[i].t << "\t" << matchPosList[i].l << "\t0\t0" << endl;
	}

	for (i = 0; i < rcMatchPosList.size(); i++) {
		dotOut << query.length - rcMatchPosList[i].q << "\t" << rcMatchPosList[i].t << "\t" << rcMatchPosList[i].l << "\t1\t1" << endl;
	}

	
	dotOut.close();
}
