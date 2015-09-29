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
	cout << "usage: dotplot query target min_k " << endl;
}

void QuickMatch(FASTASequence &target, DNASuffixArray &sa, FASTASequence &query, 
								AnchorParameters &anchorParameters, 
								vector<ChainedMatchPos> &matchPosList) {
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
			for (mp = matchLow[i]; mp < matchHigh[i]; mp++ ) {
				matchPosList.push_back(ChainedMatchPos(sa.index[mp], i, matchLength[i]));
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
	anchorParameters.minMatchLength = minK;
	anchorParameters.stopMappingOnceUnique = false;
	
	QuickMatch(target, sarray, query, anchorParameters,matchPosList);
	QuickMatch(target, sarray, queryRC, anchorParameters, rcMatchPosList);
	
	int i;
	for (i = 0; i < matchPosList.size(); i++ ){
		cout << matchPosList[i].q << "\t" << matchPosList[i].t << "\t" << matchPosList[i].l << "\t0\t0" << endl;
	}

	for (i = 0; i < rcMatchPosList.size(); i++) {
		cout << query.length - (rcMatchPosList[i].q + rcMatchPosList[i].l) << "\t" << rcMatchPosList[i].t << "\t" << rcMatchPosList[i].l << "\t0\t1" << endl;
	}
}
	
