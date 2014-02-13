#include <iostream>
#include <string>
#include "FASTAReader.h"
#include "FASTASequence.h"
#include "FileUtils.h"
#include "datastructures/suffixarray/SuffixArray.h"
#include "datastructures/suffixarray/ssort.h"
#include "algorithms/sorting/qsufsort.h"

using namespace std;

void PrintUsage() {
	cout << "usage: testSuffixArray in.fa [-slow] [-print] [-query] [-lcp] [-blt] [-ps N] [-mcilroy] [-larsson]"<<endl;
	cout << "  -slow runs the n^2 log n method for comparison purposed only" << endl;
	cout << "  -mcilroy runs the ultra fast McIlroy SA implementation" << endl;
	cout << "  -print prints the index to cout" << endl;
	cout << "  -query queries all entries in the query list against the index." << endl;
	cout << "  -lcp   query for the lcp between the entries in the list and the index." << endl;
	cout << "  -blt   Builds a lookup table to spped queries, mmm." << endl;
	cout << "  -ps    Prints the suffices up to N characters." << endl;
}

int main(int argc, char* argv[]) {

	string inName, outName;
	if (argc < 2) {
		PrintUsage();
		exit(1);
	}
	inName = argv[1];
//	outName = argv[2];
	SAType satype = manmy;

	int argi  = 2;
	int print = 0;
	string queryFileName = "";
	string lcpFileName   = "";
	int doQuery = 0;
	int doBlt   = 0;
	int doPs    = 0;
	int doLCP   = 0;
	int printSuffixLength = 0;
	int ltPrefixLength = 0;
	
	while (argi < argc ) {
		if (strcmp(argv[argi], "-slow") == 0) {
			satype = slow;
		}
		else if (strcmp(argv[argi], "-print") == 0) {
			print = 1;
		}
		else if (strcmp(argv[argi], "-query") == 0) {
			doQuery = 1;
			++argi;
			queryFileName = argv[argi];
		}
		else if (strcmp(argv[argi], "-lcp") == 0) {
			doLCP = 1;
			++argi;
			lcpFileName = argv[argi];
		}
		else if (strcmp(argv[argi], "-blt") == 0) {
			ltPrefixLength = atoi(argv[++argi]);
			doBlt = 1;
		}
		else if (strcmp(argv[argi], "-ps") == 0) {
			printSuffixLength = atoi(argv[++argi]);
			doPs = 1;
		}
		else if (strcmp(argv[argi], "-mcilroy") == 0) {
			satype = mcilroy;
		}
		else if (strcmp(argv[argi], "-larsson") == 0) {
			satype = larsson;
		}
		else {
			PrintUsage();
			cout << "bad option: " << argv[argi] << endl;
			exit(1);
		}
		++argi;
	}
	FASTAReader reader;
	reader.Init(inName);
	/*
	ofstream out;
	CriticalOpenWrite(outName, out);
	*/
	vector<int> dnaAlphabet;
	SuffixArray<Nucleotide, vector<int> > sarray, slowsarray;
	
	FASTASequence seq;
	reader.GetNext(seq);
	seq.ToTwoBit();
	if (satype == manmy) {
		sarray.MMBuildSuffixArray(seq.seq, seq.length, dnaAlphabet);
	}
	else if (satype == slow) {
		sarray.BuildSuffixArray(seq.seq, seq.length, dnaAlphabet);
	}
	else if (satype == mcilroy) {
		sarray.index = new SAIndex[seq.length+1];
		DNALength i;
		for (i = 0; i < seq.length; i++) { sarray.index[i] = seq.seq[i] + 1;}
		sarray.index[seq.length] = 0;
		ssort(sarray.index, NULL);
		for (i = 1; i < seq.length+1; i++ ){ sarray.index[i-1] = sarray.index[i];};
		sarray.length = seq.length;
	}
	else if (satype == larsson) {
		sarray.index = new SAIndex[seq.length+1];
		SAIndex *p = new SAIndex[seq.length+1];
		DNALength i;
		for (i = 0; i < seq.length; i++) { sarray.index[i] = seq.seq[i] + 1;}
		sarray.index[seq.length] = 0;
		vector<int> alphabet;
	
		//	sa.InitTwoBitDNAAlphabet(alphabet);
		//	sa.InitAsciiCharDNAAlphabet(alphabet);
		sarray.InitThreeBitDNAAlphabet(alphabet);

		sarray.LarssonBuildSuffixArray(seq.seq, seq.length, alphabet);

		for (i = 1; i < seq.length+1; i++ ){ sarray.index[i-1] = p[i+1];};
		sarray.length = seq.length;
	}

	if (doPs) {
		sarray.PrintSuffices(seq.seq, seq.length, printSuffixLength);
	}

	if (doBlt) {
		sarray.BuildLookupTable(seq.seq, seq.length, ltPrefixLength);
	}
	
	//	slowsarray.
	VectorIndex i;
	if (print) {
		//		cout << "order of suffix array: " << endl;
		for (i = 0; i < seq.length; i++) {
			cout << i << " " << sarray.index[i] << endl;
		}
	}
	if (doQuery) {
		FASTAReader queryReader;
		cout << "opening:" << queryFileName << endl;
		queryReader.Init(queryFileName);
		FASTASequence query;
		SAIndex low, high;
		int qIndex = 0;
		while (queryReader.GetNext(query)) {
			query.ToTwoBit();
			sarray.Search(seq.seq, query.seq, query.length, low, high);
			cout << qIndex << " " << low << " " << high << endl;
		}
	}

	if (doLCP) {
		FASTAReader queryReader;
		cout << "opening:" << lcpFileName << endl;
		queryReader.Init(lcpFileName);
		FASTASequence query;
		SAIndex low, high;
		int qIndex = 0;
		DNALength lcpLength = 0;
		while (queryReader.GetNext(query)) {
			query.ToTwoBit();
			sarray.SearchLCP(seq.seq, query.seq, query.length, low, high, lcpLength, 10 );
			cout << qIndex << " " << low << " " << high << " " << lcpLength << endl;
		}
	}
}
