#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>

#include "../common/FASTASequence.h"
#include "../common/FASTAReader.h"
#include "../common/DNASequence.h"
#include "../common/tuples/DNATuple.h"
#include "../common/tuples/TupleMetrics.h"
#include "../common/Types.h"

using namespace std;

int main(int argc, char* argv[]) {
  FASTAReader reader;
  if (argc < 5) {
	cout << "usage: wordCounter seqFile tupleSize tupleOutputFile posOutputFile" << endl;
	exit(1);
  }

  string fileName = argv[1];
  int    tupleSize = atoi(argv[2]);
  string tupleListName = argv[3];
	string posOutName    = argv[4];
  
	TupleMetrics tm;
  tm.Initialize(tupleSize);
  reader.Init(fileName);

  FASTASequence seq;
  reader.GetNext(seq);

  vector<CountedDNATuple> tupleList;
  CountedDNATuple tuple;
  DNALength i;
  for (i = 0; i < seq.length - tm.tupleSize + 1; i++ ) {
		if (tuple.FromStringRL((Nucleotide*) (seq.seq + i), tm)) {
			tuple.count = i;
			tupleList.push_back(tuple);
		}
  }

  std::sort(tupleList.begin(), tupleList.end());

  int t;
  int t2;
  int numTuples = tupleList.size();
  t = t2 = 0;
  int numUnique = 0;
  while (t < numTuples) {
	t2 = t;
	t2++;
	while (t2 < numTuples and tupleList[t] == tupleList[t2]) {
	  t2++;
	}
	++numUnique;
	t = t2;
  }

  ofstream countedTupleListOut;
  countedTupleListOut.open(tupleListName.c_str(), ios_base::binary);

	ofstream posOut;
	posOut.open(posOutName.c_str(), ios_base::binary);

  countedTupleListOut.write((const char*) &numUnique, sizeof(int));
  countedTupleListOut.write((const char*) &tm.tupleSize, sizeof(int));

  posOut.write((const char*) &numUnique, sizeof(int));

	//
	// Write out the tuple+counts to a file.
	//
  t = t2 = 0;
  CountedDNATuple countedTuple;
	int numMultOne = 0;
  while (t < numTuples) {
		t2 = t;
		t2++;
		while (t2 < numTuples and tupleList[t] == tupleList[t2]) {
			t2++;
		}
		countedTuple.tuple = tupleList[t].tuple;
		countedTuple.count = t2 - t;
		if (countedTuple.count == 1) ++numMultOne;
		countedTupleListOut.write((const char*) &countedTuple,sizeof(CountedDNATuple));
		
		posOut.write((char*)&countedTuple.count, sizeof(int));
		
		int tc;
		for (tc = t; tc < t2; tc++) {
			posOut.write((char*) &tupleList[tc].count, sizeof(int));
		}
		t = t2;
  }

	//
	// Write out the positions of the tuples to a file.
	//
	
	posOut.close();
	countedTupleListOut.close();

	//  cout << "found " << numUnique << " distinct " << DNATuple::TupleSize << "-mers." << endl;
	cout << numMultOne << endl;
  return 0;
}

