#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <assert.h>

#include "../common/FASTASequence.h"
#include "../common/FASTAReader.h"
#include "../common/DNASequence.h"
#include "../common/SeqUtils.h"
#include "../common/tuples/DNATuple.h"
#include "../common/tuples/TupleList.h"
#include "../common/tuples/TupleMetrics.h"

using namespace std;

int main(int argc, char* argv[]) {
  FASTAReader reader;
  if (argc < 3) {
	cout << "usage: printReadWordCount readsFile tupleFile outputFile" << endl;
	exit(1);
  }
	TupleMetrics tm;
  string readsFileName = argv[1];
  string tupleFileName = argv[2];
  
  TupleList<CountedDNATuple> tupleList;
  tupleList.InitFromFile(tupleFileName);
	tupleList.GetTupleMetrics(tm);
  cout << "read: " << tupleList.GetLength() << " tuples from a list." << endl;
  cout << "tuple size is: " << tm.tupleSize << endl;
  //
  // Open the reads file handle 
  //
  reader.Init(readsFileName);

  FASTASequence seq;
  
  while (reader.GetNext(seq)) {
	int p;
	// For now don't try and handle 'N''s.
	if (!OnlyACTG(seq))
	  continue;
	CountedDNATuple tuple;
	cout << ">" << seq.title << endl;
	for (p = 0; p < seq.length - tm.tupleSize + 1; p++) { 
	  assert(tuple.FromStringRL(&(seq.seq[p]), tm));
	  int tupleIndex = tupleList.Find(tuple);
	  if (tupleIndex >= 0){ 
		cout.width(3);
		cout << tupleList[tupleIndex].count << ",";
	  }
	  else cout << "  0,";
	  if ((p+1) % 20 == 0)
		cout << endl;
	}
	cout << endl;
  }
  return 0;
}

