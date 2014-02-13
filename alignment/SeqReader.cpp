#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>

#include "../common/FASTASequence.h"
#include "../common/FASTAReader.h"
#include "../common/DNASequence.h"
#include "../common/tuples/DNATuple.h"

using namespace std;

int main(int argc, char* argv[]) {
  FASTAReader reader;
  if (argc < 2) {
	cout << "usage: seqReader fastaFile (simply read in a fasta file for timing)" << endl;
	exit(1);
  }

  string fileName = argv[1];

  reader.Init(fileName);

  FASTASequence seq;
  int numRead = 0;
  while (reader.GetNext(seq)) { ++numRead;}

  cout << "read: " << numRead << " sequences." << endl;


  return 0;
}

