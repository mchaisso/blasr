#include <iostream>
#include <fstream>
#include "datastructures/bwt/BWT.h"
using namespace std;


int main(int argc, char* argv[]) {

	if (argc < 2) {
		cout << "usage: printpbbwt in.bwt" << endl;
		exit(1);
	}

	string bwtFileName = argv[1];

 	Bwt<PackedDNASequence, FASTASequence> bwt;
	bwt.Read(bwtFileName);
  bwt.bwtSequence.PrintUnpacked(cout);
  return 0;
}
