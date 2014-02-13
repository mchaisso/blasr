#include "datastructures/sequence/PackedDNASequence.h"
#include "FASTAReader.h"
#include <string>
#include <iostream>
using namespace std;

int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << "usage: testPackedDNA dnaseq" << endl;
		exit(0);
	}
	
	string dnaseqName = argv[1];
	FASTAReader reader; 
  reader.Init(dnaseqName);
  FASTASequence seq;
	reader.GetNext(seq);
  PackedDNASequence packedSeq;
	
	packedSeq.CreateFromDNASequence(seq);
	
	int nA, nC, nT, nG;
	nA = packedSeq.CountNuc(0,seq.length, 'A');
	nC = packedSeq.CountNuc(0,seq.length, 'C');
	nG = packedSeq.CountNuc(0,seq.length, 'G');
	nT = packedSeq.CountNuc(0,seq.length, 'T');

  cout << " nA: " << nA << " nC " << nC << " nG " << nG << " nT " << nT << endl;
	nG = packedSeq.CountNuc(1,seq.length, 'G');
	cout << "ng 1, " << seq.length << ": " << nG << endl;
	nG = packedSeq.CountNuc(1,seq.length-5, 'G');
	cout << "ng 1, " << seq.length-5 << ": " << nG << endl;
  return 0;
}
