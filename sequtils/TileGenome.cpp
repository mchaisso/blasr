#include "FASTAReader.h"
#include "FASTASequence.h"
#include "utils.h"

#include <string>
#include <sstream>
#include <algorithm>
using namespace std;
int main(int argc, char* argv[]) {
  string inFileName, outFileName;
  int length;
  if (argc < 4) {
    cout <<"usage: tileGenome input.fasta output.fasta seq_length [-stide S|-coverage C]" << endl
         << "  -stride S samples a read every S bases." << endl
         << "  -coverage C samples a read every seq_length / c bases." << endl;
    exit(1);
  }
  inFileName = argv[1];
  outFileName = argv[2];
  length = atoi(argv[3]);
  if (length == 0) {
    cout <<"ERROR, length must be greater than 0" << endl;
    exit(1);
  }
  int argi = 4;
  int stride = 0;
  float coverage = 0;

  while (argi < argc) {
    if (strcmp(argv[argi], "-stride") == 0) {
      stride = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-coverage")== 0) {
      coverage = atof(argv[++argi]);
    }
    ++argi;
  }
  
  FASTAReader reader;
  reader.Initialize(inFileName);
  ofstream outFile;
  CrucialOpen(outFileName, outFile, std::ios::out);

  FASTASequence genome;
  while (reader.GetNext(genome)) {
		if (stride == 0 and coverage == 0) {
			cout << "ERROR, must provide stride or coverage. " << endl;
			exit(1);
		}
		if (stride == 0) {
			stride = length / coverage;
		}
		
		int pos = 0;
		while (pos < genome.length) {
			int seqLength = min(genome.length - pos, (unsigned int)length);
			FASTASequence seq;
			seq.seq = &genome.seq[pos];
			seq.length = seqLength;
			stringstream titleStrm;
			titleStrm << genome.GetName() << "_" << pos << "_" << pos + seqLength;
			seq.CopyTitle(titleStrm.str());
			seq.PrintSeq(outFile);
			delete[] seq.title;
			pos += stride;
		}
	}
  outFile.close();
}

  
  
