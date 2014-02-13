#include "../common/FASTASequence.h"
#include "../common/FASTAReader.h"
#include "../common/utils.h"
#include <sstream>

using namespace std;
int main(int argc, char* argv[]) {
	string inFileName, outFileName;
	if (argc < 4) {
		cout << "usage: extractseq infile outfile start1 end1 [-chrom chr] [start2 end2 [-chrom chr] ]  [start3 end3]..."<<endl;
		exit(1);
	}
	inFileName  = argv[1];
	outFileName = argv[2];
	FASTASequence seq;
	FASTAReader reader;
  SequenceIndexDatabase<FASTASequence> seqdb;
	reader.Init(inFileName);
	reader.ReadAllSequencesIntoOne(seq, &seqdb);
	ofstream seqout;
	CrucialOpen(outFileName, seqout);
  int argi = 3;
  DNALength start, end;
  DNALength originalStart, originalEnd;
  while (argi < argc) {
		originalStart = start = atoi(argv[argi]);
		originalEnd   = end   = atoi(argv[argi+1]);
    int chrom = -1;
    bool useSeqDB = false;
    if (argi + 4 <= argc) {
      if (strlen(argv[argi+2]) > 0 and argv[argi+2][0] == '-' and strcmp(argv[argi+2], "-chrom") == 0) {
        useSeqDB = true;
        chrom = atoi(argv[argi+3]);
        argi+=2;
        start = seqdb.ChromosomePositionToGenome(chrom, start);
        end   = seqdb.ChromosomePositionToGenome(chrom,end);
      }
    }
		FASTASequence substring;
		if (start > seq.length or end > seq.length) {
			cout << "ERROR! Coordinates for substring are outside the original" << endl;
			exit(1);
		}
		substring.seq = &seq.seq[start];
		substring.length = end - start;
		
		stringstream titlestream;
		titlestream.str("");
    if (useSeqDB == false) {
      titlestream << seq.title;
    }
    else {
      titlestream << seqdb.names[chrom];
    }
    titlestream << "_" << originalStart << "_" << originalEnd;
		substring.CopyTitle(titlestream.str());
		substring.PrintSeq(seqout);
		argi+=2;
	}
}

