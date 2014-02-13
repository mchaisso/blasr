#include <string>

#include "../common/FASTASequence.h"
#include "../common/FASTAReader.h"
#include "../common/SeqUtils.h"
#include "../common/qvs/QualityValue.h"
#include "../common/cmpseq/CompressedSequence.h"
#include "../common/utils.h"


void PrintUsage() {
		cout << "Usage: compressseq inFile outFile [-revindex] [-4bit] [-indexStride S]" << endl;
		cout << " -revindex  Creates an index to map from compressed sequence to original. " << endl << endl
				 << " -indexStride S Specify how frequently to log an index point (G=1000)." <<endl << endl
				 << " -4bit       Create a reversibly compressed sequence.  Otherwise, all homopolymers." <<endl
				 << "             will be condensed into single-bases and it is necessary to have the " << endl
				 << "             full original genome to align reads." << endl;
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		PrintUsage();
		exit(1);
	}
	string inFileName, outFileName;
	inFileName = argv[1];
	outFileName = argv[2];
	int argi = 3;
	int doBuildIndex = 0;
	int do4BitCompression   = 0;
	string indexName = "";
	int doCondense   = 1;
	int indexStride  = 1000;
	while (argi < argc) {
		if (strcmp(argv[argi], "-revindex") == 0 ) {
			doBuildIndex = 1;
		}
		else if (strcmp(argv[argi], "-4bit")  == 0) {
			do4BitCompression = 1;
			doCondense = 0;
		}
		else if (strcmp(argv[argi], "-stride")  == 0) {
			indexStride = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			cout << "Invalid option: " << argv[argi] << endl;
			exit(1);
		}
		++argi;
	}
	
	FASTAReader reader;
	reader.Init(inFileName);
	ofstream outFile;
					 
	
	CompressedSequence<FASTASequence> seq;
	reader.GetNext(seq);
	seq.SetHasTitle();
	
	if (doBuildIndex) {
		int run = 0;
		if (do4BitCompression) {
			run = 15;
		}
		seq.ToUpper();
		seq.SetHasIndex();
		seq.BuildReverseIndex(run, indexStride);
	}

	if (doCondense) {
		FASTASequence fastaRead;
		seq.CondenseHomopolymers();
		fastaRead = seq;
		ofstream outFile;
		CrucialOpen(outFileName, outFile);
		fastaRead.PrintSeq(outFile);
	}

	if (do4BitCompression) {
		seq.FourBitCompressHomopolymers();
		seq.Write(outFileName);
	}
	return 0;
}
