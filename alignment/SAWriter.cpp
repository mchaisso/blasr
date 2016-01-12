#include <vector>
#include <string>
#include "../common/datastructures/suffixarray/SuffixArray.h"
#include "../common/FASTASequence.h"
#include "../common/FASTAReader.h"
#include "../common/NucConversion.h"
#include "../common/Types.h"
#include "../common/datastructures/suffixarray/ssort.h"
#include "../common/algorithms/sorting/qsufsort.h"
#include "../common/algorithms/sorting/Karkkainen.h"
#include "../common/cmpseq/CompressedSequence.h"


void PrintUsage() {
	cout << "usage: sawriter saOut fastaIn [fastaIn2 fastaIn3 ...] [-blt p] [-larsson] [-4bit] [-manmy] [-kar]" << endl;
  cout << "   or  sawriter fastaIn  (writes to fastIn.sa)." << endl;
	cout << "       -blt p      Build a lookup table on prefixes of length 'p'. This speeds " << endl
			 << "                   up lookups considerably (more than the LCP table), but misses matches " << endl
			 << "                   less than p when searching." << endl;
	cout << "       -4bit       Read in (one) fasta file as a compressed sequence file." << endl;
	cout << "       -larsson  (default)  Uses the method of Larsson and Sadakane to build the array." << endl;
	cout << "       -mamy      Uses the method of MAnber and MYers to build the array (slower than larsson, " << endl
			 << "                   and produces the same result. This is mainly for double checking"<<endl
			 << "                   the correctness of larsson)." << endl
			 << "       -kark       Use Karkkainen DS3 method for building the suffix array.  This will probably be more "<<endl
			 << "                   slow than larsson, but takes only an extra N/(sqrt 3) extra space." << endl
			 << "       -mafe       (disabled for now!) Use the lightweight construction algorithm from Manzini and Ferragina" << endl
			 << "       -welter     Use lightweight (sort of light) suffix array construction.  This is a bit more slow than" << endl
			 << "                   normal larsson." << endl
			 << "       -welterweight N use a difference cover of size N for building the suffix array.  Valid values are 7,32,64,111, and 2281." << endl;


}

int main(int argc, char* argv[]) {

	if (argc < 2) {
		PrintUsage();
		exit(1);
	}
	int argi = 1;
	string saFile = argv[argi++];
	vector<string> inFiles;
	
	int doBLT = 1;
	int bltPrefixLength = 13;
	int parsingOptions = 0;
	SAType saBuildType = larsson;
	int read4BitCompressed  = 0;
	int diffCoverSize = 0;
	while (argi < argc) {
		if (strlen(argv[argi]) > 0 and
				argv[argi][0] == '-'){ 
			parsingOptions = 1;
		}
		if (!parsingOptions) {
			inFiles.push_back(argv[argi]);
		}
		else {
			if (strcmp(argv[argi], "-blt") == 0) {
				doBLT = 1;
        if (argi < argc - 1) {
          bltPrefixLength = atoi(argv[++argi]);
          if (bltPrefixLength == 0) {
            cout << argv[argi] << " is not a valid lookup table length." << endl;
            exit(1);
          }
        }
        else {
          cout << "Please specify a lookup table length." << endl;
          exit(1);
        }
			}
			else if (strcmp(argv[argi], "-mamy") == 0) {
				saBuildType = manmy;
			}
			else if (strcmp(argv[argi], "-larsson") == 0) {
				saBuildType = larsson;
			}
			else if (strcmp(argv[argi], "-mcilroy") == 0) {
				saBuildType = mcilroy;
			}
			else if (strcmp(argv[argi], "-slow") == 0) {
				saBuildType = slow;
			}
			else if (strcmp(argv[argi], "-kark") == 0) {
				saBuildType = kark;
			}
			else if (strcmp(argv[argi], "-mafe") == 0) {
				saBuildType = mafe;
			}
			else if (strcmp(argv[argi], "-welter") == 0) {
				saBuildType = welter;
			}
			else if (strcmp(argv[argi], "-welterweight") == 0) {
        if (argi < argc-1) {
          diffCoverSize = atoi(argv[++argi]);
        }
        else {
          cout << "Please specify a difference cover size.  Valid values are 7,32,64,111, and 2281.  Larger values use less memory but may be slower." << endl;
          exit(1);
        }
        if ( ! (diffCoverSize == 7 or 
                diffCoverSize == 32 or
                diffCoverSize == 64 or 
                diffCoverSize == 111 or
                diffCoverSize == 2281) ) {
          cout << "The difference cover size must be one of 7,32,64,111, or 2281." << endl;
          cout << "Larger numbers use less space but are more slow." << endl;
          exit(1);
        }
			}
			else if (strcmp(argv[argi], "-4bit") == 0) {
				read4BitCompressed = 1;
			}
			else {
				PrintUsage();
				cout << "ERROR, bad option: " << argv[argi] << endl;
				exit(1);
			}
		}
		++argi;
	}
  
  if (inFiles.size() == 0) {
    //
    // Special use case: the input file is a fasta file.  Write to that file + .sa
    //
    inFiles.push_back(saFile);
    saFile = saFile + ".sa";
  }
  
	VectorIndex inFileIndex;
	FASTASequence seq;
	CompressedSequence<FASTASequence> compSeq;

	if (read4BitCompressed == 0) {
		for (inFileIndex = 0; inFileIndex < inFiles.size(); ++inFileIndex) {
			FASTAReader reader;
			reader.Init(inFiles[inFileIndex]);
			reader.SetSpacePadding(111);
			if (saBuildType == kark) {
				//
				// The Karkkainen sa building method requires a little extra
				// space at the end of the dna sequence so that counting may
				// be done mod 3 without adding extra logic for boundaries.
				//
			}
  
			if (inFileIndex == 0) {
				reader.ReadAllSequencesIntoOne(seq);
				reader.Close();
			}
			else {
				while(reader.ConcatenateNext(seq)) {
					cout << "added " << seq.title << endl;
				}
			}
		}
		seq.ToThreeBit();
		//seq.ToUpper();
	}
	else {
		assert(inFiles.size() == 1);
		cout << "reading compressed sequence." << endl;
		compSeq.Read(inFiles[0]);
		seq.seq = compSeq.seq;
		seq.length = compSeq.length;
		compSeq.RemoveCompressionCounts();
		cout << "done." << endl;
	}

  //
  // For now, do not allow creation of suffix arrays on sequences > 4G.
  //
  if (seq.length >= UINT_MAX) {
    cout << "ERROR, references greater than " << UINT_MAX << " bases are not supported." << endl;
    cout << "Consider breaking the reference into multiple files, running alignment. " << endl;
    cout << "against each file, and merging the result." << endl;
    exit(1);
  }
	vector<int> alphabet;
	
	SuffixArray<Nucleotide, vector<int> >  sa;
	//	sa.InitTwoBitDNAAlphabet(alphabet);
	//	sa.InitAsciiCharDNAAlphabet(alphabet);
  sa.InitThreeBitDNAAlphabet(alphabet);

	if (saBuildType == manmy) {
		sa.MMBuildSuffixArray(seq.seq, seq.length, alphabet);
	}
	else if (saBuildType == mcilroy) {
		sa.index = new SAIndex[seq.length+1];
		DNALength i;
		for (i = 0; i < seq.length; i++) { sa.index[i] = seq.seq[i] + 1;}
		sa.index[seq.length] = 0;
		ssort(sa.index, NULL);
		for (i = 1; i < seq.length+1; i++ ){ sa.index[i-1] = sa.index[i];};
		sa.length = seq.length;
	}
	else if (saBuildType == larsson) {
		sa.LarssonBuildSuffixArray(seq.seq, seq.length, alphabet);
	}
	else if (saBuildType == kark) {
		sa.index = new SAIndex[seq.length];
		seq.ToThreeBit();
		DNALength p;
		for (p = 0; p < seq.length; p++ ){ seq.seq[p]++; }
		KarkkainenBuildSuffixArray<Nucleotide>(seq.seq, sa.index, seq.length, 5);
		sa.length = seq.length;
	}
	else if (saBuildType == mafe) {
		//		sa.MaFeBuildSuffixArray(seq.seq, seq.length);
		
	}
	else if (saBuildType == welter) {
		if (diffCoverSize == 0) {
			sa.LightweightBuildSuffixArray(seq.seq, seq.length);
		}
		else {
			sa.LightweightBuildSuffixArray(seq.seq, seq.length, diffCoverSize);
		}
	}
	if (doBLT) {
		sa.BuildLookupTable(seq.seq, seq.length, bltPrefixLength);
	}
	sa.Write(saFile);

	return 0;

}
