#include "../common/datastructures/suffixarray/SuffixArray.h"
#include "../common/datastructures/suffixarray/SuffixArrayTypes.h"
#include "../common/FASTASequence.h"
#include "../common/FASTAReader.h"
#include "../common/NucConversion.h"


void PrintUsage() {
	cout << "usage:  saquery reference_file ref_suffixarray_file queryfile [options]" << endl;
	cout << "  -print   Print all locations of hits." << endl << endl;
	cout << "  -max N   Lookup locations if there are less than N hits." << endl << endl
			 << "  -count   Only count the number of times a pattern appears, rather than finding locations." << endl << endl;
}

using namespace std;
int main(int argc, char* argv[]) {

	if (argc < 4) {
		PrintUsage();
		exit(1);
	}
	int argi = 1;
	string refFileName = argv[argi++];
	string saFile = argv[argi++];
	string queryFileName = argv[argi++];
	int printPos = 0;
	int count =0;
	bool countOnly = false;
  bool printCount = false;
	while (argi < argc) {
		if (strcmp(argv[argi], "-print") == 0) {
			printPos = 1;
		}
    else if (strcmp(argv[argi], "-printCount") == 0) {
      printCount = true;
    }
		else if (strcmp(argv[argi], "-max") == 0) {
			count = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-count") == 0) {
			countOnly = true;
		}
		else {
			cout << "ERROR, bad option: " << argv[argi] << endl;
			exit(1);
		}
		++argi;
	}
	FASTAReader refReader;
	VectorIndex inFileIndex;
	FASTASequence reference;

	refReader.Initialize(refFileName);
	refReader.ReadAllSequencesIntoOne(reference);

	DNASuffixArray sa;
	sa.Read(saFile);
	FASTAReader queryReader;
	queryReader.Initialize(queryFileName);
	FASTASequence query;
	SAIndex low, high;
	int qIndex = 0;
	DNALength lcpLength;
	vector<DNALength> positions;
	vector<DNALength> lowMatchBound, highMatchBound;
	
	while (queryReader.GetNext(query)) {
		query.ToUpper();
		positions.clear();
		lowMatchBound.clear();
		highMatchBound.clear();
		sa.StoreLCPBounds(&reference.seq[0], (long int) reference.length, 
											&query.seq[0], query.length,
											true, 0,
											lowMatchBound, highMatchBound);
		DNALength matchIndex;
		DNALength saLow = lowMatchBound[lowMatchBound.size()-1];
		DNALength saHigh = highMatchBound[highMatchBound.size()-1];
    if (printCount and saHigh - saLow < count) {
      cout << saHigh - saLow << endl;
    }
		if (countOnly == false and (count == 0 or saHigh - saLow < count)) {
			
			for (matchIndex = saLow; matchIndex <= saHigh; matchIndex++) {
				positions.push_back(sa.index[matchIndex]);
				if (printPos) {
					cout << sa.index[matchIndex] << " ";
				}
			}
		}
		//		cout << "matched " << positions.size()  << " positions." << endl;
		if (printPos) {
			cout << endl;
		}
	}
}

	
	
	

