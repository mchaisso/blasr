#include <string>
#include "../common/datastructures/suffixarray/SharedSuffixArray.h"
#include "../common/datastructures/suffixarray/SuffixArrayTypes.h"
#include "../common/DNASequence.h"
#include "../common/FASTAReader.h"

int main(int argc, char* argv[]) {

	bool ps = false;
	string refName;
	string saFile;
	if (argc < 2) {
		cout << "usage: saprinter file.sa [-ps genome.fa]" << endl;
		cout << " -ps prints 50 bases of suffices" << endl;
		exit(1);
	}
	saFile = argv[1];
	int argi = 2;
	while (argi < argc) {
		if (strcmp(argv[argi], "-ps") == 0) {
			refName = argv[++argi];
			ps = true;
		}
		++argi;
	}

	DNASuffixArray sarray;
	sarray.Read(saFile);
	FASTASequence genome;
	if (ps){ 
		FASTAReader reader;
		reader.Initialize(refName);
		reader.GetNext(genome);
	}
		
	DNASequence seq;
	unsigned int i;
	cout << sarray.length << endl;
	for (i = 0; i < sarray.length; i++ ){
		cout << sarray.index[i];
		if (ps) {
			if (sarray.length - sarray.index[i] < 50) {
				seq.length = sarray.length - sarray.index[i];
			}
			else {
				seq.length = 50;
			}
			seq.seq = &genome.seq[sarray.index[i]];
			cout << " ";
			seq.PrintSeq(cout);
		}
		else {
			cout << endl;
		}
	}

	return 0;
}
