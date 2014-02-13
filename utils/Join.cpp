#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include "utils.h"
#include "utils/StringUtils.h"
#include <assert.h>

using namespace std;

// This is a sort of an implementation of a join of two text files based on keywords in each text file.

class Line {
public:
	vector<string> values;
	Line& operator=(Line&rhs) {
		values = rhs.values;
		return *this;
	}
	friend ostream& operator<<(ostream &out, Line &l) {
		int vIndex;
		for (vIndex = 0; vIndex < l.values.size(); vIndex++) {
			out << l.values[vIndex] << " ";
		}
		return out;
	}
};

class SortLines {
public:
	int k;
	int operator()(const Line &a, const Line&b) const {
		return a.values[k].compare(b.values[k]) < 0;
	}
};


void ReadFile(ifstream &in, vector<Line> &lines) {
	while (in) {
		string linestring;
		std::getline(in,linestring);
		Line newline;
		stringstream valuestream(linestring);
		while(valuestream) {
			string val;
			valuestream >> val;
			if (val.size() > 0) {
				newline.values.push_back(val);
			}
		}
		int nValues = newline.values.size();
		if (nValues > 0) 
			lines.push_back(newline);
		assert(lines[lines.size() -1].values.size() > 0);
	}
}
		

int main(int argc, char* argv[]) {
	
	string aFileName, bFileName, aOutFileName, bOutFileName;

	if (argc < 5) {
		cout << "usage: join afile bfile aout bout [-ka a_key_index] [-kb b_key_index]"<<endl;
		return 1;
	}
	aFileName = argv[1];
	bFileName = argv[2];
	aOutFileName = argv[3];
	bOutFileName = argv[4];

	int akey = 0;
	int bkey = 0;
	
	int argi = 5;
	while (argi < argc) {
		if (strcmp(argv[argi], "-ka") == 0) {
			akey = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-kb") == 0) {
			bkey = atoi(argv[++argi]);
		}
		argi++;
	}
	
	SortLines asorter, bsorter;
	asorter.k = akey;
	bsorter.k = bkey;

	ifstream aFile, bFile;

	CrucialOpen(aFileName, aFile);
	CrucialOpen(bFileName, bFile);

	ofstream aOutFile, bOutFile;
  CrucialOpen(aOutFileName, aOutFile, std::ios::out);
  CrucialOpen(bOutFileName, bOutFile, std::ios::out);

	vector<Line> aLines, bLines;
	ReadFile(aFile, aLines);
	ReadFile(bFile, bLines);
	std::sort(aLines.begin(), aLines.end(), asorter);
	std::sort(bLines.begin(), bLines.end(), bsorter);

	int aIndex, bIndex;
	aIndex = 0; bIndex = 0;
	int nALines, nBLines;
	nALines = aLines.size();
	nBLines = bLines.size();
	while (aIndex < nALines and bIndex < nBLines) {
		int cmp = aLines[aIndex].values[akey].compare(bLines[bIndex].values[bkey]);
		if (cmp == 0) {
			aOutFile << aLines[aIndex] << endl;
			bOutFile << bLines[bIndex] << endl;
			++aIndex;
			++bIndex;
		}
		else if (cmp < 0) {
			++aIndex;
		}
		else {
			++bIndex;
		}
	}
	return 0;
}
	
