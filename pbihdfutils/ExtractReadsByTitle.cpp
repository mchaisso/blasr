#include "data/hdf/HDFPlsReader.h"
#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFBasWriter.h"
#include "data/hdf/PlatformId.h"
#include "utils/StringUtils.h"
#include "utils/FileOfFileNames.h"
#include "files/WriterAgglomerate.h"

#include "utils.h"
#include "SMRTSequence.h"
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include "datastructures/reads/BaseFile.h"
using namespace std;



class Read {
public: 
	int x, y;
	string movieName;
	void InitFromReadTitle(string title) {
		// Parse read titles in the format:
		// x10_y22_0971105-0001_m100508_081612_Cog_p1_b10/0_1181
		stringstream titlestrm(title);
		char c;
		titlestrm.get(c);
		titlestrm >> x;
		titlestrm.get(c);
		titlestrm.get(c);
		titlestrm >> y;
		titlestrm.get(c);
		string titleSuffix;
		titlestrm >> titleSuffix;
		int p, p2;
		// advance to '_' preceeding movie name.
		for(p = 0; p < titleSuffix.size() && titleSuffix[p] != '_'; p++);
		p++;
		movieName.assign(titleSuffix, p, titleSuffix.size() - p);
	}
};

class OrderByMovie {
public:
	int operator()(const Read &lhs, const Read &rhs) {
		return lhs.movieName < rhs.movieName;
	}
};

class OrderByMovieAndCoordinate {
public:
	int operator()(const Read &lhs, const Read &rhs) {
		if (lhs.movieName == rhs.movieName) {
			if (lhs.x == rhs.x) {
				return lhs.y < rhs.y;
			}
			else {
				return lhs.x < rhs.x;
			}
		}
		else {
			return lhs.movieName < rhs.movieName;
		}
	}
};


int main(int argc, char* argv[]) {
	
	string plsFofnName, readNamesFileName, outFileName;

	if (argc < 4) {
		cout << "usage: extractReadsByTitle in.{pls.h5|fofn} titles.txt out.pls.h5|out.fasta " << endl;
		cout << " in.{pls.h5|fofn} is either a .pls/.bas file or file of file names." << endl
				 << " titles.txt is a file of read names to extrace. " << endl
				 << " out.pls.h5 is where they all go." << endl;
		exit(1);
	}
	plsFofnName = argv[1];
	readNamesFileName = argv[2];
	outFileName = argv[3];
	
	vector<string> plsFileNames;
	
	if (FileOfFileNames::IsFOFN(plsFofnName)) {
		FileOfFileNames::FOFNToList(plsFofnName, plsFileNames);
	}
	else {
		plsFileNames.push_back(plsFofnName);
	}

	std::vector<string> readNames;
	ifstream readNamesFile;
	CrucialOpen(readNamesFileName, readNamesFile);
	while(readNamesFile) {
		string readName;
		getline(readNamesFile, readName);
		if (readName.size() > 0) {
			readNames.push_back(readName);
		}
	}
	
	vector<Read> reads;
	reads.resize(readNames.size());
	int i;
	for (i = 0; i < reads.size(); i++) {
		reads[i].InitFromReadTitle(readNames[i]);
	}

	sort(reads.begin(), reads.end(), OrderByMovieAndCoordinate());

	//
	// Now process reads in hdf files and output the ones in the read title list.
	//
	
	FileType fileType;
	BaseSequenceIO::DetermineFileTypeByExtension(outFileName, fileType);
	
	
	HDFBasWriter writer;
	ofstream seqOut;
	if (fileType == HDFBase) {
		writer.Initialize(outFileName, "ex_movie", "ex_runcode");
	}
	else if (fileType == Fasta || fileType == Fastq) {
		CrucialOpen(outFileName, seqOut);
	}
	int f;
	FASTQSequence seq;
	Read seqInfo;
	
	for (f = 0; f < plsFileNames.size(); f++) {
		HDFBasReader reader;
		reader.Initialize(plsFileNames[f]);
		int numExtracted = 0;
		while(reader.GetNext(seq)) {
			Read seqInfo;
			seqInfo.InitFromReadTitle(seq.title);
			vector<Read>::iterator searchIt;
			searchIt = lower_bound(reads.begin(), reads.end(), seqInfo, OrderByMovieAndCoordinate());
			int searchIndex = searchIt - reads.begin();
			if (searchIt != reads.end() and 
					searchIt->movieName == seqInfo.movieName and
					searchIt->x == seqInfo.x and
					searchIt->y == seqInfo.y) {
				if (fileType == HDFBase) {
					writer.Write(seq);
				}
				else {
					seq.PrintSeq(seqOut);
				}
				++numExtracted;
			}
		}
		cout << numExtracted <<" " << f << " " << plsFileNames[f] << endl;
	}
}
				 
