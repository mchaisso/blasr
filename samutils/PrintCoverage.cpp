#include "algorithms/alignment/readers/sam/SAMReader.h"
#include "datastructures/alignmentset/SAMAlignment.h"
#include "datastructures/alignmentset/SAMToAlignmentCandidateAdapter.h"
#include "FASTAReader.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "algorithms/alignment/printers/StickAlignmentPrinter.h"
#include "CommandLineParser.h"
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <pthread.h>

void Increment(unsigned int &value, int increment = 1) {
	value += increment;
	/*	int bigval = value;
	if (increment + bigval >= 255) {
		value = 255;
	}
	else {
		value+=increment;
		}*/
}

bool FindReadIndex(string &read, int &start, int &end) {
	start = read.find('/');
	if (start != read.npos) {
		end = read.find('/', start+1);
		if (end != read.npos) {
			return true;
		}
	}
	return false;
}

bool CompareReadNumbers(string &read1, string &read2) {
	int r1s, r1e, r2s, r2e;
	if (FindReadIndex(read1, r1s, r1e) and 
			FindReadIndex(read2, r2s, r2e)) {
		string n1= read1.substr(r1s, r1e-r1s);
		string n2 =read2.substr(r2s, r2e-r2s);
		return  n1 == n2;
	}
	else {
		return false;
	}
}

void WriteConsensus(string fileName, vector< unsigned int > &counts, int binSize) {
	ofstream outFile;
	outFile.open(fileName.c_str(), std::ios::out|std::ios::binary);
	outFile.write((const char*) &binSize, sizeof(int));
	int length = counts.size();
	outFile.write((const char*) &length, sizeof(int));
	int i;

	outFile.write((const char*) &counts[0], sizeof(unsigned int)*length);
	outFile.close();
}

class Data {
public:
	string samFileName;
	vector<vector<unsigned int> > *refCounts;
	bool unique;
  map<string,int> *refToIndex;
	int binSize;
};

void StoreLengths(Data* data) {
	string samFileName = (data->samFileName);
	cerr << samFileName << endl;
	vector<vector<unsigned int > >& refCounts = *data->refCounts;
	bool unique = data->unique;
  map<string,int> &refToIndex = *data->refToIndex;
	int binSize = data->binSize;
	SAMReader<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> samReader;

	samReader.Initialize(samFileName);
	AlignmentSet<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> alignmentSet;
	samReader.ReadHeader(alignmentSet);
			
	SAMAlignment samAlignment;  

	int index = 0;
	bool sameAsPrev = false;
	string prevRead= "";
	while (samReader.GetNextAlignment(samAlignment)) {

			
		vector<AlignmentCandidate<> > convertedAlignments;
				
		int c;
		//
		// Only count primary alignments.
		//
		if (samAlignment.flag & 256 != 0) {
			continue;
		}
		if (samAlignment.mapQV < 10) {
			continue;
		}
		//
		// Check to see if coverage of subreads is considered.
		//
		if (! unique or  // if not unique, all subreads count.
				samAlignment.qName == prevRead  or // process multiple disjoint hits from the same read.
				(unique and CompareReadNumbers(samAlignment.qName, prevRead) == false)) {
					

			int refIndex = refToIndex[samAlignment.rName];

			DNALength tStart, tEnd;
			tStart = samAlignment.pos;
			tEnd   = samAlignment.pos + samAlignment.tLen;
        
			DNALength tPos = tStart - 1; // SAM is 1 based.

			vector<int> opLengths;
			vector<char> ops;
			samAlignment.cigar.Vectorize(opLengths, ops);
			DNALength middle = (tEnd - tStart) / 2 + tStart;
			for (tPos = tStart - 1; tPos < middle; tPos++) {
				// possible race condition, but not worth locking.
				refCounts[refIndex][tPos / binSize] +=1;
			}
			for (tPos = middle; tPos < tEnd - 1; tPos++) {
				// possible race condition, but not worth locking.
				refCounts[refIndex][tPos / binSize] += 1;
			}
		}
		prevRead = samAlignment.qName;
	}
}


int main(int argc, char* argv[]) {

	CommandLineParser clp;
	vector<string> samFileNames;
	string outDir;
	int binSize = 10;
	string referenceFileName = "";
	bool unique = false;
	int maxThreads = 1;
  clp.RegisterStringListOption("sam", &samFileNames, "Alignments.");
	clp.RegisterStringOption("outDir", &outDir, "Make this directory for output.");
	clp.RegisterIntOption("bin", &binSize, "Count by this bin size.", CommandLineParser::PositiveInteger);
	clp.RegisterFlagOption("unique", &unique, "Consider only one subread from a read.");
	clp.RegisterIntOption("j", &maxThreads, "Max threads to launch.", CommandLineParser::PositiveInteger);
  clp.ParseCommandLine(argc, argv);
	
	
	mkdir(outDir.c_str(), S_IRWXU |
				S_IRGRP | S_IXGRP |
				S_IROTH | S_IXOTH);
	
	bool setupRef = true;
	int s;
	vector<string> refNames;
	vector<vector<unsigned int > > refCounts;

  int i;
  map<string,int> refToIndex;
	long nBases = 0;
	s = 0;
	while (s < samFileNames.size()) {


		cerr << samFileNames[s] << endl;
		if (setupRef) {

			SAMReader<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> samReader;
			string samFileName = samFileNames[s];
			samReader.Initialize(samFileName);
			AlignmentSet<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> alignmentSet;
			samReader.ReadHeader(alignmentSet);
		
			string prevRead = "";
			setupRef = false;
			int ri;
			refCounts.resize(alignmentSet.references.size());
			for (ri = 0; ri < alignmentSet.references.size(); ri++) {
				refNames.push_back(alignmentSet.references[ri].sequenceName);
				refCounts[ri].resize(alignmentSet.references[ri].length / binSize + alignmentSet.references[ri].length % binSize, 0); // base
				refToIndex[alignmentSet.references[ri].sequenceName] = ri;
			}
		}
		int ti;
		int nThreads = min((int) samFileNames.size() - s, maxThreads);
		pthread_t *threads = new pthread_t[nThreads];		
		pthread_attr_t *threadAttr = new pthread_attr_t[nThreads];
		vector<Data> threadData(nThreads);
		int procIndex;
		for (procIndex = 0; procIndex < nThreads; procIndex++ ){
			pthread_attr_init(&threadAttr[procIndex]);
			threadData[procIndex].samFileName = samFileNames[s];
			s++;
			threadData[procIndex].unique = unique;
			threadData[procIndex].refCounts = &refCounts;
			threadData[procIndex].refToIndex = &refToIndex;
			threadData[procIndex].binSize = binSize;
			pthread_create(&threads[procIndex], &threadAttr[procIndex], (void* (*)(void*))StoreLengths, &threadData[procIndex]);
			
		}

		for (procIndex = 0; procIndex < nThreads; procIndex++ ){
			pthread_join(threads[procIndex], NULL);
		}

	}
	
	int r;
	for (r = 0; r < refCounts.size(); r++) {
		string outFileName = outDir + "/" + refNames[r] + ".data";
		int p;
		for (p = 0; p < refCounts[r].size(); p++) {
			refCounts[r][p] /= binSize;
		}
		WriteConsensus(outFileName, refCounts[r], binSize);
	}
	return 0;
}
