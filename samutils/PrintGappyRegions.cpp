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

bool GetReadIndex(string &read, int &readIndex) {
	int start, end;
	if (FindReadIndex(read, start, end)) {
		readIndex = atoi(&read.c_str()[start]);
		return true;
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



class Event {
public:
	int refIndex;
	int position;
	int type;
	int readIndex;
	int operator<(const Event &rhs) const { 
		if (type != rhs.type) {
			return type < rhs.type;
		}
		else {
			if (refIndex != rhs.refIndex) {
				return refIndex < rhs.refIndex;
			}
			else {
				if (position != rhs.position) {
					return position < rhs.position;
				}
				else {
					return readIndex < rhs.readIndex;
				}
			}
		}
	}
};

class Data {
public:
	string samFileName;
	bool unique;
  map<string,int> *refToIndex;
	int binSize;
	vector<Event> *events;
	int readIndexOffset;
	float minInsRatio, minDelRatio;
};


void StoreConsensus(Data* data) {
	string samFileName = (data->samFileName);
	cerr << samFileName << endl;
	vector<Event> &events = *data->events;
	int readIndexOffset=data->readIndexOffset;
	float minInsRatio = data->minInsRatio;
	float minDelRatio = data->minDelRatio;
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
		if (unique and samAlignment.mapQV < 10) {
			continue;
		}

		int readIndex;
		GetReadIndex(samAlignment.rName, readIndex);
		
		//
		// Check to see if coverage of subreads is considered.
		//
		vector<vector<unsigned int> > counts;
		counts.resize(4);
		int nBins = samAlignment.tLen / binSize + samAlignment.tLen % binSize;
	  for (c = 0; c < 4; c++) {
			counts[c].resize(nBins, 0);
		}
		
		if (! unique or  // if not unique, all subreads count.
				samAlignment.qName == prevRead  or // process multiple disjoint hits from the same read.
				(unique and CompareReadNumbers(samAlignment.qName, prevRead) == false)) {
					

			int refIndex = refToIndex[samAlignment.rName];

			DNALength tStart, tEnd;
			tStart = samAlignment.pos;
			tEnd   = samAlignment.pos + samAlignment.tLen;
        
			DNALength tPos = 0;

			vector<int> opLengths;
			vector<char> ops;
			samAlignment.cigar.Vectorize(opLengths, ops);

			// Compact ins or del stretches.
			vector<int> compOpLengths;
			vector<int> compOps;
			for (c = 0; c < ops.size(); c++) {
				if (ops[c] == 'M') {
					int d;
					d = c + 1;
					int nMatch = 0;
					int nGap = opLengths[d];
					while (d < ops.size() - 3 and ops[d+2] == ops[d] and ops[d+1] == 'M' and opLengths[d+1] < 5) {
						opLengths[c]  += opLengths[d+1];
						opLengths[c+1] += opLengths[d+2];
						d += 2;
					}
					if (d > c + 1) {
						compOpLengths.push_back(opLengths[c]);
						compOps.push_back(ops[c]);
						compOpLengths.push_back(opLengths[c+1]);
						compOps.push_back(ops[c+1]);
						c = d - 1;
					}
					else {
						compOps.push_back(ops[c]);
						compOpLengths.push_back(opLengths[c]);
					}
				}
			}
			//			ops = compOps;
			cout << "THIS CODE IS NOT DONE"<< endl;
			assert(0);
			opLengths = compOpLengths;
			
			for (c = 0; c < ops.size(); c++) {
				
				//
				// Although an alignment may begin with an insertion or
				// deletion, they are not counted here, just advance past any
				// gaps in the beginning.
				//
				int gi = 0;

				if (ops[c] == 'M' or ops[c] == 'D') {
					int i;
					int whichArray = 0;
					if (ops[c] == 'D') whichArray = 1;
					for (i = 0; i < opLengths[c]; i++) {
						assert(tPos < alignmentSet.references[refIndex].length);
						Increment(counts[whichArray][tPos / binSize], 1);
						++tPos;
					}
				}
				else if (ops[c] == 'I') {
					int i;
					assert(tPos < alignmentSet.references[refIndex].length);
					Increment(counts[2][tPos / binSize], opLengths[c]);
				}
			}
			
			int i, j;
			for (i = 0; i < nBins; i++) {
				float total = 0;
				for (j = 0; j < 3; j++) {
					total += counts[j][i];
				}
				Event e;
				e.refIndex = refIndex;
				e.position = (tPos + i*binSize) / binSize;
				e.readIndex = readIndexOffset + readIndex;
				if (counts[1][i] / total > minInsRatio) {
					e.type = 1;
					events.push_back(e);
				}
				else if (counts[2][i] / total > minDelRatio) {
					e.type = 2;
					events.push_back(e);
				}
			}
		}
		prevRead = samAlignment.qName;
	}
}


int main(int argc, char* argv[]) {

	CommandLineParser clp;
	vector<string> samFileNames;
	int binSize = 10;
	string referenceFileName = "";
	bool unique = false;
	int maxThreads = 1;
	float minInsRatio = 0.5;
	float minDelRatio = 0.5;

  clp.RegisterStringListOption("sam", &samFileNames, "Alignments.");
	clp.RegisterIntOption("bin", &binSize, "Count by this bin size.", CommandLineParser::PositiveInteger);
	clp.RegisterFlagOption("unique", &unique, "Consider only one subread from a read.");
	clp.RegisterIntOption("j", &maxThreads, "Max threads to launch.", CommandLineParser::PositiveInteger);
	clp.RegisterFloatOption("i", &minInsRatio, "Minimum insertion ratio to store.", CommandLineParser::PositiveFloat);
	clp.RegisterFloatOption("d", &minDelRatio, "Minimum deletion ratio to store.", CommandLineParser::PositiveFloat);
  clp.ParseCommandLine(argc, argv);
	
	
	
	bool setupRef = true;
	int s;
	vector<string> refNames;
	vector<vector<vector<unsigned int> > > refCounts;

  int i;
  map<string,int> refToIndex;
	long nBases = 0;
	s = 0;
	vector<Event> events;
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
				refCounts[ri].resize(3);
				refCounts[ri][0].resize(alignmentSet.references[ri].length / binSize + alignmentSet.references[ri].length % binSize, 0); // base
				refCounts[ri][1].resize(alignmentSet.references[ri].length / binSize + alignmentSet.references[ri].length % binSize, 0); // insertion
				refCounts[ri][2].resize(alignmentSet.references[ri].length / binSize + alignmentSet.references[ri].length % binSize, 0); // deletion
				refToIndex[alignmentSet.references[ri].sequenceName] = ri;
			}
		}


		int ti;
		int nThreads = min((int) samFileNames.size() - s, maxThreads);
		pthread_t *threads = new pthread_t[nThreads];		
		pthread_attr_t *threadAttr = new pthread_attr_t[nThreads];
		vector<Data> threadData(nThreads);
		int procIndex;
		vector<vector<Event> > tEvents(nThreads);

		int offset;
		for (procIndex = 0; procIndex < nThreads; procIndex++ ){
			pthread_attr_init(&threadAttr[procIndex]);
			threadData[procIndex].samFileName = samFileNames[s];
			s++;
			threadData[procIndex].unique = unique;
			threadData[procIndex].events = &tEvents[procIndex];
			threadData[procIndex].minInsRatio = minInsRatio;
			threadData[procIndex].minDelRatio = minDelRatio;
			threadData[procIndex].refToIndex = &refToIndex;
			threadData[procIndex].binSize = binSize;


			pthread_create(&threads[procIndex], &threadAttr[procIndex], (void* (*)(void*))StoreConsensus, &threadData[procIndex]);
			
		}
		// Append the results
		for (procIndex = 0; procIndex < nThreads; procIndex++ ){
			pthread_join(threads[procIndex], NULL);
			events.insert(events.end(), tEvents[procIndex].begin(), tEvents[procIndex].end());			
			cout << events.size() << endl;
		}
	}

	sort(events.begin(), events.end());
	return 0;
}
