#include <assert.h>
#include <string>
#include <iostream>
#include <vector>

#include "FASTAReader.h"
#include "FASTASequence.h"
#include "algorithms/alignment.h"
#include "CommandLineParser.h"
using namespace std;


void SplitRead(FASTASequence &read, int readStart, int readLength, 
							 FASTASequence &ad1, FASTASequence &ad2, 
							 int indelCost,
							 vector<int> &passStarts,
							 vector<int> &passLengths,
							 vector<int> &lastAdapters,
							 int lastAdapter,
							 vector<int> &scoreMat,
							 vector<Arrow> &pathMat, 
							 float minIdentity, int minLength) {
	
	/*
	 *
	 * Align each adapter to the reads. Pick the higher scoring one, 
	 * and divide the read into two halves.  If there is no adapter found
	 * in either half, add it to the list of split reads.
	 *
	 */
	int ad1Score, ad2Score;
	FASTASequence readSubseq;
	readSubseq.seq = &read.seq[readStart];
	readSubseq.length = readLength;
	DistanceMatrixScoreFunction<FASTASequence, FASTASequence> distScoreFn;
	distScoreFn.del = indelCost;
	distScoreFn.ins = indelCost;
	distScoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);

	MatchedAlignment ad1Alignment, ad2Alignment;
	ad1Score = SWAlign(ad1, readSubseq, scoreMat, pathMat, 
										 ad1Alignment, distScoreFn, QueryFit);
	ad1Alignment.tStart = ad1Alignment.qStart = 0;
	cout << "adapter 1: " << endl;
	StickPrintAlignment(ad1Alignment, ad1, readSubseq, cout);
	ComputeAlignmentStats(ad1Alignment, ad1.seq, readSubseq.seq, distScoreFn); 

	cout << endl;
	cout << "adapter 2: " << endl;
	ad2Score = SWAlign(ad2, readSubseq, scoreMat, pathMat, 
										 ad2Alignment, distScoreFn, QueryFit);
	ad2Alignment.tStart = ad2Alignment.qStart = 0;
	StickPrintAlignment(ad2Alignment, ad2, readSubseq, cout);
	ComputeAlignmentStats(ad2Alignment, ad2.seq, readSubseq.seq, distScoreFn);
	cout << ad1Alignment.pctSimilarity << " " << ad2Alignment.pctSimilarity << endl;
	int adStart, adEnd;
	adStart = adEnd = -1;
	if (ad1Alignment.pctSimilarity > ad2Alignment.pctSimilarity) {
		if (ad1Alignment.blocks.size() > 0) {
	  	if (ad1Alignment.pctSimilarity / 100 > minIdentity) {
				cout << " splitting on ad1." << endl;
				adStart = ad1Alignment.tPos;
				int lastBlock = ad1Alignment.blocks.size() - 1;
				adEnd   = ad1Alignment.tPos + ad1Alignment.blocks[lastBlock].tPos + ad1Alignment.blocks[lastBlock].length;
				lastAdapter = 1;
			}
		}
	}
	else {
		if (ad2Alignment.pctSimilarity / 100 > minIdentity) {
			if (ad2Alignment.blocks.size() > 0) {
				cout << " splitting on ad2." << endl;
				int lastBlock = ad2Alignment.blocks.size() - 1;
				adStart = ad2Alignment.tPos;
				adEnd   = ad2Alignment.tPos + ad2Alignment.blocks[lastBlock].tPos + ad2Alignment.blocks[lastBlock].length;
				lastAdapter = 2;
			}
		}
	}
 
	if (adStart >= 0 and adEnd >= 0) {
		int firstHalfStart, firstHalfLength, secondHalfStart, secondHalfLength;
		firstHalfStart = readStart;
		firstHalfLength = adStart - readStart;
		secondHalfStart = adEnd+ 1;
		secondHalfLength = readLength - secondHalfStart;
		cout << "split coords: " << firstHalfStart << " " << firstHalfLength << " " << secondHalfStart <<" " << secondHalfLength << endl;
		if (firstHalfLength > minLength) {
			SplitRead(read, firstHalfStart, firstHalfLength, ad1, ad2, indelCost, passStarts, passLengths, lastAdapters,
								lastAdapter,
								scoreMat, pathMat, minIdentity, minLength);
		}
		
		if (secondHalfLength > minLength) {
			SplitRead(read, secondHalfStart, secondHalfLength, ad1, ad2, indelCost, passStarts, passLengths, lastAdapters, 
								lastAdapter,
								scoreMat, pathMat, minIdentity, minLength);
		}
	}
	else {
		passStarts.push_back(readStart);
		passLengths.push_back(readLength);
		lastAdapters.push_back(lastAdapter);
	}
}


int main(int argc, char* argv[]) {
	string ad1File, ad2File, readsFile, readsOutFile;

	FASTAReader ad1Reader;
	FASTAReader ad2Reader;
	FASTAReader reader;


	CommandLineParser cl;
	float minPctSimilarity = 0.60;
	int indel = 3;
	int minLength = 10;
	cl.RegisterStringOption("ad1", &ad1File, "FASTA file with the first adapter");
	cl.RegisterStringOption("ad2", &ad2File, "FASTA file with the second adapter");
	cl.RegisterStringOption("reads", &readsFile, "FASTA file with SMRTBell reads");
	cl.RegisterStringOption("readsout", &readsOutFile, "output file for split reads");
	cl.RegisterPreviousFlagsAsHidden();
	cl.RegisterFloatOption("pctSim", &minPctSimilarity, "Minimum percent similarity to trigger a match to an adapter.", 
												 CommandLineParser::PositiveFloat);
	cl.RegisterIntOption("indel", &indel, "Penalty for indel (positive)", CommandLineParser::NonNegativeInteger);
	cl.RegisterIntOption("minLength", &minLength, "Minimum length pass to retain.", CommandLineParser::PositiveInteger);
	vector<string> opts;
	cl.ParseCommandLine(argc, argv, opts);

	/*
	 * Open all the required files, quitting if they are unavailable.
	 */

	ad1Reader.Init(ad1File);
	ad2Reader.Init(ad2File);
	reader.Init(readsFile);

	ofstream splitOut;
	CrucialOpen(readsOutFile, splitOut);

	FASTASequence ad1, ad2;
	ad1Reader.GetNext(ad1);
	ad2Reader.GetNext(ad2);

	FASTASequence read;
	vector<int> scoreMat;
	vector<Arrow> pathMat;
	int readIndex = 0;
	while(reader.GetNext(read)) {
		read.ToUpper();
		//
		// Do a fitting sequence alignment to match one of the two 
		// adapters into the read.
		//
		vector<int> passStarts, passLengths, la;
		read.PrintSeq(cout);
		SplitRead(read, 0, read.length, ad1, ad2,
							indel, 
							passStarts, passLengths,la, 0,
							scoreMat, pathMat, minPctSimilarity, minLength);
		int i;
		for (i = 0; i < passStarts.size(); i++) {
			cout << "read: " << readIndex << " pass: " << i << " " << passStarts[i] << " " << passLengths[i] << " " << la[i] << endl;
		}
		++readIndex;
	}
}

	
