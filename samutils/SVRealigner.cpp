#include "algorithms/alignment/readers/sam/SAMReader.h"
#include "datastructures/alignmentset/SAMAlignment.h"
#include "datastructures/alignmentset/SAMToAlignmentCandidateAdapter.h"
#include "FASTAReader.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "algorithms/alignment/printers/StickAlignmentPrinter.h"
#include "algorithms/alignment/ScoreMatrices.h"
#include "algorithms/alignment/SWAlign.h"
#include "algorithms/alignment/DistanceMatrixScoreFunction.h"
#include "CommandLineParser.h"
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include "datastructures/alignment/AlignmentBlock.h"
#include "datastructures/alignment/AlignmentContext.h"
#include "algorithms/alignment/printers/SAMPrinter.h"

using namespace std;
typedef AlignmentCandidate<FASTASequence,FASTASequence> FastaAlignmentCandidate;
bool NotInGap(Block &b1, Block &b2, int maxGap, int minLength=0) {
	return (b2.length >= minLength and 
					b2.qPos - (b1.qPos + b1.length) < maxGap and
					b2.tPos - (b1.tPos + b1.length) < maxGap);
}

int main(int argc, char* argv[]) {
  string samFileName, genomeFileName;
  string outFileName = "/dev/stdout";
  string alignmentsFileName = "";
  int minGapLength = 50;
	int minMatchLength = 20;
  CommandLineParser clp;
	int maxAlignScore = 50;
	int maxRealignLength = 500;
	string trimmedCoordinatesFileName = "";
	vector<string> samHeader;
	map<string, int> fivePrimeTrim, threePrimeTrim, fullTrim;

  clp.RegisterStringOption("genome", &genomeFileName, "Genome.", true);
  clp.RegisterStringOption("sam", &samFileName, "Alignments.", true);
  clp.RegisterPreviousFlagsAsHidden();
  clp.RegisterStringOption("out", &outFileName, 
													 "Output file. Default to stdout", false);
  clp.RegisterIntOption("match", &minMatchLength, "Minimum ", 
												CommandLineParser::NonNegativeInteger, false);
  clp.RegisterIntOption("maxAlignScore", &maxAlignScore, "Maximum score of pairwise alignment to accept change. ", 
												CommandLineParser::NonNegativeInteger, false);

  clp.RegisterIntOption("gap", &minGapLength, "Minimum ", 
												CommandLineParser::NonNegativeInteger, false);
	clp.RegisterStringOption("trimmed", &trimmedCoordinatesFileName, "A file with coordinates of low-quality sequences to remove.");
  clp.ParseCommandLine(argc, argv);

  FASTAReader fastaReader;

	
	if (trimmedCoordinatesFileName != "") {
		ifstream trimmedIn(trimmedCoordinatesFileName.c_str());
		while (trimmedIn) {
			string chrom, category;
			int start, end, length;
			if (!(trimmedIn >> chrom >> start >> end >> length >> category)) {
				break;
			}
			if (category == "5p") {
				fivePrimeTrim[chrom] = end;
			}
			if (category == "3p") {
				threePrimeTrim[chrom] = start;
			}
			if (category == "full"){ 
				fullTrim[chrom] =1;
			}
		}
	}

  fastaReader.Initialize(genomeFileName);
	DistanceMatrixScoreFunction<DNASequence, DNASequence> scoreFn(LocalAlignLowMutationMatrix, 10, 10);


  ostream *outPtr;
  ofstream outFile;
  ofstream alignmentsOut, gapFastaOut;
	CrucialOpen(outFileName, outFile, std::ios::out);

	int samFileIndex;
  ifstream roiIn;
  vector<FASTASequence> references;

  fastaReader.ReadAllSequences(references);
  int i;
  map<string,int> refToIndex;
  for (i = 0; i < references.size(); i++) {
    refToIndex[references[i].title] = i;
  }

	vector<string> header;
	ifstream samFileHeader(samFileName.c_str());
	string headerLine;
	while(samFileHeader) {
		string line;
		getline(samFileHeader, line);
		if (line.size() > 0 and line[0] == '@') {
			header.push_back(line);
		}
		else {
			break;
		}
	}
	samFileHeader.close();


	SAMReader<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> samReader;
	samReader.Initialize(samFileName);
	cerr << samFileName << endl;
	AlignmentSet<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> alignmentSet;
	samReader.ReadHeader(alignmentSet);

	SAMAlignment samAlignment;  

	int index = 0;
	int lineNumber;
	for (lineNumber=0; lineNumber < header.size(); lineNumber++) {
		outFile << header[lineNumber] << endl;
	}
	while (samReader.GetNextAlignment(samAlignment)) {
		if (samAlignment.rName == "*") {
			continue;
		}
		cerr << samAlignment.qName << endl;
		vector<AlignmentCandidate<> > convertedAlignments;
		//
		// Convert the SAM alignment into a block structure.  This lets
		// the match blocks be added and deleted arbitrarily.  A new CIGAR
		// string is newly built from the blocks, so there is no need to
		// keep track of the CIGAR operations while modifying alignments.
		//
		SAMAlignmentsToCandidates(samAlignment, 
															references, refToIndex,
															convertedAlignments, false);

		int i=0;
		for (i = 0; i < convertedAlignments.size(); i++) {
			//
			// Look for a region of high identity.
			//
			vector<int> realignedGapStart, realignedGapEnd;
			vector<AlignmentCandidate<> > replacementAlignments;

			AlignmentCandidate<> &aln =convertedAlignments[i];
			int b = 0;
			while (b < aln.blocks.size() - 1) {
				int blockStart = b;
				while (b < aln.blocks.size() - 1 and 
							 NotInGap(aln.blocks[b], aln.blocks[b+1], minGapLength)) {
					b++;
				}

				int blockEnd = b;
				int gapStart = b;
				//
				// Advance until the next contiguous sequence.
				while (b < aln.blocks.size() - 1 and 
							 NotInGap(aln.blocks[b], aln.blocks[b+1], minGapLength, minMatchLength) == false) {
					b++;
				}

				int b2 = b;
				int gapEnd = b;

				while (b2 < aln.blocks.size() -1 and 
							 NotInGap(aln.blocks[b2], aln.blocks[b2+1], minGapLength)) {
					b2++;
				}
				// Advance to past the next gap.
				int nextBlockEnd = b2;
				while (b2 < aln.blocks.size() - 1 and 
							 NotInGap(aln.blocks[b2], aln.blocks[b2+1], minGapLength, minMatchLength) == false) {
					b2++;
				}
				int nextGapEnd = b2;




				if (gapEnd - gapStart > 1) {

					//
					// There is opportunity to smush the blocks in the middle
					// into a single sequence that is aligned to either the
					// prefix or suffix of the larger gapped sequence.
					//
						
					int qGapSeqStart = aln.blocks[gapStart].qPos + aln.blocks[gapStart].length;
					int qGapSeqEnd   = aln.blocks[gapEnd].qPos;
					int tGapSeqStart = aln.blocks[gapStart].tPos +  aln.blocks[gapStart].length;
					int tGapSeqEnd   = aln.blocks[gapEnd].tPos;
					int tGapLength = tGapSeqEnd - tGapSeqStart;
					int qGapLength = qGapSeqEnd - qGapSeqStart;
						
					DNASequence qGapSeq, tGapSeq;

					qGapSeq.seq = &aln.qAlignedSeq.seq[qGapSeqStart];
					qGapSeq.length = qGapLength;
						
					tGapSeq.seq = &aln.tAlignedSeq.seq[tGapSeqStart];
					tGapSeq.length = tGapLength;
						
					DNASequence queryPrefixSeq, querySuffixSeq, targetPrefixSeq, targetSuffixSeq;
					
					int tPrefixOffset, tSuffixOffset;
					int qPrefixOffset, qSuffixOffset;
					if (qGapLength > tGapLength) {
						//
						// Fit T into Q
						//
						int fitLength;
						fitLength = min(qGapLength, (int)(tGapLength * 1.1));
						int prefixScore=0,suffixScore=0;

						queryPrefixSeq.seq=&qGapSeq.seq[0];
						queryPrefixSeq.length = fitLength;
						qPrefixOffset = 0;

						querySuffixSeq.seq=&qGapSeq.seq[qGapSeq.length - fitLength];
						querySuffixSeq.length = fitLength;
						qSuffixOffset = qGapSeq.length - fitLength;

						targetPrefixSeq.seq = tGapSeq.seq;
						targetPrefixSeq.length = tGapSeq.length;
						tPrefixOffset = 0;
						targetSuffixSeq.seq = tGapSeq.seq;
						targetSuffixSeq.length = tGapSeq.length;
						tSuffixOffset = 0;
					}
					else {
						// Fit Q into T
						int fitLength = min(tGapLength, (int)(qGapLength * 1.1));
						targetPrefixSeq.seq=tGapSeq.seq;
						targetPrefixSeq.length = fitLength;
						tPrefixOffset= 0;
						targetSuffixSeq.seq=&tGapSeq.seq[tGapSeq.length - fitLength];
						targetSuffixSeq.length = fitLength;
						tSuffixOffset= tGapSeq.length - fitLength;
						queryPrefixSeq.seq = qGapSeq.seq;
						queryPrefixSeq.length = qGapSeq.length;
						qPrefixOffset = 0 ;
						querySuffixSeq.seq = qGapSeq.seq;
						querySuffixSeq.length = qGapSeq.length;
						qSuffixOffset = 0;
					}
					int prefixScore=maxAlignScore,suffixScore=maxAlignScore;
					vector<int> scoreMat;
					vector<Arrow> pathMat;
					AlignmentCandidate<> prefixAlignment, suffixAlignment;
					if (queryPrefixSeq.length < maxRealignLength and targetPrefixSeq.length < maxRealignLength) {
						prefixScore = SWAlign(queryPrefixSeq, targetPrefixSeq, scoreMat, pathMat, prefixAlignment, scoreFn, Global);
					}
					if (querySuffixSeq.length < maxRealignLength and targetSuffixSeq.length < maxRealignLength) {
						suffixScore = SWAlign(querySuffixSeq, targetSuffixSeq, scoreMat, pathMat, suffixAlignment, scoreFn, Global);
					}
					StickPrintAlignment(prefixAlignment, queryPrefixSeq, targetPrefixSeq, cerr);
					StickPrintAlignment(suffixAlignment, querySuffixSeq, targetSuffixSeq, cerr);
					prefixAlignment.tAlignedSeq.seq = targetPrefixSeq.seq;
					prefixAlignment.qAlignedSeq.seq = queryPrefixSeq.seq;
					suffixAlignment.tAlignedSeq.seq = targetSuffixSeq.seq;
					suffixAlignment.qAlignedSeq.seq = queryPrefixSeq.seq;
					prefixAlignment.tAlignedSeq.length = targetPrefixSeq.length;
					prefixAlignment.qAlignedSeq.length = queryPrefixSeq.length;
					suffixAlignment.tAlignedSeq.length = targetSuffixSeq.length;
					suffixAlignment.qAlignedSeq.length = queryPrefixSeq.length;
					
					

					if (prefixScore < maxAlignScore or suffixScore < maxAlignScore) {
						realignedGapStart.push_back(gapStart);
						realignedGapEnd.push_back(gapEnd);

						if (prefixScore <= suffixScore ) {
							int pa;
							for (pa = 0; pa < prefixAlignment.blocks.size(); pa++) {
								prefixAlignment.blocks[pa].qPos += qPrefixOffset;
								prefixAlignment.blocks[pa].tPos += tPrefixOffset;
							}
							replacementAlignments.push_back(prefixAlignment);
						}

						else if (suffixScore < prefixScore ) {
							int sa;
							for (sa = 0; sa < suffixAlignment.blocks.size(); sa++) {
								suffixAlignment.blocks[sa].qPos += qSuffixOffset;
								suffixAlignment.blocks[sa].tPos += tSuffixOffset;
							}
							replacementAlignments.push_back(suffixAlignment);
						}
					}
				}
			}

			aln.gaps.clear();
			if (replacementAlignments.size() > 0) {
				vector<int> newReplacementBlockIndex;
				int r=0;
				int n=0;
				for (b = 0; b < aln.blocks.size(); b++) {
					aln.blocks[n] = aln.blocks[b];

					if (r < realignedGapStart.size() and b == realignedGapStart[r]) {
						newReplacementBlockIndex.push_back(n);
						assert(realignedGapEnd[r]-1>b);
						b += realignedGapEnd[r] - realignedGapStart[r]-1;
						r++;
					}
					n++;
				}
				int nDeleted = b - n;
				aln.blocks.resize(n);
				//
				// Now add the replacement blocks.
				//
				vector<Block> newBlocks;
				r = 0;
				for (b = 0; b < aln.blocks.size(); b++ ) {
					newBlocks.push_back(aln.blocks[b]);
					if (r < newReplacementBlockIndex.size() and b == newReplacementBlockIndex[r]) {
						int tOffset = 0;
						int qOffset = 0;

						tOffset = aln.blocks[b].tPos + aln.blocks[b].length;
						qOffset = aln.blocks[b].qPos + aln.blocks[b].length;

						for (n = 0; n < replacementAlignments[r].blocks.size(); n++) {
							replacementAlignments[r].blocks[n].tPos += tOffset;
							replacementAlignments[r].blocks[n].qPos += qOffset;
							newBlocks.push_back(replacementAlignments[r].blocks[n]);
						}
						r++;
					}
				}
				aln.blocks = newBlocks;

			}

			//
			// Now process the first and last block, since this tends to be wrong with -alignContigs.
			//
			int last = aln.blocks.size() - 1;
			bool endsNotInGap = NotInGap(aln.blocks[last-1], aln.blocks[last], minGapLength, minMatchLength);
			while (last > 0 and aln.blocks[last].length < minMatchLength) {
				last--;
			}
			if (threePrimeTrim.find(aln.qName) != threePrimeTrim.end()) {
				int initialLast = last;
				int trimQPos = threePrimeTrim[aln.qName];
				while (last > 0 and aln.blocks[last].qPos > trimQPos) {
					last--;
				}
				cerr << "Used trimmed file to remove 3p: " << initialLast - last << " " << aln.blocks[initialLast].tPos - aln.blocks[last].tPos  << " " << trimQPos << endl; ;
			}
			if (aln.blocks.size() >= 2 and (endsNotInGap == false or last < aln.blocks.size()-1)) {
				aln.blocks.resize(last);
			}
			
			bool startsNotInGap =  NotInGap(aln.blocks[0], aln.blocks[1], minGapLength, minMatchLength);
			int startBlock = 0;
			while (startBlock < aln.blocks.size() and aln.blocks[startBlock].length < minMatchLength) {
				startBlock ++;
			}
			if (startBlock == 0 and startsNotInGap == false) {
				startBlock = 1;
			}
			if (fivePrimeTrim.find(aln.qName) != fivePrimeTrim.end()) {
				int trimStart = fivePrimeTrim[aln.qName];
				int initialStartBlock = startBlock;
				cerr << "5p: " << trimStart << " " 
						 << aln.blocks[startBlock].qPos << " " << aln.blocks[0].qPos << " " 
						 << aln.blocks[startBlock].qPos - aln.blocks[0].qPos << " " << endl;

				while (startBlock < aln.blocks.size() and aln.blocks[startBlock].qPos < trimStart) {
					startBlock++;
				}
				cerr << "Used trimmed file to remove 5p: " << startBlock - initialStartBlock << endl;
			}
				
			if (aln.blocks.size() >= 2 and startBlock > 0) {
				for (b = 0; b < aln.blocks.size()-startBlock; b++) {
					aln.blocks[b] = aln.blocks[b+startBlock];
				}
				cerr << "Trimming " << startBlock << endl;
				aln.blocks.resize(aln.blocks.size()-startBlock);
			}


			SMRTSequence seq;
			seq.seq = (Nucleotide*) samAlignment.seq.c_str();
			seq.length = samAlignment.seq.size();
			seq.CopyTitle(convertedAlignments[i].qName);
			AlignmentContext defaultAlignmentContext;
			defaultAlignmentContext.readGroupId = samAlignment.rg;
			if (defaultAlignmentContext.readGroupId == "") {
				defaultAlignmentContext.readGroupId = "temp";
			}
			SupplementalQVList emptyQVList;			
			SAMOutput::PrintAlignment(convertedAlignments[i],
																seq,
																outFile,
																defaultAlignmentContext,
																emptyQVList,
																SAMOutput::soft);
		}
	}
	outFile.close();
}
