#include "../common/algorithms/alignment/AffineKBandAlign.h"
#include "../common/algorithms/alignment/SWAlign.h"
#include "../common/algorithms/alignment/AlignmentUtils.h"
#include "../common/algorithms/alignment/AlignmentPrinter.h"
#include "../common/algorithms/alignment/DistanceMatrixScoreFunction.h"
#include "../common/datastructures/alignment/Path.h"
#include "../common/datastructures/alignment/AlignmentCandidate.h"
#include "../common/algorithms/alignment/ScoreMatrices.h"
#include "../common/SMRTSequence.h"
#include "../common/FASTQSequence.h"
#include "../common/FASTAReader.h"
#include "../common/FASTASequence.h"
#include "../common/tuples/DNATuple.h"
#include "../common/tuples/TupleMetrics.h"
#include "../common/CommandLineParser.h"
#include "../common/files/ReaderAgglomerate.h"
#include "../common/files/WriterAgglomerate.h"
#include "../common/utils/FileOfFileNames.h"
#include "../common/SMRTSequence.h"
#include <string>
#include <set>

using namespace std;

class SeqKeyword {
public:
	DNATuple tuple;
	int readPos;
	int operator<(const SeqKeyword &rhs) const {
		if (tuple == rhs.tuple) {
			return readPos < rhs.readPos;
		}
		else {
			return tuple < rhs.tuple;
		}
	}
	SeqKeyword & operator=(const SeqKeyword &rhs) {
		tuple     = rhs.tuple;
		readPos   = rhs.readPos;
		return *this;
	}
};

class CompareKeywordsByTupleOnly {
public:
	int operator()(const SeqKeyword &lhs, const SeqKeyword &rhs) const {
		return lhs.tuple < rhs.tuple;
	}
};

typedef AlignmentCandidate<FASTQSequence, FASTQSequence> FastqAlignment;

int main(int argc, char* argv[]) {
	
	string adaptersFileName, readsFileName, outputFileName;
	// seed alignment matches for speeding local alignment.
	TupleMetrics tm;
	tm.tupleSize = 4;
	float maxInsRate = 0.2;
	CommandLineParser clp;
	bool useComplement        = false;
	bool splitReadsAtAdaptors = false;
	int  gapInit = 1;
	int  gapExt  = 1;
	int  maxScore = -21;
	int  match = -1;
	clp.SetProgramName("ourmatch");
	clp.SetProgramSummary("Count the number of occurrences of every k-mer in a file.");
	clp.RegisterStringOption("reads",  &readsFileName,  "The reads to screen.");
	clp.RegisterStringOption("adapters", &adaptersFileName, "The file of adapters to screen for.");
	clp.RegisterStringOption("readsout", &outputFileName, "The screened reads.");
	clp.RegisterPreviousFlagsAsHidden();
	clp.RegisterIntOption("gap_init", &gapInit, "Gap initialization penalty.", CommandLineParser::NonNegativeInteger);
	clp.RegisterIntOption("gap_ext", &gapExt, "Gap extension penalty. For now, non-affine alignment is used gap_ext is ignored.", 
												CommandLineParser::NonNegativeInteger);
	clp.RegisterIntOption("match", &match, "Match bonus.", CommandLineParser::Integer);
	clp.RegisterIntOption("maxScore", &maxScore, "Maximum alignment score.", CommandLineParser::Integer);
	clp.RegisterIntOption("wordsize", &tm.tupleSize, "Size of keyword to speed alignment", 
												CommandLineParser::NonNegativeInteger);
	clp.RegisterFloatOption("maxInsRate", &maxInsRate, "Maximum insertion rate",
													CommandLineParser::NonNegativeFloat);
	clp.RegisterFlagOption("comp", &useComplement, "Allow adapter to match in the complementary strand.");
	clp.RegisterFlagOption("split", &splitReadsAtAdaptors, "Split reads that have adaptors.");
	vector<string> leftovers;
	clp.ParseCommandLine(argc, argv, leftovers);
	
	int i;

	for(i =0;i<5;i++) { CrossMatchMatrix[i][i] = match;}

	//
	// Process the reads into a vector of read keywords
	//
	vector<string> readsFileNames;	
	ReaderAgglomerate reader;
	ReaderAgglomerate adapterReader;
	WriterAgglomerate<FASTQSequence> writer;

  bool optionIsFOFN;
  if (FileOfFileNames::IsFOFN(readsFileName)) {
    FileOfFileNames::FOFNToList(readsFileName, readsFileNames);
  }
  else {
    readsFileNames.push_back(readsFileName);
  }
 
	reader.Initialize(readsFileNames[0]);
	adapterReader.Initialize(adaptersFileName);
	writer.Initialize(outputFileName);
	if (reader.GetFileType() == HDFBase or reader.GetFileType() == HDFPulse) {
		if (writer.GetFileType() == HDFBase or writer.GetFileType() == HDFPulse) 
			writer.hdfWriter.AddRunInfo(reader.hdfBasReader.GetMovieName(), reader.hdfBasReader.GetRunCode());
	}

		
	DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distScoreFn;
	distScoreFn.del = 5;
	distScoreFn.ins = 5;
	distScoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);

	ofstream maskedReadsOut;
	// Determine the output file type.

	CrucialOpen(outputFileName, maskedReadsOut, std::ios::out);
	//
	// Read all adapters into a table.
	//

	FASTASequence adapter;
	FASTQSequence seqRC;

	int readIndex = 0;
	vector<FASTASequence> adaptors;
	vector<vector<SeqKeyword> >  adaptorKeywords;
	SeqKeyword keyword;
	vector<SeqKeyword> emptySeqKeywordVect;
	set<int> foundAdapters;
	while(adapterReader.GetNext(adapter)) {
		DNALength pos;
		if (adapter.length < tm.tupleSize) 
			continue;
		adaptors.push_back(adapter);
		adaptorKeywords.push_back(vector<SeqKeyword>() );
		for (pos = 0; pos < adapter.length - tm.tupleSize + 1; pos++) {
			keyword.tuple.FromStringLR(&adapter.seq[pos], tm);
			keyword.readPos = pos;
			adaptorKeywords[readIndex].push_back(keyword);
		}
		readIndex++;
		if (useComplement) {
			adaptorKeywords.push_back(vector<SeqKeyword>() );
			adapter.MakeRC(seqRC);
			adaptors.push_back(seqRC);
			for (pos = 0; pos < seqRC.length - tm.tupleSize + 1; pos++) {
				keyword.tuple.FromStringLR(&seqRC.seq[pos], tm);
				keyword.readPos = pos;
				adaptorKeywords[readIndex].push_back(keyword);
			}
			readIndex++;
		}
	}

	//
	// Sort the index into each adaptor.
	//
	int adaptorIndex;
	for (adaptorIndex = 0; adaptorIndex < adaptorKeywords.size(); adaptorIndex++) {
		std::sort(adaptorKeywords[adaptorIndex].begin(),
							adaptorKeywords[adaptorIndex].end());
	}
	
  std::vector<int> prevAlignedGenomePos;
  std::vector<int> readOptScore;
  std::vector<FastqAlignment > optAlignment;
	std::vector<int> optGenomeAlignPos;
	std::vector<int> optGenomeAlignLength;

	SMRTSequence read;
	vector<bool> alignedAt;
	std::vector<SeqKeyword>::iterator keyIt, endIt;
	int numAdaptors = adaptorKeywords.size();
	vector<int> scoreMat;
	vector<Arrow> pathMat;
	Alignment alignment;
	readIndex = 0;
	int nAdapters = 0;
	int nScreenedReads = 0;
	bool adapterFound;
	DistanceMatrixScoreFunction<FASTASequence, FASTASequence> scoreFn;
	scoreFn.InitializeScoreMatrix(CrossMatchMatrix);
	int numPrinted = 0;
  int readsFileIndex = 0;
  for (readsFileIndex = 0; readsFileIndex < readsFileNames.size(); readsFileIndex++) {
	reader.Initialize(readsFileNames[readsFileIndex]);
	while (reader.GetNext(read)) {
		adapterFound = false;
		alignedAt.resize(read.length);
		int readPos;
		int readAlignStart, readAlignEnd, readAlignLength;
		CompareKeywordsByTupleOnly comp;
		vector<int> maskStart, maskEnd;
		for (adaptorIndex = 0; adaptorIndex < numAdaptors; adaptorIndex++) {
			std::fill(alignedAt.begin(), alignedAt.end(), false);
			if (tm.tupleSize < read.length) {
				for (readPos = 0; readPos < read.length - tm.tupleSize + 1; readPos++) {
					assert(readPos >= 0);
					assert(readPos <= ((int)read.length) - tm.tupleSize);
					keyword.tuple.FromStringLR(&read.seq[readPos], tm);
					keyIt    = lower_bound(adaptorKeywords[adaptorIndex].begin(), adaptorKeywords[adaptorIndex].end(), keyword, comp);
					endIt    = upper_bound(adaptorKeywords[adaptorIndex].begin(), adaptorKeywords[adaptorIndex].end(), keyword, comp);
					for (; keyIt != endIt; keyIt++) {
						//					cout << " rp: " << readPos << " key read pos: " << (*keyIt).readPos << endl;
						readAlignStart = readPos - (*keyIt).readPos;
						int fixedReadAlignStart = readAlignStart;
						//
						// If the anchor is way too far at the beginning of the
						// sequence, just skip this one.
						//
						if (readAlignStart < 0)
							continue;
						readAlignEnd   = readAlignStart + adaptors[adaptorIndex].length;
						//
						// Adjust the allowed boundaries according 
						// to indel rates.
						//
						int insWiggle =  maxInsRate * adaptors[adaptorIndex].length;

						// adjust front and end wiggle according to the position of the anchor.
						// sequence, just skip this one.
						//
						if (readAlignStart < 0)
							continue;
						readAlignEnd   = readAlignStart + adaptors[adaptorIndex].length;
						//
						// Adjust the allowed boundaries according 
						// to indel rates.
						//
						int frontInsWiggle, endInsWiggle;

						// adjust front and end wiggle according to the position of the anchor.
						frontInsWiggle = insWiggle * ((1.0*(*keyIt).readPos) / adaptors[adaptorIndex].length);
						endInsWiggle   = insWiggle - frontInsWiggle;
					
						if (readAlignStart - frontInsWiggle < 0) {
							frontInsWiggle = readAlignStart;
						}
					
						readAlignStart -= frontInsWiggle;

						//
						// Do not repeat the same alignments.
						//
						if (alignedAt[readAlignStart] == true) {
							continue;
						}

						alignedAt[readAlignStart] = true;
					
						// adjust the end wiggle
						if (readAlignEnd + endInsWiggle > read.length) {
							endInsWiggle = read.length - readAlignEnd;
						}
						readAlignEnd += endInsWiggle;
						readAlignLength = readAlignEnd - readAlignStart;
					
						FASTASequence readSubstring;
						assert(readAlignStart >= 0);
						if (readAlignLength > read.length) continue;
						readSubstring.title = read.title;
						readSubstring.titleLength = read.titleLength;
						readSubstring.seq = &read.seq[readAlignStart];
						readSubstring.length = readAlignLength;
						int alignScore;
						alignment.qPos = 0;
						alignment.tPos = 0;
						/*
						alignScore =  KBandAlign(adaptors[adaptorIndex], readSubstring,
																		 CrossMatchMatrix, 1,1, 10,
																		 scoreMat, pathMat, 
																		 alignment, scoreFn, LocalBoundaries);
						*/
						alignScore = SWAlign(adaptors[adaptorIndex], readSubstring,
																 scoreMat, pathMat, 
																 alignment, distScoreFn, LocalBoundaries);

						alignment.tPos += readAlignStart;

						if (alignScore < maxScore) {
							maskStart.push_back(alignment.tPos);
							maskEnd.push_back(alignment.tPos + alignment.tLength);
							readPos += adaptors[adaptorIndex].length + 1;
							nAdapters++;
							adapterFound = true;
							foundAdapters.insert(adaptorIndex);
						}
						//
						// If the alignment score is sufficiently good, record this
						// as a position to be masked.
						//
					}
				}
			}
		}

		if (adapterFound) {
			++nScreenedReads;
		}
		int maskIndex, maskPos;
		int curReadPos = 0;
		FASTQSequence substr;

		// 
		// Merge overlapping mask positions.
		//
		int m, n;
		if (maskStart.size() > 0) {
			for (m = 0; m < maskStart.size()- 1; m++ ) {
				for (n = m + 1; n < maskStart.size(); ){ 
					if ((maskStart[m] >= maskStart[n] and maskStart[m] <= maskEnd[n]) or
							(maskEnd[m] >= maskStart[n] and maskEnd[m] <= maskEnd[n])) {
						if (maskStart[m] > maskStart[n]) {
							maskStart[m] = maskStart[n];
						}
						if (maskEnd[m] < maskEnd[n]) {
							maskEnd[m] = maskEnd[n];
						}
						maskStart.erase(maskStart.begin() + n);
						maskEnd.erase(maskEnd.begin() + n);
					}
					else {
						n++;
					}
				}
			}
		}

		for (maskIndex =0 ; maskIndex < maskStart.size(); maskIndex++) {
			
			if (splitReadsAtAdaptors) {
				substr.ReferenceSubstring(read, curReadPos, maskStart[maskIndex] - curReadPos);
				if (substr.length > 0) {
					//  Assign the title of the substring.
					// If this is the first pass, just print that out.
          // if this is not, the x and y coordinates need to be modified
          // so that they reflect the fact that these are reads from the
          // same template, just different passes around the bell.
          //
          // Form a specific title.
          //
          // parse an ASTRO title, and yank out the x and y coordinates
          // >x1_y2_2000088-0001_m091023_075136_Twe_p1_b20
          DNALength t;
          int numUnderscores = 0;
          int firstUnderscore = 0;
          for (t = 0; t < read.titleLength; t++ ){
            if (read.title[t] == '_') {
              numUnderscores++;
					    if (numUnderscores == 1) {
								firstUnderscore = t;
							}
						}
						if (numUnderscores == 2) { break; }
					}
					if (t == read.titleLength) {
						cout << "ERROR, the read titles must be in ASTRO format right now. This "<< endl
								 << "one doesn't seem to be: " << read.title << endl;
						exit(1);
					}
					string titleSuffix(&read.title[t], read.titleLength - t + 1);
					int xpos = atoi(&read.title[1]);
					if (xpos == 0){ 
						cout << "ERROR, the read titles either indicate a read has a ZMW coordinate of 0, or" << endl
								 << "this " << read.title << " was not parsed correctly. Both are not ok." << endl;
						exit(1);
					}
					xpos += 1000 * (maskIndex + 1);
					int ypos = atoi(&read.title[firstUnderscore+2]);
					stringstream newtitle;
					newtitle << "x" << xpos << "_y" << ypos << titleSuffix;
					substr.CopyTitle(newtitle.str());
					curReadPos = maskEnd[maskIndex]+1;
					cout << "Writing substr from " << read.title << endl;
					//					writer.Write(substr);
				}
			}
			else {
				for (maskPos = maskStart[maskIndex]; maskPos < maskEnd[maskIndex]; maskPos++) {
					read.seq[maskPos] = 'X';
				}
			}
		}
	
		if (splitReadsAtAdaptors) {
			// Write out the last substring (or the entire read if there are
			// no adapaters in it.
			if (curReadPos < read.length)  {
				substr.ReferenceSubstring(read, curReadPos, read.length - curReadPos);
				if (substr.length > 0) {
					substr.CopyTitle(read.title);
					cout << "writing seq for " << read.title << endl;
					writer.Write(substr);
					++numPrinted;
				}
			}
		}
		else {
			//
			// Not splitting the reads.  Just write out the masked read.
			//
			writer.Write(read);
		}
		++readIndex;
	}
}
	cout << "screened " << nAdapters << " from " << nScreenedReads << " ";
	set<int>::iterator setit, setend;
	setend = foundAdapters.end();
	setit = foundAdapters.begin();
	while (setit != setend) {
		cout << *setit << " ";
		++setit;
	}
	cout << endl;
	
	return 0;
}
