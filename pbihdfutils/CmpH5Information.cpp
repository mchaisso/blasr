#include "data/hdf/HDFCmpFile.h"
#include "datastructures/alignment/CmpFile.h"
#include "CommandLineParser.h"
#include "datastructures/alignment/ByteAlignment.h"
#include "algorithms/alignment/printers/StickAlignmentPrinter.h"
#include "PrintAlignment.h"
#include <string>
#include <map>


using namespace std;
using namespace __gnu_cxx;

void AppendBinnedErrorRate(vector<unsigned char> &byteAlignment, vector<vector<float> > &bins) {
	int b, bi;
	int binWidth = byteAlignment.size() / bins.size();
	int pos;
	for (b = 0; b < bins.size(); b++) { 
		pos = b * binWidth;
		int nMatch = 0;
		for (bi = 0; bi < binWidth; bi++) {
			if (IsMatch(byteAlignment, pos+bi)) {
				nMatch++;
			}
		}
		bins[b].push_back(((1.0)*nMatch)/binWidth);
	}
}

void AddDistancesBetweenErrors(ByteAlignment &byteAlignment, map<int,int> &hist) {
	int c, lastError;
	lastError = 0;
	for (c = 0; c < byteAlignment.size(); c++) {
		if (IsMatch(byteAlignment, c) == false) {
			int d = c - lastError;
			if (hist.find(d) == hist.end()) {
				hist[d] = 1;
			}
			else {
				hist[d]++;
			}
			lastError = c;
		}
	}
}

void StoreMatchCounts(vector<unsigned char>& byteAlignment, 
                      vector<int> &matchVect,
                      vector<int> &insVect,
                      vector<int> &delVect,
                      vector<int> &mismatchVect) {


  if (byteAlignment.size() > matchVect.size()) {
    cout << "resizing to " << byteAlignment.size() << endl;
    int curVectSize = matchVect.size();
    int s= byteAlignment.size();
    matchVect.resize(s);
    insVect.resize(s);
    delVect.resize(s);
    mismatchVect.resize(s);
    std::fill(&matchVect[curVectSize], &matchVect[matchVect.size()], 0);
    std::fill(&insVect[curVectSize], &insVect[insVect.size()], 0);
    std::fill(&delVect[curVectSize], &delVect[delVect.size()], 0);
    std::fill(&mismatchVect[curVectSize], &mismatchVect[mismatchVect.size()], 0);
  }
  int i;
  for (i = 0; i < byteAlignment.size(); i++) {
    int q, t;
    q = QueryChar[byteAlignment[i]];
    t = RefChar[byteAlignment[i]];
    if (q == ' ') {
      delVect[i]++;
    }
    else if (t == ' ') {
      insVect[i]++;
    }
    else if (q != t) {
      mismatchVect[i]++;
    }
    else {
      matchVect[i]++;
    }
  }
}

void StoreMovingAverage(vector<unsigned char> &byteAlignment, vector<float> &movingAverage, int window=100, int stride=25) {
  int i;
  int nMatch = 0;
  if (byteAlignment.size() < window) {
    movingAverage.resize(0);
    return;
  }
  else {
    movingAverage.resize(byteAlignment.size() / window);
  }
  for (i = 0; i < window and i < byteAlignment.size(); i++) {
    
    if (QueryChar[byteAlignment[i]] == RefChar[byteAlignment[i]]) {
      nMatch++;
    }
  }
  movingAverage[0] = ( (1.0)*nMatch/window);
  for (i = window; i < byteAlignment.size(); i++) {
    // Move past the old window
    if (QueryChar[byteAlignment[i-window]] == RefChar[byteAlignment[i-window]]) {
      nMatch--;
    }
    // Advance to current
    if (QueryChar[byteAlignment[i]] == RefChar[byteAlignment[i]]) {
      nMatch++;
    }
    movingAverage[i-window+1] = ( (1.0)*nMatch/window);
  }
}

void StoreErrorRateMatrix(vector<unsigned char> &byteAlignment,  vector<vector<int>  > &cor,   vector<vector<int>  > &incor) {
  int lengthBin = min(byteAlignment.size()/1000, cor.size()-1);
  int i;
  if (byteAlignment.size() > cor[0].size()) {
    int prevLen = cor[0].size();
    for (i = 0; i < cor.size(); i++) {
      cor[i].resize(byteAlignment.size());
      incor[i].resize(byteAlignment.size());
      fill(&cor[i][prevLen], &cor[i][cor[i].size()], 0);
      fill(&incor[i][prevLen], &incor[i][incor[i].size()], 0);
    }
  }
  for (i = 0; i < byteAlignment.size(); i++) {
    if (QueryChar[byteAlignment[i]] == RefChar[byteAlignment[i]]) {
      cor[lengthBin][i]++;
    }
    else {
      incor[lengthBin][i]++;
    }
  }
}

void StoreBinnedErrorRate(vector<unsigned char> &byteAlignment, vector<float> &bins) {
	int b, bi;
	int binWidth = byteAlignment.size() / bins.size();
	int pos;
	int templateLength = CountBasesInQuery(byteAlignment);
	// int type ceil.
	int basesPerBin = templateLength / bins.size() + (templateLength % bins.size() ? 1 : 0);
	
	int binStart, binEnd;
	binStart = 0;
	int binIndex = 0;
	while (binStart < byteAlignment.size()) {
		binEnd = binStart;
		int binBaseCount = 0;
		int nBinIns = 0, nBinDel = 0, nBinMismatch = 0, nBinMatch = 0;
		while (binEnd < byteAlignment.size() and binBaseCount < basesPerBin) {
			if (QueryChar[byteAlignment[binEnd]] != ' ') { binBaseCount++; }
			binEnd++;
		}
		CountStats(byteAlignment, nBinMatch, nBinMismatch, nBinDel, nBinIns, binStart, binEnd);
		bins[binIndex] += 1.0 - (1.0*nBinIns + nBinDel + nBinMismatch)/binBaseCount;
		binStart = binEnd;
		++binIndex;
	}
}

int main(int argc, char* argv[]) {
	string cmpFileName;

	CommandLineParser clp;
	bool printTotalAlignedBases = false;
	bool printNReads = false;
	bool printGlobalAccuracy = false;
	bool printNMatches = false;
	bool printAverageAccuracy = false;
	bool printDistBetweenErrorsHist = false;
	float identityCutoff = 0.0;
	bool printBinnedErrorRate = false;
	int minAlignLength = 0;
	bool printBreakdown = false;
	int nBins = 20;
	float totalPercentIdentity = 0;
	string binnedErrorDistributionFileName = "";
  bool countBasesByMovie = false;
  string matchGapFileName = "";
  string readMatchFileName = "";
  int matchGapK = 15;
  string matchRunFileName = "";
  string lengthsFileName = "";
  bool printAverageLength = false;
  bool printMovingAverage = false;
  string movingAverageFileName = "";
  string matchCountFileName = "";
  string errorMatrixName = "";
  bool discardFirstMatch = false;
  bool discardLastMatch  = false;
	clp.RegisterStringOption("cmpH5File", &cmpFileName, "Input cmp.h5 file.", true);
	clp.RegisterPreviousFlagsAsHidden();
	clp.RegisterFlagOption("nAlignments",            &printNReads, "Print the total number of reads.", false);
	clp.RegisterFlagOption("nAlignedBases", &printTotalAlignedBases, "Print the total number of aligned bases", false);
	clp.RegisterFloatOption("identityCutoff",   &identityCutoff, "Print the total number of aligned bases", CommandLineParser::PositiveFloat, false);
	clp.RegisterFlagOption("binnedErrorRate", &printBinnedErrorRate, "Divide each read into N bins, and count the error rate in each bin", false);
	clp.RegisterIntOption("nBins", &nBins, "The number of bins for binned error rate.", CommandLineParser::PositiveInteger, false);
	clp.RegisterIntOption("minAlignLength", &minAlignLength, "Disregard alignments less than n length.", CommandLineParser::PositiveInteger, false);
  clp.RegisterStringOption("matchCount", &matchCountFileName, "Print the count of matches/mismatches/ins/del per pos to file.", false);
	clp.RegisterFlagOption("globalAccuracy", &printGlobalAccuracy, "Print the accuracy across all sequences.", false);
	clp.RegisterFlagOption("averageAccuracy", &printAverageAccuracy, "Print average accuracy of reads.", false);
	clp.RegisterFlagOption("nMatches", &printNMatches, "Print the number of bases matched.", false);
	clp.RegisterFlagOption("breakdown", &printBreakdown, "Print insertion/deletion/mismatch breakdown.", false);
	clp.RegisterFlagOption("errdist", &printDistBetweenErrorsHist, "Print a histogram of distance between errors.", false);
	clp.RegisterFlagOption("bymovie", &countBasesByMovie, "Count the number of bases aligned in each movie.", false);
  clp.RegisterStringOption("matchRun", &matchRunFileName, "Print lengths of runs of matches to file.", false);
  clp.RegisterFlagOption("discardFirstMatch", &discardFirstMatch, "Do not print the first run of matches.", false);
  clp.RegisterFlagOption("discardLastMatch", &discardLastMatch, "Do not print the last run of matches.", false);
  clp.RegisterStringOption("lengths", &lengthsFileName, "Print all subread lengths to a file.", false);
  clp.RegisterFlagOption("averageLength", &printAverageLength, "Print the average subread length.", false);
  clp.RegisterStringOption("printMovingAverageErrorRate", &movingAverageFileName, "Print moving average of accuracy.", false);
	clp.RegisterStringOption("binnedErrorDistribution", &binnedErrorDistributionFileName, 
													 "Print all binned error rates to a file", false);
  clp.RegisterStringOption("matchGap", &matchGapFileName,
                           "Print all matches and the gap between them to file.", false);
  clp.RegisterStringOption("matchesPerRead", &readMatchFileName,
                           "Print statistics for the number of matches found in a read", false);
  clp.RegisterIntOption("matchGapK", &matchGapK,
                        "(15) The minimum word size to match.",
                        CommandLineParser::PositiveInteger, false);
  clp.RegisterStringOption("errmat", &errorMatrixName, "Store error rates by read length", false);
	clp.ParseCommandLine(argc, argv);
  vector<vector<int>  > matchMatrix;
  vector<vector<int>  > errorMatrix;

  matchMatrix.resize(20);
  errorMatrix.resize(20);

	map<int,int> distHist;
  map<string,int> basesByMovie, readsByMovie;

	CmpFile cmpFile;
	
	/*
	 * These readers pull information from the same pls file.
	 */
	HDFCmpFile<CmpAlignment> cmpReader;

	if (cmpReader.Initialize(cmpFileName, H5F_ACC_RDONLY) == 0) {
		cout << "ERROR, could not open the cmp file." << endl;
		exit(1);
	}
  ofstream movingAverageFile;
  if (movingAverageFileName != "") {
    printMovingAverage= true;
    CrucialOpen(movingAverageFileName, movingAverageFile, std::ios::out);
  }
  ofstream matchCountFile;
  if (matchCountFileName != "") {
    CrucialOpen(matchCountFileName, matchCountFile, std::ios::out);
  }

	cmpReader.Read(cmpFile);
	int alignmentIndex;
	long nAlignments = 0;
	long totalAlignedBases = 0;
	long totalMatchedBases = 0;
	long totalAlignedLength = 0;
	long totalInsertion = 0, totalDeletion = 0, totalMismatch = 0;
	ofstream accuracyBinsFile, matchGapFile, readMatchFile, matchRunFile, lengthsFile, errMatFile;
  
	vector<vector<float> > accuracyDistributionBins;

  vector<int> mc, ic, dc, mmc;

	if (binnedErrorDistributionFileName != "") {
		CrucialOpen(binnedErrorDistributionFileName, accuracyBinsFile, std::ios::out);
		accuracyDistributionBins.resize(nBins);
	}
	
  if (lengthsFileName != "") {
    CrucialOpen(lengthsFileName, lengthsFile, std::ios::out);
  }

  if (errorMatrixName != "") {
    CrucialOpen(errorMatrixName, errMatFile, std::ios::out);
  }

	vector<float> accuracyBins;
	if (nBins > 0) {
		accuracyBins.resize(nBins);
	}
  
  if (matchRunFileName != "") {
    CrucialOpen(matchRunFileName, matchRunFile, std::ios::out);
    matchRunFile << "run_length gap_length" << endl;
  }

  if (matchGapFileName != "") {
    CrucialOpen(matchGapFileName, matchGapFile, std::ios::out);
    matchGapFile  << "match_length gap_length" << endl;
  }
	
  if (readMatchFileName != "") {
    CrucialOpen(readMatchFileName, readMatchFile, std::ios::out);
    readMatchFile << "read_length nmatches" << endl;
  }
  
	for (alignmentIndex = 0; alignmentIndex < cmpFile.alnInfo.alignments.size(); alignmentIndex++) {

		vector<unsigned char> byteAlignment;
		UInt offsetBegin, offsetEnd;
		
		offsetBegin = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetBegin();
		offsetEnd   = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetEnd();
		int subreadLength = (cmpFile.alnInfo.alignments[alignmentIndex].GetQueryEnd() - 
                         cmpFile.alnInfo.alignments[alignmentIndex].GetQueryStart());

    if (lengthsFileName != "") {
      lengthsFile << subreadLength << endl;
    }

		int alignedSequenceLength = offsetEnd - offsetBegin;
		string alignedSequence;
		if (alignedSequenceLength >= 0) {
			alignedSequence.resize(alignedSequenceLength);
			byteAlignment.resize(alignedSequenceLength);
		}
	
		// Read the alignment string.  All alignments 
		//
		// Alignments are grouped by ref group id then movie id.
		//
		int refGroupId = cmpFile.alnInfo.alignments[alignmentIndex].GetRefGroupId();
        int alnGroupId = cmpFile.alnInfo.alignments[alignmentIndex].GetAlnGroupId();
			
		//
		// Now locate where this movie is stored.
		//
		if (cmpReader.refGroupIdToArrayIndex.find(refGroupId) == cmpReader.refGroupIdToArrayIndex.end()) {
			cout << "ERROR!  An alignment " << alignmentIndex << " is specified with reference group " << endl
					 << refGroupId << " that is not found as an alignment group." << endl;
			exit(1);
		}
		int refGroupIndex    = cmpReader.refGroupIdToArrayIndex[refGroupId];
        string readGroupName = cmpReader.alnGroupIdToReadGroupName[alnGroupId];
        int readGroupIndex   = cmpReader.refAlignGroups[refGroupIndex]->experimentNameToIndex[readGroupName];

		cmpReader.refAlignGroups[refGroupIndex]->readGroups[readGroupIndex]->alignmentArray.Read(offsetBegin, 
																																														 offsetEnd, 
																																														 &byteAlignment[0]);
    /*
    if (matchRunFileName != "") {
      PrintAlignment(byteAlignment, matchRunFile);
    }
    */
    int n1 = 0, ngt1 = 0;
    if (matchGapFileName != "" or readMatchFileName != "" or matchRunFileName != "") {
      //
      // Print the matches and the gaps in this alignment to 
      // the file.
      //
      int i = 0;
      int alignEnd = byteAlignment.size() - matchGapK;
      bool matchFound = false;
      int prevMatchEnd = 0;
      int prevMatchLength = 0;
      int nMatches = 0;
      int prevMatchAnyEnd = 0;
      while (i < alignEnd) {
        // Find the first matching character.
        int nonMatchLength = 0;
        while (i < alignEnd and 
               (QueryChar[byteAlignment[i]] == ' ' or
                RefChar[byteAlignment[i]] == ' ' or 
                QueryChar[byteAlignment[i]] != RefChar[byteAlignment[i]])) {
          i++;
          nonMatchLength++;
        }
        if (i >= alignEnd) {
          break;
        }
        // find the end of this match
        int matchStart = i;
        while (i < alignEnd and QueryChar[byteAlignment[i]] != ' ' 
               and RefChar[byteAlignment[i]] != ' ' 
               and QueryChar[byteAlignment[i]] == RefChar[byteAlignment[i]]) {
          i++;
        }
        int matchEnd = i;
               
        if (matchRunFileName != "") {
          bool printAlignment = true;
          if (discardFirstMatch and matchStart == 0) {
            printAlignment = false;
          }
          if (discardLastMatch and matchEnd == alignEnd) {
            printAlignment = false;
          }
          if (printAlignment) {
            matchRunFile << matchEnd - matchStart << " " << nonMatchLength << endl;
          }
          if  (matchEnd - matchStart == 1) {
            n1++;
          }
          else { ngt1++;}
        }
        prevMatchAnyEnd = i;
        // If this match counts as an anchor, process it.
        if (i - matchStart >= matchGapK) {
          // Processing starts by looking to see if a previous anchor was found
          // and if so, printing that and the gap to the current anchor.
          //
          if (matchFound == true) {
            if (matchGapFileName != "") {
              matchGapFile << prevMatchLength << " " << matchStart - prevMatchEnd << endl;
            }
            ++nMatches;
          }
          // Processing ends by storing the length of the match, and where
          // it ended so that the next iteration can use it.
          matchFound = true;
          prevMatchLength = i - matchStart;
          prevMatchEnd    = i;
        }
      }
      if (readMatchFileName != "") {
        readMatchFile << byteAlignment.size() << " " << nMatches << endl;
      }
    }
    if (matchRunFileName != "") {
      //      matchRunFile << "n1: " << n1 << " ngt1 " << ngt1 << endl;
    }

    if (countBasesByMovie) {
      string movieName = cmpFile.movieInfo.name[readGroupIndex];
      if (basesByMovie.find(movieName) ==
          basesByMovie.end()) {
        readsByMovie[movieName] = 0;
        basesByMovie[movieName] = 0;
      }
      basesByMovie[movieName]+= offsetEnd - offsetBegin;
      readsByMovie[movieName]++;
    }
      
		if (printDistBetweenErrorsHist) {
			AddDistancesBetweenErrors(byteAlignment, distHist);
		}

		float readPctIdentity;
		readPctIdentity = ComputePacBioAccuracy(byteAlignment);
		
		if (readPctIdentity < identityCutoff) {
			continue;
		}

		if (byteAlignment.size() < minAlignLength) {
			continue;
		}

		++nAlignments;
    //
    // several stats use total aligned length
    int qStart = cmpFile.alnInfo.alignments[alignmentIndex].GetQueryStart();
    int qEnd   = cmpFile.alnInfo.alignments[alignmentIndex].GetQueryEnd();
    totalAlignedBases += qEnd - qStart;

    if (errorMatrixName != "") {
      StoreErrorRateMatrix(byteAlignment, matchMatrix, errorMatrix);
    }

		if (printBinnedErrorRate) {
			StoreBinnedErrorRate(byteAlignment, accuracyBins);
		}

		if (binnedErrorDistributionFileName != "") {
			AppendBinnedErrorRate(byteAlignment, accuracyDistributionBins);
		}
		
		int nMatch, nMismatch, nIns, nDel;
		CountStats(byteAlignment, nMatch, nMismatch, nIns, nDel);
		float total = nMatch + nMismatch + nIns + nDel;

		totalMatchedBases    += nMatch;
		totalMismatch        += nMismatch;
		totalInsertion       += nIns;
		totalDeletion        += nDel;
		//		totalMatchedBases    += CountNMatches(byteAlignment);
		totalAlignedLength   += byteAlignment.size();
		totalPercentIdentity += ComputePercentIdentity(byteAlignment);

    if (printMovingAverage) {
      vector<float> movingAverage;
      StoreMovingAverage(byteAlignment, movingAverage);
      int i;
      for (i = 0; i < movingAverage.size(); i++) {
        movingAverageFile << i << " " << movingAverage[i] << endl;
      }
    }

    if (matchCountFileName != "") {
      StoreMatchCounts(byteAlignment, mc, ic, dc, mmc);
    }
	}

  if (countBasesByMovie) {
    map<string,int>::iterator mapIt;
    for (mapIt = basesByMovie.begin(); mapIt != basesByMovie.end(); ++mapIt) {
      cout << mapIt->first << " " << readsByMovie[mapIt->first] << " " << mapIt->second << endl;;
    }
  }

  if (matchCountFileName != "") {
    int i;
    for (i = 0; i < mc.size(); i++) {
      matchCountFile << mc[i] << " " << ic[i] << " " << dc[i] << " " << mmc[i] << endl;
    }
    matchCountFile.close();
  }

	if (printNReads) {
		cout << "nAlignments\t"<< nAlignments << endl;
	}

	if (printTotalAlignedBases) {
		cout << "totalAlignedBases\t" << totalAlignedBases << endl;
	}

	if (printBreakdown) {
		float totalTemplate = totalMatchedBases + totalMismatch + totalDeletion + totalInsertion;
		cout << "M_MM_I_D " << totalMatchedBases << " " << totalMismatch << " " << totalInsertion << " " << totalDeletion << " "
				 << totalMatchedBases / totalTemplate << " " 
				 << totalMismatch     / totalTemplate << " " 
				 << totalInsertion    / totalTemplate << " " 
				 << totalDeletion     / totalTemplate << endl;
	}

	if (printNMatches) {
		cout << "totalMatches\t" << totalMatchedBases << endl;
	}

	if (printGlobalAccuracy) {
		cout << "globalAccuracy\t"<< ((float)totalMatchedBases) / totalAlignedLength << endl;
	}

	if (printAverageAccuracy) {
		cout << "averageAccuracy\t" << totalPercentIdentity / nAlignments << endl;
	}

  if (printAverageLength) {
    cout << "averageAlignmentLength\t" << totalAlignedBases / (float)nAlignments << endl;
  }

  
	if (printBinnedErrorRate) {
		int i;
		for (i = 0; i < nBins; i++) {
			accuracyBins[i] /= nAlignments;
			cout << accuracyBins[i] << " ";
		}
		cout << endl;
	}

	if (binnedErrorDistributionFileName != "") {
		int i, b;
		for (i = 0; i < accuracyDistributionBins[0].size(); i++) {
			for (b = 0; b < nBins; b++) {
				accuracyBinsFile << accuracyDistributionBins[b][i] << " ";
			}
			accuracyBinsFile << endl;
		}
		accuracyBinsFile.close();
	}

  if (lengthsFileName != "") {
    lengthsFile.close();
  }
    
	if (printDistBetweenErrorsHist) {
		map<int,int>::iterator histIt;
		for (histIt = distHist.begin(); histIt != distHist.end(); ++histIt) {
			cout << "hist " << histIt->first << " " << histIt->second << endl;
		}
	}
	if (printMovingAverage){ 
    movingAverageFile.close();
  }
  if (errorMatrixName != "") {
    int i, j;
    for (i = 0; i< matchMatrix.size(); i++) {
      for (j = 0; j < matchMatrix[i].size(); j++) {
        if ((matchMatrix[i][j] + errorMatrix[i][j]) > 0) {
          errMatFile << ((float)matchMatrix[i][j]) / ((matchMatrix[i][j] + errorMatrix[i][j])) << " ";
        }
        else {
          errMatFile << " 0 ";
        }
      }
      errMatFile << endl;
    }
    errMatFile.close();
  }

	return 0;
}
