#include "data/hdf/HDFPlsReader.h"
#include "data/hdf/HDFRegionTableReader.h"
#include "datastructures/reads/RegionTable.h"
#include "datastructures/reads/ReadInterval.h"
#include "files/ReaderAgglomerate.h"
#include "utils/FileOfFileNames.h"
#include "utils/RegionUtils.h"
#include "SMRTSequence.h"
#include "utils.h"
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include "CommandLineParser.h"

using namespace std;

/* 
   Some notes about this program:
     It is in mid-refactor state.  


*/


int main(int argc, char* argv[]) {

	string plsFileName, fastaOutName;
	vector<string> plsFileNames;
	bool trimByRegion, maskByRegion;
	trimByRegion = false;
	maskByRegion = false;
	int argi = 3;
	RegionTable regionTable;
	string regionsFOFNName = "";
	vector<string> regionFileNames;
	bool splitSubreads = true;
	int minSubreadLength = 0;
	bool addSimulatedData = false;
	bool printSimulatedCoordinate = false;
	bool printSimulatedSequenceIndex = false;
  bool printFastq = false;
  bool printCcs   = false;
  int  lineLength = 50;
  int minReadScore = 0;
  bool encodeCCSPassesInTitle = false;
  vector<int> holeNumbers;
  CommandLineParser clp;
  bool printOnlyBest = false;
  clp.SetProgramName("pls2fasta");
  clp.RegisterStringOption("file.pls.h5", &plsFileName, "Input pls/bas.h5 file.", true);
  clp.RegisterStringOption("out.fasta", &fastaOutName, "Output fasta/fastq file.", true);
  clp.RegisterPreviousFlagsAsHidden();
  clp.RegisterFlagOption("trimByRegion", &trimByRegion, "Trim away low quality regions.");
  clp.RegisterFlagOption("maskByRegion", &maskByRegion, "Mask low quality regions with 'N'.");
  clp.RegisterStringOption("regionTable", &regionsFOFNName, "Optional HDF file with a /PulseData/Regions dataset.");
  clp.RegisterIntOption("minSubreadLength", &minSubreadLength, "Do not write subreads less than the specified length.", CommandLineParser::PositiveInteger);
  clp.RegisterFlagOption("noSplitSubreads", &splitSubreads, "Do not split reads on adapter sequences.");
  clp.RegisterIntListOption("holeNumber", &holeNumbers, "Only print this hole number (or list of numbers).");
  clp.RegisterFlagOption("fastq", &printFastq, "Print in FASTQ format with quality.");
  clp.RegisterFlagOption("ccs", &printCcs, "Print de novo CCS sequences");
  clp.RegisterFlagOption("passesInTitle", &encodeCCSPassesInTitle, "Append /N_passes/ to ccs sequence title, where"
                         " N is the number of passes.");
  clp.RegisterIntOption("lineLength", &lineLength, "Specify fasta/fastq line length", CommandLineParser::PositiveInteger);
  clp.RegisterIntOption("minReadScore", &minReadScore, "Minimum read score to print a read.  The score is "
                        "a number between 0 and 1000 and represents the expected accuracy percentage * 10. "
                        "A typical value would be between 750 and 800.  This does not apply to ccs reads.", CommandLineParser::NonNegativeInteger);
  clp.RegisterFlagOption("best", &printOnlyBest, "If a CCS sequence exists, print this.  Otherwise, print the longest"
                         "subread.  This does not support fastq.");
  clp.SetProgramSummary("Converts bas.h5 files to fasta or fastq files. Although fasta files are provided"
                        " with every run, they are not trimmed nor split into subreads. This program takes "
                        "additional annotation information, such as the subread coordinates and high quality regions "
                        "and uses them to create fasta sequences that are substrings of all bases called. Most of the time "
                        "you will want to trim low quality reads, so you should specify -trimByRegion.");
                        
  clp.ParseCommandLine(argc, argv);

	if (trimByRegion and maskByRegion) {
		cout << "ERROR! You cannot both trim and mask regions. Use one or the other." << endl;
		exit(1);
	}
		 
  if (printFastq) {
    // Setting lineLength to 0 flags to print on one line.
    lineLength = 0;
  }

	if (FileOfFileNames::IsFOFN(plsFileName)) {
		FileOfFileNames::FOFNToList(plsFileName, plsFileNames);
	}
	else {
		plsFileNames.push_back(plsFileName);
	}
	if (regionsFOFNName == "") {
		regionFileNames = plsFileNames;
	}
	else {
		if (FileOfFileNames::IsFOFN(regionsFOFNName)) {
			FileOfFileNames::FOFNToList(regionsFOFNName, regionFileNames);
		}
		else {
			regionFileNames.push_back(regionsFOFNName);
		}
	}



	ofstream fastaOut;
	CrucialOpen(fastaOutName, fastaOut);
	int plsFileIndex;
	HDFRegionTableReader hdfRegionReader;
  sort(holeNumbers.begin(), holeNumbers.end());
	for (plsFileIndex = 0; plsFileIndex < plsFileNames.size(); plsFileIndex++) {

		if (trimByRegion or maskByRegion or splitSubreads) {
			hdfRegionReader.Initialize(regionFileNames[plsFileIndex]);
			hdfRegionReader.ReadTable(regionTable);
			regionTable.SortTableByHoleNumber();
		}
    HDFCCSReader<CCSSequence> ccsReader;		
		HDFBasReader smrtReader;
		if (printCcs) {
			ccsReader.Initialize(plsFileNames[plsFileIndex]);
		}
		else {
			smrtReader.Initialize(plsFileNames[plsFileIndex]);
		}

		vector<ReadInterval> subreadIntervals;;
    SMRTSequence seq;
    CCSSequence  ccsSeq;

		while (true) {
			if (printCcs == true or printOnlyBest == true) {
				if (ccsReader.GetNext(ccsSeq) == false) {
					break;
				}
				else {
					seq = ccsSeq.unrolledRead;
				}
			}
			else {
				if (smrtReader.GetNext(seq) == false) {
					break;
				}
			}

      if (holeNumbers.size() != 0 and 
          binary_search(holeNumbers.begin(), holeNumbers.end(), seq.zmwData.holeNumber) == false) {
        continue;
      }

      if (encodeCCSPassesInTitle) {
				assert(printCcs);
        string title = ccsSeq.title;
        stringstream titleStrm;
        titleStrm << title << "/"<<ccsSeq.numPasses << "_passes/";
        ccsSeq.CopyTitle(titleStrm.str());
      }

		  if (printCcs == true or (printOnlyBest and ccsSeq.length > 0)) {
        //
        // The reason for checking to see if the length is greater
        // than 0 again is in case -ccs flag is specified, but the ccs
        // sequence is empty, nothing should be printed.  When
        // printing only the best this block should only be entered if
        // a ccs sequence exists, because the rest of the loop is
        // skipped after this.

        if (ccsSeq.length > 0) {
          if (printFastq == false) {
            ccsSeq.PrintSeq(fastaOut);
          }
          else {
            ccsSeq.PrintFastq(fastaOut, lineLength);
          }
        }
				ccsSeq.Free();
				seq.Free();
        continue;
      }	


      if (seq.length == 0) {
        continue;
      }


      //
      // Determine the high quality boundaries of the read.  This is
      // the full read is no hq regions exist, or it is stated to
      // ignore regions.
      //
      DNALength hqReadStart, hqReadEnd;
      int hqRegionScore;
      if (GetReadTrimCoordinates(seq, seq.zmwData, regionTable, hqReadStart, hqReadEnd, hqRegionScore) == false or 
          (trimByRegion == false and maskByRegion == false)) {
        hqReadStart = 0;
        hqReadEnd   = seq.length;
      }
      
      //
      // Mask off the low quality portions of the reads.
      //
			if (maskByRegion) {
        if (hqReadStart > 0) {
          fill(&seq.seq[0], &seq.seq[hqReadStart], 'N');
        }
        if (hqReadEnd != seq.length) {
          fill(&seq.seq[hqReadEnd], &seq.seq[seq.length], 'N');
        }
			}


      //
      // Now possibly print the full read with masking.  This could be handled by making a 
      // 
			if (splitSubreads == false) {
        ReadInterval wholeRead(0, seq.length);
        // The set of subread intervals is just the entire read.
        subreadIntervals.clear();
        subreadIntervals.push_back(wholeRead);
			}
			else {
				//
				// Print subread coordinates no matter whether or not reads have subreads.
				//
				subreadIntervals.clear(); // clear old, new intervals are appended.
				CollectSubreadIntervals(seq, &regionTable, subreadIntervals);
      }
      //
      // Output all subreads as separate sequences.
      //
      int intvIndex;
      SMRTSequence bestSubreadSequence;
      int bestSubreadScore = -1;
      int bestSubreadIndex = 0;
      int bestSubreadStart = 0, bestSubreadEnd = 0;
      SMRTSequence bestSubread;
      for (intvIndex = 0; intvIndex < subreadIntervals.size(); intvIndex++) {
        SMRTSequence subreadSequence, subreadSequenceRC;
					
        subreadSequence.subreadStart = subreadIntervals[intvIndex].start;
        subreadSequence.subreadEnd   = subreadIntervals[intvIndex].end;
          
        // 
        // When trimming by region, only output the parts of the
        // subread that overlap the hq region.
        //
        if (trimByRegion == true) {
          subreadSequence.subreadStart = max((DNALength) subreadIntervals[intvIndex].start, hqReadStart);
          subreadSequence.subreadEnd   = min((DNALength) subreadIntervals[intvIndex].end, hqReadEnd);
        }

        if (subreadSequence.subreadStart >= subreadSequence.subreadEnd or 
            subreadSequence.subreadEnd - subreadSequence.subreadStart <= minSubreadLength) {
          //
          // There is no high qualty portion of this subread. Skip it.
          //
          continue;
        }

        if (hqRegionScore < minReadScore) {
          continue;
        }

        //
        // Print the subread, adding the coordinates as part of the title.
        //
        subreadSequence.ReferenceSubstring(seq, subreadSequence.subreadStart, 
                                           subreadSequence.subreadEnd - subreadSequence.subreadStart);
        stringstream titleStream;
        titleStream << seq.title;
        if (splitSubreads) {
          //
          // Add the subread coordinates if splitting on subread.
          //
          titleStream << "/" 
                      << subreadSequence.subreadStart
                      << "_" << subreadSequence.subreadEnd;
        }
          
        subreadSequence.CopyTitle(titleStream.str());

        //
        // Eventually replace with WriterAgglomerate.
        //
        if (printOnlyBest == false) {
          if (subreadSequence.length > 0) {
            if (printFastq == false) {
              ((FASTASequence*)&subreadSequence)->PrintSeq(fastaOut);
            }
            else {
              subreadSequence.PrintFastq(fastaOut, lineLength);
            }
          }
          delete[] subreadSequence.title;
        }
        else {
          int subreadWeightedScore = subreadSequence.length * hqRegionScore;
          if (subreadWeightedScore > bestSubreadScore) {
            bestSubreadIndex = intvIndex;
            bestSubread = subreadSequence;
            bestSubreadScore = subreadWeightedScore;
          }
        }
      }

      if (printOnlyBest) {
        if (bestSubreadScore >= 0) {
          if (printFastq == false) {
            bestSubread.PrintSeq(fastaOut);
          }
          else {
            bestSubread.PrintFastq(fastaOut, bestSubread.length);
          }
          bestSubread.Free();
        }
      }
      ccsSeq.Free();
      seq.Free();
    }
		if (printCcs or printOnlyBest) {
			ccsReader.Close();
		}
		else {
			smrtReader.Close();
		}
    hdfRegionReader.Close();
  }
}
