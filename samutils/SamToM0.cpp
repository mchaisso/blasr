#include "algorithms/alignment/readers/sam/SAMReader.h"
#include "datastructures/alignmentset/SAMAlignment.h"
#include "datastructures/alignmentset/SAMToAlignmentCandidateAdapter.h"
#include "FASTAReader.h"
//#include "datastructures/sequence/FAITable.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "algorithms/alignment/printers/StickAlignmentPrinter.h"
#include "CommandLineParser.h"
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

class Interval {
public:
  int index;
  unsigned int start;
  unsigned int end;
  Interval(int i, unsigned int st, unsigned int en): index(i), start(st), end(en) {}
  Interval() {}
  bool operator<(const Interval &rhs) const {
    return start < rhs.start;
  }
};


typedef set<Interval> IntervalSet;
typedef map<string, IntervalSet> RegionMap;


void ReadROIFile(istream &roiIn, RegionMap &regions) {
  int lineNumber = 0;

  while(roiIn){
    string chr;
    unsigned int start, end;
    if ( (roiIn >> chr >> start >> end) ) {
      regions[chr].insert(Interval(lineNumber, start,end));
      ++lineNumber;
    }
  }
}


bool FindRegion(RegionMap &regions, string chr, unsigned int start, unsigned int end,
                Interval &foundRegion) {
  if (regions.find(chr) == regions.end()) {
    return false;
  }
  else {
    //
    // Scan this region for a match.  It is expected that start <
    // intv.start and end > intv.end.
    //
    IntervalSet::iterator it;
    Interval tmpIntv;
    tmpIntv.start = start;
    it = regions[chr].lower_bound(tmpIntv);
    while (it != regions[chr].end() and
           (*it).start < end) {
      if ((*it).start >= start and (*it).start < end) {
        foundRegion = (*it);
        return true;
      }
      ++it;
    }
    return false;
  }
}

int main(int argc, char* argv[]) {
  string samFileName, genomeFileName;
  string outFileName = "/dev/stdout";
  string alignmentsFileName = "";
  int minGapLength = 100;
  string gapFastaFileName = "";
  string regionOfInterestFileName = "";
  CommandLineParser clp;
  int minMapQV = 0;
	vector<string> samFileNames;
	int excludeFlag = 0;
	int expand = 0;
	int maxMerge = 0;
	int format = 0;
	int minLength = 0;
	bool forward = false;
  clp.RegisterStringOption("genome", &genomeFileName, "Genome.", true);
  clp.RegisterStringListOption("sam", &samFileNames, "Alignments.", true);
  clp.RegisterPreviousFlagsAsHidden();
  clp.RegisterStringOption("out", &outFileName, "Output file. Default to stdout", false);
  clp.RegisterIntOption("format", &format, "Format (supports 0 only, formats 1,4,5 to be added)", CommandLineParser::NonNegativeInteger, false);
  clp.RegisterFlagOption("forward", &forward, "Print ref in forward direction", false);	
  clp.RegisterIntOption("minLength", &minLength, "Skip alignments under this length.", CommandLineParser::NonNegativeInteger, false);
  clp.ParseCommandLine(argc, argv);

  FASTAReader fastaReader;
	//	FAITable reference;
  fastaReader.Initialize(genomeFileName);

  ostream *outPtr;
  ofstream outFile;
  ofstream alignmentsOut, gapFastaOut;
	CrucialOpen(outFileName, outFile, std::ios::out);

  if (gapFastaFileName != "") {
    CrucialOpen(gapFastaFileName, gapFastaOut, std::ios::out);
  }
	int samFileIndex;
  ifstream roiIn;
  RegionMap regions;
  vector<FASTASequence> references;


  fastaReader.ReadAllSequences(references);
	
  int i;
  map<string,int> refToIndex;
  for (i = 0; i < references.size(); i++) {
		string name = references[i].GetName();
		references[i].CopyTitle(name.c_str());
    refToIndex[references[i].title] = i;
  }
	for (samFileIndex=  0; samFileIndex < samFileNames.size(); samFileIndex++) {
		string samFileName = samFileNames[samFileIndex];
	
		SAMReader<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> samReader;
		samReader.Initialize(samFileName);
		cerr << samFileName << endl;
		AlignmentSet<SAMFullReferenceSequence, SAMReadGroup, SAMPosAlignment> alignmentSet;
		samReader.ReadHeader(alignmentSet);

		SAMAlignment samAlignment;  

		int index = 0;

		while (samReader.GetNextAlignment(samAlignment)) {
			if (samAlignment.rName == "*") {
				continue;
			}
			if (samAlignment.flag & excludeFlag != 0) {
				continue;
			}
			if (forward) {
				samAlignment.flag = 0;
			}
			vector<AlignmentCandidate<> > convertedAlignments;
			SAMAlignmentsToCandidates(samAlignment, 
																references, refToIndex,
																convertedAlignments, false, false);
			
			int a;
			for (a = 0; a < convertedAlignments.size(); a++) {
				if (format == 0) {
					if (convertedAlignments[a].qAlignedSeq.length >= minLength) {
						StickPrintAlignment(convertedAlignments[a], convertedAlignments[a].qAlignedSeq, convertedAlignments[a].tAlignedSeq, outFile, convertedAlignments[a].qAlignedSeqPos,convertedAlignments[a].tAlignedSeqPos  );
					}
				}
				else {
					cout << "Other formats than 0 coming soon!" << endl;
				}
			}
    }
  }
}
