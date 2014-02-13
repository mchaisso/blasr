#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFRegionTableReader.h"
#include "SMRTSequence.h"
#include "algorithms/alignment/AffineGuidedAlign.h"
#include "algorithms/alignment/SDPAlign.h"
#include "algorithms/alignment/DistanceMatrixScoreFunction.h"
#include "algorithms/alignment/IDSScoreFunction.h"
#include "algorithms/alignment/ScoreMatrices.h"
#include "utils/FileOfFileNames.h"
#include "utils/RegionUtils.h"
#include "CommandLineParser.h"
#include "datastructures/reads/ReadInterval.h"
#include "datastructures/reads/RegionTable.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "data/hdf/HDFRegionTableWriter.h"
#include "algorithms/alignment/printers/StickAlignmentPrinter.h"
#include "statistics/statutils.h"
#include <sys/stat.h>

class Metrics {
public:
  int numFiltered;
  vector<int> criticalAccuracy;
  vector<int> criticalLength;
  vector<int> gapLength;
  vector<int> nSubreads;
  string movieName;
  Metrics() {
    numFiltered = 0;
  }

  void Print(ostream &out) {
    out << movieName << endl;
    out.precision(2);
    float meanCA, varCA, meanCL, varCL;
    float meanGL, varGL;
    MeanVar(criticalAccuracy, meanCA, varCA);
    MeanVar(criticalLength, meanCL, varCL);
    MeanVar(gapLength, meanGL, varGL);
    
    out << "critical_accuracy\t" << criticalAccuracy.size() << "\t"
        << ios::fixed  << meanCA << "\t"<< varCA << endl;
    out << "critical_length\t" << criticalLength.size() << "\t"
        << ios::fixed  << meanCL << "\t"<< varCL << endl;
    out << "gap_length\t" << gapLength.size() << "\t"
        << ios::fixed  << meanGL << "\t"<< varGL << endl;
    
    int numMultiPass = 0;
    int i;
    for (i = 0; i < nSubreads.size(); i++) {
      if (nSubreads[i] > 1) {
        numMultiPass++;
      }
    }
    out << "multi_pass\t"<<numMultiPass << "\t" << nSubreads.size() << endl;
  }
};

class MappingBuffers {
public:
  vector<int> hpInsScoreMat, insScoreMat;
  vector<int> kbandScoreMat;
  vector<Arrow> hpInsPathMat, insPathMat;
  vector<Arrow> kbandPathMat;
  vector<int>   scoreMat;
  vector<Arrow> pathMat;
  vector<int>  affineScoreMat;
  vector<Arrow> affinePathMat;
  TupleList<PositionDNATuple> sdpCachedTargetTupleList;
  TupleList<PositionDNATuple> sdpCachedTargetPrefixTupleList;
  TupleList<PositionDNATuple> sdpCachedTargetSuffixTupleList;
  std::vector<int> sdpCachedMaxFragmentChain;
  vector<double> probMat;
  vector<double> optPathProbMat;
  vector<float>  lnSubPValueMat;
  vector<float>  lnInsPValueMat;
  vector<float>  lnDelPValueMat;
  vector<float>  lnMatchPValueMat;
};


int AlignmentHasIndel(AlignmentCandidate<> aln, int indelLength) {
  int g, gi;
  if (aln.gaps.size() == 0) {
    return false;
  }
  for (g = 1; g < aln.gaps.size() - 1; g++) {
    for (gi = 0; gi < aln.gaps[g].size(); gi++) {
      if (aln.gaps[g][gi].length > indelLength) {
        return aln.gaps[g][gi].length;
      }
    }
  }
  return false;
}
  

int main(int argc, char* argv[]) {
  string inFileName;
  string regionOutFileName;
  string regionInFileName = "";
  int sdpTupleSize = 12;
  CommandLineParser clp;
  int insertion = 20;
  int deletion = 20;
  int minSubreadLength = 100;
  string criticalAccuracyFileName = "";
  string alignmentsFileName = "";
  int    criticalAccuracyMinLength = 0;
  int    criticalAccuracyMinAccuracy = 0;
  float  criticalAccuracyLengthFraction = 0.0;
  bool   filterByCriticalAccuracy = false;
  bool   computeCriticalAccuracy = false;
  string filteredReadNameFileName = "";
  bool   printOnlyFiltered = false;
  int    indelLength = 0;
  bool   filterByIndel = false;
  bool   printRegions = false;
  string regionOutDirName = "";
  string regionFofnFileName = "";
  string metricsFileName = "";
  clp.RegisterStringOption("inFile", &inFileName, "Input bas.h5 or fofn.", true);
  clp.RegisterPreviousFlagsAsHidden();
  clp.RegisterStringOption("regionDir", &regionOutDirName, "Write new region tables into this directory.", false);
  clp.RegisterStringOption("regionFofn", &regionFofnFileName, "Write a region table FOFN .", false);
  clp.RegisterStringOption("regionTable", &regionInFileName, "Alternative region table file to read.");
  clp.RegisterIntOption("sdpTupleSize", &sdpTupleSize, "Word size to generate alignment guide.", CommandLineParser::PositiveInteger);
  clp.RegisterIntOption("minSubreadLength", &minSubreadLength, "Do not try and filter subreads less than this length.", CommandLineParser::PositiveInteger);
  clp.RegisterStringOption("printAlignments", &alignmentsFileName, "Print alignments to this file.", false);
  clp.RegisterFlagOption("printOnlyFiltered", &printOnlyFiltered, "When printing alignments, only print those of sequences that have large abberations.", false);
  clp.RegisterStringOption("ca", &criticalAccuracyFileName, "Measure the critical accuracy table and write it to this file.", false);
  clp.RegisterIntOption("filterByIndel", &indelLength, "Flag reads that have an indel greater than this length.", CommandLineParser::PositiveInteger, false);
  clp.RegisterIntOption("filterByCriticalAccuracy", &criticalAccuracyMinAccuracy, 
                        "Set the accuracy for filtering by critical accuracy.  This should be used in conjunction with the "
                        "-criticalAccuracyLength flag.  Let a be the critical accuracy-accuracy (in [0,100]), and "
                        "l be the critical accuracy-length, if a read has a stretch longer than l of accuracy a (or less) "
                        "it is flagged as filtered.", CommandLineParser::PositiveInteger, false);
  clp.RegisterIntOption("filterByCriticalAccuracyLength", &criticalAccuracyMinLength,
                        "Set the length of a critical accuracy.", CommandLineParser::PositiveInteger,false);
  clp.RegisterFloatOption("filterByCriticalAccuracyFraction", &criticalAccuracyLengthFraction,
                          "Compute the critical accuracy length as a fraction of the subread length, in [0,1]).",
                          CommandLineParser::PositiveFloat, false);
  clp.RegisterStringOption("filteredReadNames", &filteredReadNameFileName, 
                           "Print the names of reads that are filtered to this file.", false);
  clp.RegisterStringOption("metricsFile", &metricsFileName, 
                           "Print some statistics about the filtering.", false);

  clp.ParseCommandLine(argc, argv);
  if (criticalAccuracyFileName != "") {
    filterByCriticalAccuracy = true;
  }
  if (criticalAccuracyMinAccuracy > 0) {
    filterByCriticalAccuracy = true;
    computeCriticalAccuracy  = true;
    if (criticalAccuracyMinAccuracy < 0 or criticalAccuracyMinAccuracy > 100) {
      cout << "ERROR, the critical accuracy must be between 0 and 100 (and an integer)." << endl;
      exit(1);
    }
    if (criticalAccuracyLengthFraction > 1 or criticalAccuracyLengthFraction < 0) {
      cout << "ERROR, the critical accuracy fraction must be between 0 and 1. (0 is off)."; 
      exit(1);
    }
  }

  if (regionOutDirName == "" and regionFofnFileName != "") {
    cout << "ERROR. I a region fofn is specified, a directory to write the output must be specified, and it cannot be \".\"" << endl;
    exit(1);
  }

  if (regionOutDirName != "") {
    if (regionOutDirName == ".") {
      cout << "ERROR.  Do not use '.' for an output directory.  Think of a a better place" << endl
           << "to put region files." << endl;
      exit(1);
    }

    printRegions = true;
    int retval;
    retval = mkdir(regionOutDirName.c_str(), S_IRWXU);
    if (retval != 0) {
      cout << "ERROR, could not create the directory " << regionOutDirName << endl;
      exit(1);
    }
  }
  vector<string> readFileNames, regionFileNames;

  MappingBuffers mappingBuffers;

  if (regionInFileName == "") {
    regionInFileName = inFileName;
  }

  ofstream caFile;
  if (criticalAccuracyFileName != "") {
    CrucialOpen(criticalAccuracyFileName, caFile, std::ios::out);
    computeCriticalAccuracy = true;
  }

  ofstream alnFile;
  if (alignmentsFileName != "") {
    CrucialOpen(alignmentsFileName, alnFile, std::ios::out);
  }

  ofstream filteredReadNameFile;
  if (filteredReadNameFileName != "") {
    CrucialOpen(filteredReadNameFileName, filteredReadNameFile, std::ios::out);
  }

  ofstream regionFofnFile;
  if (regionFofnFileName != "") {
    CrucialOpen(regionFofnFileName, regionFofnFile, std::ios::out);
  }
  ofstream metricsFile;
  if (metricsFileName != "") {
    CrucialOpen(metricsFileName, metricsFile, std::ios::out);
  }
  
  DistanceMatrixScoreFunction<SMRTSequence, SMRTSequence> distScoreFn;
  distScoreFn.del = deletion;
  distScoreFn.ins = insertion;
  distScoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);


  FileOfFileNames::StoreFileOrFileList(inFileName, readFileNames);
  if (regionInFileName != "") {
    FileOfFileNames::StoreFileOrFileList(regionInFileName, regionFileNames);
  }
  else {
    regionFileNames = readFileNames;
  }
  
  int f;
  for (f = 0; f < readFileNames.size(); f++) {

    HDFBasReader reader;
    reader.InitializeDefaultIncludedFields();
    reader.Initialize(readFileNames[f]);
    Metrics metrics;

    HDFRegionTableReader regionTableReader;

    if (regionTableReader.Initialize( regionFileNames[f]) == false) {
      cout <<"ERROR, could not read region table " << inFileName << endl;
      exit(1);
    }
    
    HDFRegionTableWriter regionTableWriter;
    if (printRegions) {
      string movieName = reader.GetMovieName();
      regionOutFileName = regionOutDirName + "/" + movieName + ".rgn.h5";
      regionTableWriter.Create(regionOutFileName);
      if (regionFofnFileName != "") {
        regionFofnFile << regionOutFileName << endl;
      }
    }

    
    SMRTSequence smrtRead;
    RegionTable regionTable;
    regionTableReader.ReadTable(regionTable);
    IDSScoreFunction<DNASequence, FASTQSequence> idsScoreFn;
    idsScoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);
    idsScoreFn.ins = 4;
    idsScoreFn.del = 5;
    idsScoreFn.affineExtend = 20;
    idsScoreFn.substitutionPrior = 20;
    idsScoreFn.globalDeletionPrior = 13;
    if (criticalAccuracyFileName != "") {
      caFile << "read length";
      int i;
      for (i = 0; i < 100; i++) {
        caFile << " " << i;
      }
      caFile << endl;
    }
    while (reader.GetNext(smrtRead)) {

      vector<ReadInterval> subreadIntervals;
      int hqStart, hqEnd;
      int score;
      LookupHQRegion(smrtRead.zmwData.holeNumber, regionTable, hqStart, hqEnd, score);
      smrtRead.lowQualityPrefix = hqStart;
      smrtRead.lowQualitySuffix = smrtRead.length - hqEnd;
      CollectSubreadIntervals(smrtRead, &regionTable, subreadIntervals);
      if (hqStart > 0 and hqEnd > 0) {
        metrics.nSubreads.push_back(subreadIntervals.size());
      }
      int intvIndex;
      for (intvIndex = 1; intvIndex < subreadIntervals.size(); intvIndex++) {

        if (subreadIntervals[intvIndex-1].end < hqStart) {
          continue;
        }
        if (subreadIntervals[intvIndex].start > hqEnd) {
          continue;
        }
      
        DNALength prevStart = min(max(subreadIntervals[intvIndex-1].start, hqStart), hqEnd);
        DNALength prevEnd   = min(max(subreadIntervals[intvIndex-1].end, hqStart), hqEnd);
        DNALength curStart  = min(max(subreadIntervals[intvIndex].start, hqStart), hqEnd);
        DNALength curEnd    = min(max(subreadIntervals[intvIndex].end, hqStart), hqEnd);

        if (prevEnd - prevStart < minSubreadLength) {
          continue;
        }
        if (curEnd - curStart < minSubreadLength) {
          continue;
        }
      
        SMRTSequence prevSubread, curSubread, curSubreadRC;
        smrtRead.MakeSubreadAsReference(prevSubread, prevStart, prevEnd); 
        smrtRead.SetSubreadTitle(prevSubread, prevStart, prevEnd); 
        smrtRead.MakeSubreadAsReference(curSubread, curStart, curEnd); 
        smrtRead.SetSubreadTitle(curSubread, curStart, curEnd); 

        curSubread.MakeRC(curSubreadRC);
        smrtRead.SetSubreadTitle(curSubreadRC, curStart, curEnd); 
        Alignment alignmentGuide;
        int alignScore;
        alignScore = SDPAlign(prevSubread, curSubreadRC, distScoreFn, sdpTupleSize, 
                              5, 5, 1.30, 
                              alignmentGuide, Local, false, false);

        int b;
        for (b = 0; b < alignmentGuide.blocks.size(); b++) {
          alignmentGuide.blocks[b].qPos += alignmentGuide.qPos;
          alignmentGuide.blocks[b].tPos += alignmentGuide.tPos;
        }
        alignmentGuide.tPos = 0;
        alignmentGuide.qPos = 0;

        AlignmentCandidate<> refinedAlignment;
      
        AffineGuidedAlign(prevSubread, curSubreadRC, alignmentGuide, 
                          idsScoreFn, 15,
                          mappingBuffers, 
                          refinedAlignment, Global, false);


        string curStr, alnStr, prevStr;
        CreateAlignmentStrings(refinedAlignment, prevSubread.seq, curSubreadRC.seq, curStr, alnStr, prevStr, prevSubread.length, curSubread.length);

        vector<float> accBins;
        vector<int>   accLengths;
        accLengths.resize(100);
        fill(accLengths.begin(), accLengths.end(), 0);
      
      
        int i, w;
        //
        // This step can take a while, so only do it if specified as an option.
        //
        bool readIsGood = true;

        if (computeCriticalAccuracy) {

          if (criticalAccuracyLengthFraction > 0) {
            criticalAccuracyMinLength = max(prevSubread.length, curSubreadRC.length) * criticalAccuracyLengthFraction;
          }
          for (w = alnStr.size(); w > 20; w-=10) {
            int nMatch=0, nMismatch=0, nIns=0, nDel=0;
        
            //
            // Prepare the first window
            //
            for (i = 0; i < w; i++) {
              if (curStr[i] == '-') {
                nIns++;
              }
              else if (prevStr[i] == '-') {
                nDel++;
              }
              else if (toupper(curStr[i]) == toupper(prevStr[i])) {
                nMatch++;
              }
              else {
                nMismatch++;
              }
            }
            float acc;
            acc = nMatch / (1.0* nMatch + nMismatch + nIns + nDel);
            int index = min((int)(100*acc),99);
            if (accLengths[index] < w) {
              accLengths[index] = w;
            }

            for (i = 1; i + w < alnStr.size(); i++) {
              assert(nMatch >= 0);
              assert(nIns >= 0);
              assert(nDel >= 0);
              assert(nMismatch >= 0);
              //
              // Try a new window.
              //

              // remove previous
              if (curStr[i-1] == '-') {
                nIns--;
              }
              else if (prevStr[i-1] == '-'){ 
                nDel--;
              }
              else if (toupper(curStr[i-1]) == toupper(prevStr[i-1])) {
                nMatch--;
              }
              else {
                nMismatch--;
              }
              // add next base
              if (curStr[i+w-1] == '-') {
                nIns++;
              }
              else if (prevStr[i+w-1] == '-') {
                nDel++;
              }
              else if (toupper(curStr[i+w-1]) == toupper(prevStr[i+w-1])) {
                nMatch++;
              }
              else {
                nMismatch++;
              }
              acc = nMatch / (1.0* nMatch + nMismatch + nIns + nDel);
              int index = min((int)(100*acc),99);
              assert(index < accLengths.size());
              if (accLengths[index] < w) {
                accLengths[index] = w;
              }
            }
          }
        }

        if (computeCriticalAccuracy and criticalAccuracyFileName != "") {
          caFile << prevSubread.title << " " << prevSubread.length;
          for (i = 0; i < 100; i++) {
            caFile << " " << accLengths[i];
          }
          caFile << endl;
        }

        if (filterByCriticalAccuracy) {
          int a;
          for (a = 0; a < criticalAccuracyMinAccuracy and readIsGood; a++) {
            if (accLengths[a] >= criticalAccuracyMinLength) {
              readIsGood = false;
              metrics.criticalAccuracy.push_back(a);
              metrics.criticalLength.push_back(accLengths[a]);
            }
          }
        }

        if (filterByIndel) {
          int indelLength = 0;
          if ((indelLength = AlignmentHasIndel(refinedAlignment, indelLength))) {
            readIsGood = false;
            metrics.gapLength.push_back(indelLength);
          }
        }


        if (readIsGood == false) {
          if (filteredReadNameFileName != "") {
            filteredReadNameFile << prevSubread.title << endl;
            filteredReadNameFile << curSubreadRC.title << endl;
          }
          
          //
          // Later when the region table is written for this read,
          // setting these to 0 makes the read get skipped by
          // subsequenct alignment/analysis programs.
          //
          hqStart = 0;
          hqEnd   = 0;
        }

        
        //
        // Printing the alignemnts is optional, and either all
        // alignments are printed (just to peruse), or alignments for
        // gapped/bad reads are printed, to QC manually.
        //
        if (alignmentsFileName != "") {
          if (printOnlyFiltered == false or readIsGood == false ) {
            refinedAlignment.qName = prevSubread.title;
            refinedAlignment.tName = curSubreadRC.title;
            refinedAlignment.qLength = prevSubread.length;
            refinedAlignment.tLength = curSubreadRC.length;
            StickPrintAlignment(refinedAlignment, prevSubread, curSubreadRC, alnFile);
          }
        }

        curSubread.Free();
        curSubreadRC.Free();
        prevSubread.Free();
      }
      smrtRead.Free();


      if (printRegions) {
        vector<int> regionIndices;
        // 
        // Default parameters for this function collect all regions.
        //
        CollectRegionIndices(smrtRead, regionTable, regionIndices);
        int r;
        for (r = 0; r < regionIndices.size(); r++) {
          if (regionTable.GetType(regionIndices[r]) == HQRegion) {
            regionTable.SetStart(regionIndices[r], hqStart);
            regionTable.SetEnd(regionIndices[r], hqEnd);
            regionTable.SetScore(regionIndices[r], 0);
            break;
          }
        }
      }
    }

    if (metricsFileName != "") {
      metrics.Print(metricsFile);
    }

    if (printRegions) {
      int r;
      regionTableWriter.WriteRows(regionTable.table);
      regionTableWriter.Finalize(regionTable.columnNames,
                                 regionTable.regionTypes,
                                 regionTable.regionDescriptions,
                                 regionTable.regionSources);
      regionTableWriter.Close();
    }

  }
  
  //
  // These files are written to across all reads.  Don't close them
  // until the last bas.h5 file has been processed.
  //
  if (criticalAccuracyFileName != "") {
    caFile.close();
  }
  if (filteredReadNameFileName != "") {
    filteredReadNameFile.close();
  }
  if (metricsFileName != "") {
    metricsFile.close();
  }
}
