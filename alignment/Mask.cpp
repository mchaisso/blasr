#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFRegionTableReader.h"
#include "data/hdf/HDFRegionTableWriter.h"
#include "FASTAReader.h"
#include "algorithms/alignment/ExtendAlign.h"
#include "tuples/TupleMatching.h"
#include "CommandLineParser.h"
#include "utils.h"
#include "utils/RegionUtils.h"
#include "tuples/HashedTupleList.h"
#include "tuples/DNATuple.h"
#include "tuples/TupleMetrics.h"
#include "algorithms/alignment/KBandAlign.h"
#include "algorithms/alignment/IDSScoreFunction.h"
#include "algorithms/alignment/printers/StickAlignmentPrinter.h"
#include <set>
#include <vector>
#include <algorithm>

using namespace std;


vector<set<DNATuple> > adapterMatchActiveSet;

template<typename T_Iterator, typename T_Comparison>
bool GetPos(T_Iterator it, T_Iterator end, DNALength &pos, T_Comparison &comp) {
  if (it == end) {
    return false;
  }
  bool minSet = false;
  while (it != end) {
    if (minSet == false) {
      pos = (*it).pos;
      minSet = true;
    }
    else {
      pos = comp((*it).pos, (const DNALength) pos);
    }
    ++it;
  }
}

template<typename T>
class Min {
public:
  T operator()(const T&lhs, const T&rhs) const {
    return lhs < rhs ? lhs : rhs;
  }
};

template<typename T>
class Max {
public:
  T operator()(const T&lhs, const T&rhs) const {
    return lhs > rhs ? lhs : rhs;
  }
};
    
template <typename T, typename T_Iterable>
bool GetMinPos(T_Iterable &iterable, T &minPos) {
  Min<T> m;
  return GetPos(iterable.begin(), iterable.end(), minPos, m);
}

template <typename T, typename T_Iterable>
bool GetMaxPos(T_Iterable &iterable, T &maxPos) {
  Max<T> m;
  return GetPos(iterable.begin(), iterable.end(), maxPos, m);
}

template <typename T_clearable>
class ClearIt: public unary_function<T_clearable, void> {
public:
  void operator()(T_clearable &c) { c.clear(); }
};

//
// Search for adapters that flank the current insert.
//

void FindAdjacentAdapterIndices(RegionTable &regionTable, vector<int> &regionIndices,
                                int insertIndex, int &frontAdapterIndex, int &backAdapterIndex) {
  
  //
  // Look to see if there is an adapter in front of the current
  // insert.
  //

  if (insertIndex == 0) {
    // 
    // If the insert is the first region, there is no adapter in front
    // of it.
    //
    frontAdapterIndex = -1;
  }
  else {
    //
    // If the region before this one is an adapter, make sure it is
    // flanking.
    //
    if (regionTable.GetType(regionIndices[insertIndex-1]) == Adapter) {
      if (regionTable.GetEnd(regionIndices[insertIndex-1]) == regionTable.GetStart(regionIndices[insertIndex])) {
        frontAdapterIndex = regionIndices[insertIndex - 1];
      }
      else {
        frontAdapterIndex = -1;
      }
    }
  }

  //
  //  Look to see if htere is an adapter after the current insert.
  //
  if (insertIndex == regionIndices.size() - 1) {
    backAdapterIndex = -1;
  }
  else {
    if (regionTable.GetType(regionIndices[insertIndex+1]) == Adapter) {
      if (regionTable.GetEnd(regionIndices[insertIndex]) == regionTable.GetStart(regionIndices[insertIndex+1])) {
        backAdapterIndex = regionIndices[insertIndex+1];
      }
      else {
        backAdapterIndex = -1;
      }
    }
  }
}

class ModifiedIterval {
public:
  typedef enum E_ModificationType { Truncation, Replacement } ModificationType;
  vector<RegionAnnotation> regions;
};

void MapAdaptersToRead(SMRTSequence &seq, vector<FASTASequence> &adapters,
                       vector<HashedTupleList<DNATuple> > &adapterTupleLists,
                       IDSScoreFunction<DNASequence, FASTQSequence> &idsScoreFn,
                       vector<Arrow> &pathMat,
                       vector<int>   &scoreMat,
                       int maxScore,
                       int minWords=2, 
                       int k=10,
                       float delta=1.20) {
  vector<int> nMatches;
  vector<DNALength> first, last; 
  vector<set<DNATuple> > allMatches;
  int a;
  
  nMatches.resize(seq.length);
  first.resize(seq.length);
  last.resize(seq.length);
  allMatches.resize(seq.length);
  int nAdapters = adapters.size();

  TupleMetrics tm;
  tm.tupleSize = 5;

  for (a = 0; a < nAdapters; a++) {
    int prefixLength = min(adapters[a].length, seq.length);
    if (tm.tupleSize < prefixLength) {
      prefixLength -= tm.tupleSize;
    }
    else {
      prefixLength = 0;
    }
      
    fill(nMatches.begin(), nMatches.end(), 0);
    fill(first.begin(), first.end(), adapters[a].length);
    fill(last.begin(), last.end(), 0);
    for_each(allMatches.begin(), allMatches.end(), ClearIt<set<DNATuple> >());
    //
    // Add the matches for the first full length match on the sequence
    //
    bool firstIsSet = false;
    int i;
    for (i = 0; i < prefixLength; i++) {
      DNATuple tuple;
      tuple.FromStringLR(&seq.seq[i], tm);
      int hashPos, location;
      if (adapterTupleLists[a].Find(tuple, hashPos, location)) {
        tuple.pos = adapterTupleLists[a].hashTable[hashPos].tupleList[location].pos;
        adapterMatchActiveSet[a].insert(tuple);
      }
    }
    nMatches[prefixLength] = adapterMatchActiveSet[a].size();
    GetMinPos(adapterMatchActiveSet[a],first[prefixLength]);
    GetMaxPos(adapterMatchActiveSet[a],last[prefixLength]);
    //
    // Add the matches for the remainder of the sequence, kicking out any 
    // match at the beginning of the sequence.

    int fullLength = seq.length;
    if (seq.length < tm.tupleSize) {
      fullLength = 0;
    }
    else {
      fullLength = seq.length - tm.tupleSize + 1;
    }
      
    // 
    // Start at the length of the adapter in case the sequence is
    // short, so that it is skipped.
    //
    DNATuple trailingTuple, leadingTuple;
    int adapterLength = adapters[a].length;
    for (i = adapterLength; i < fullLength; i++) {
      trailingTuple.FromStringLR(&seq.seq[i-adapterLength], tm);
      leadingTuple.FromStringLR(&seq.seq[i], tm);
      multiset<DNATuple>::iterator trailingTupleIt = adapterMatchActiveSet[a].find(trailingTuple);
      while (trailingTupleIt != adapterMatchActiveSet[a].end()) {
        adapterMatchActiveSet[a].erase(trailingTupleIt);
        trailingTupleIt = adapterMatchActiveSet[a].find(trailingTuple);
      }
      int hashPos, location;
      if (adapterTupleLists[a].Find(leadingTuple, hashPos, location)) {
        while(location < adapterTupleLists[a].hashTable[hashPos].tupleList.size() and 
              adapterTupleLists[a].hashTable[hashPos].tupleList[location].tuple == leadingTuple.tuple) {
          leadingTuple.pos = adapterTupleLists[a].hashTable[hashPos].tupleList[location].pos;
          adapterMatchActiveSet[a].insert(leadingTuple);
          location++;
        }
      }
      //allMatches[i].insert(adapterMatchActiveSet[a].begin(), adapterMatchActiveSet[a].end());
      GetMinPos(adapterMatchActiveSet[a], first[i]);
      GetMaxPos(adapterMatchActiveSet[a], last[i]);
      nMatches[i] = adapterMatchActiveSet[a].size();
    }

    bool allMatchesExhausted = false;
    // 
    // Now greedily try and align portions of the read to adapters.
    //
    if (nMatches.size() == 0) {
      allMatchesExhausted = true;
    }

    while (!allMatchesExhausted) {
      int maxIndex, maxValue;
      vector<int>::iterator minIt = max_element(nMatches.begin(), nMatches.end());
      if (minIt == nMatches.end()) {
        allMatchesExhausted = true;
        break;
      }
      if ( (*minIt) <= minWords ) {
        allMatchesExhausted = true;
        break;
      }
        
      //
      // Found a match, try aligning the region it is in.
      // 
      SMRTSequence readSeq;
      DNALength matchPos = minIt - nMatches.begin();
      DNALength readStart, readEnd;
      //
      // Find the starting position of where the pairwise alignment
      // for where the adapter alignment should start.
      //
      if (matchPos >= (int)first[matchPos]*delta) {
        readStart = matchPos -first[matchPos]*delta;
      }
      else {
        readStart = 0;
      }

      DNALength trailingAdapter = adapters[a].length - last[matchPos];
      if ( matchPos + trailingAdapter * delta > seq.length) {
        readEnd = seq.length;
      }
      else {
        readEnd = matchPos + trailingAdapter * delta;
      }
      readSeq.ReferenceSubstring(seq, readStart, readEnd - readStart);
        
      Alignment alignment;
      int qvAwareScore;
      qvAwareScore = KBandAlign(readSeq, adapters[a], SMRTDistanceMatrix, 
                                idsScoreFn.ins, idsScoreFn.del, k,
                                scoreMat, pathMat,
                                alignment, idsScoreFn, Global);
      //        cout << "Pos: " << matchPos << " read pos: " << readStart << " n anchors: " << nMatches[matchPos] << endl;
      ComputeAlignmentStats(alignment, readSeq.seq, adapters[a].seq, idsScoreFn);
      if (alignment.score < maxScore) {
        alignment.qName = seq.title;
        alignment.tName = adapters[a].title;
              
        alignment.qLength = readSeq.length;
        alignment.tLength = adapters[a].length;
        //        StickPrintAlignment(alignment, readSe, adapters[a], cout, readInterval.start + readStart, 0);
        DNALength i;
        for (i = readStart; i < readEnd; i++) {
          nMatches[i] = 0;
          seq.seq[i] = 'N';
        }
      }
      else {
        nMatches[matchPos] = 0;
      }
    }
  }

  //
  // Reset match set for the next sequence
  //
  for (a = 0; a < nAdapters; a++) {
    adapterMatchActiveSet[a].clear();
  }

}


void PrintDot(int progress, ostream &out, int step = 1000, int line=50000) {
  if (progress % step == step-1 and progress > 0) {
    out << "." << flush;
    if (progress % line == line-1 and progress > 0) {
      out << endl;
    }
  }
}


int main(int argc, char* argv[]) {
  string sourceReadsFileName;
  string adaptersFileName;
  CommandLineParser clp;

  //
  // Parameters for masking anywhere in the read.
  //
  int minWords = 2;
  int k = 10;
  int maxScore = -50;
  float delta = 1.20;
  //
  // Alignment parameters if quality values are missing
  //

  int insertion = 3;
  int deletion  = 5;

  //
  // When flagged as true, only align to the ends of a read that are
  // adjacent to an adapter.
  //
  bool findFlushWithEnd = true;

  //
  // For testing, only mask a subset of reads.
  //
  vector<int> readIndices;

  // 
  // If set, write annotations of alignments to this file.
  //
  string annotationsFileName = "";
  vector<string> metricsList;
  int verbosity = 0;

  TupleMetrics tm;
  tm.tupleSize = 5;
  bool maskCCS = false;
  bool maskBest = false;
  bool writeRegionTable = true;
  string fastaOutputFileName = "";
  string outputRegionTableFileName = "";

  clp.RegisterStringOption("sourceReads", &sourceReadsFileName, "Input bas.h5 file.", true);
  clp.RegisterStringOption("adapters", &adaptersFileName, 
                           "Multi fasta file of adapters.  These should be short, < 100 base sequences.", true);
  clp.RegisterPreviousFlagsAsHidden();
  clp.RegisterIntOption("maxScore", &maxScore, "Mask adapters if the pairwise alignment score is less than maxScore."
                        "A score of -50 filters out many false positives.  More negative scores mean higher accuracy alignments.",
                        CommandLineParser::Integer);
  clp.RegisterStringOption("regionTable", &outputRegionTableFileName, "Write the new modified rgn.h5 table.");
  clp.RegisterStringOption("annotations", &annotationsFileName, "Specifies a file to write descriptions of adapter alignment hits."
                           "By default, the only description printed is the masked sequence with the higest scoring hit and "
                           "the score. Additional metrics may be added using the -metrics option. A description of "
                           "the columns is as follows: \n"
                           "\treadName\tThe name of the subread.\n" 
                           "\tfrontMaskIndex\tThe index of the adapter found at the front\n"
                           "\t\t\tof the read.  The index is zero based, and values >= \n"
                           "\t\t\t#adapters are reverse complement. A value of -1 means no \n"
                           "\t\t\tadapter found.\n"
                           "\tfrontTopScore\tThe score of the alignment to the adapter. \n"
                           "\t\t\t0 for no alignment.\n"
                           "\tbackMaskIndex,backTopScore The same as frontMaskIndex, \n"
                           "\t\t\tbut for the end of the read.\n"
                           "\tforwardAdapter\t1 if an adapter is present at the beginning\n"
                           "\t\t\tof this subread, 0 if not. Similar for endAdapter.\n"
                           "\tfrontAlignment\tThe pairwise alignment used to mask the\n"
                           "\t\t\tbeginning of the read. Similar for backAlignment.\n"); 

  clp.RegisterStringListOption("metrics", &metricsList, "Add metrics to print to annotations. The possible metrics are:\n"
                               "\tadapters\tPrint a column of 1/0 values indicating the\n"
                               "\t\t\tpresence of a forward and reverse hit.\n"
                               "\talignment\tPrint the pairwise alignment of the top hit (in \n\t\t\tone string).\n"
                               "\tcoordinates\tPrint the start/end of unmasked sequence for every\n"
                               "\t\t\tsubread or ccs sequence.\n"
                               "\t Values should be separated by spaces:\n"
                               "\t\t\t  \"-metrics adapterts coordinates\"\n");
  clp.RegisterFlagOption("anywhere", &findFlushWithEnd, "Allow adapters to be found anywhere, not just at beginning and ending of subreads.");
  clp.RegisterFlagOption("ccs", &maskCCS, "Mask circular consensus sequences instead of raw.");
  clp.RegisterFlagOption("best", &maskBest, "Mask ccs if it exists, and the full length (or longest) subread otherwise.  The type of read that is given will be annotated in a 'best' column of the annotation table.");
  clp.RegisterStringOption("fasta", &fastaOutputFileName, "Print the masked sequences to a file.");
  clp.RegisterIntOption("word", &tm.tupleSize, "Seed alignments on words of size 'word'.", CommandLineParser::PositiveInteger);
  clp.RegisterIntOption("minWords", &minWords, "Trigger alignments when minWords are matched.", CommandLineParser::PositiveInteger);
  clp.RegisterIntOption("k", &k, "K-band size for pairwise alignment for adapter alignment.", CommandLineParser::PositiveInteger);
  clp.RegisterFloatOption("delta", &delta, "Search length*delta for an adapter hit.", CommandLineParser::PositiveFloat);
  clp.RegisterIntListOption("readIndex", &readIndices, "Only mask the specified read indices.");
  clp.RegisterIntOption("v", &verbosity, "Debugging verbosity. 0 = none (default), 1 = little, 2 = show full alignments.", CommandLineParser::NonNegativeInteger);


  clp.SetProgramSummary("Filter non-pacbio adapters from reads.  If you have ligated extra adapters "
                        "to templates before sequencing, but do not want them to show up in alginments"
                        " use this to filter them out.  The original bas.h5 files are not modified. "
                        "Instead, a new region table file is created that has redefined insert intervals "
                        "that do not contain adapters.");
  clp.ParseCommandLine(argc, argv);
  
  ofstream annotationsFile;
  bool printAnnotations = annotationsFileName != "";
  
  if (maskCCS == true) {
    writeRegionTable = false;
  }

  if (annotationsFileName != "") {
    CrucialOpen(annotationsFileName, annotationsFile, std::ios::out);
  }
  vector<FASTASequence> adapters;
  FASTAReader fastaReader;
  fastaReader.Init(adaptersFileName);
  HDFBasReader reader;
  HDFBasReader ccsReader;

	reader.InitializeDefaultIncludedFields();  
  if (maskCCS) {
    reader.SetReadBasesFromCCS();
  }
  reader.Initialize(sourceReadsFileName);
  if (maskBest) {
    ccsReader.SetReadBasesFromCCS();
    ccsReader.Initialize(sourceReadsFileName);
  }

  fastaReader.ReadAllSequences(adapters);
  int nAdapters = adapters.size();
  adapters.resize(nAdapters*2);
  int i;
  for (i = 0; i < nAdapters; i++) {
    adapters[i].MakeRC(adapters[nAdapters+i]);
  }
  nAdapters = adapters.size();

  ofstream fastaOutFile;
  if (fastaOutputFileName != "") {
    CrucialOpen(fastaOutputFileName, fastaOutFile, std::ios::out);
  }

  SMRTSequence seq;

  //
  // This is for aligning to the engire read.  Not necessary until
  // that code is written.
  //

  vector<Arrow> pathMat;
  vector<int>   scoreMat;

  HDFRegionTableReader	regionTableReader;
  HDFRegionTableWriter  regionTableWriter;
	RegionTable regionTable;

  if (regionTableReader.Initialize(sourceReadsFileName) == 0) {
    cout <<" Could not find a region table in " << sourceReadsFileName << endl;
  }

  regionTableReader.ReadTable(regionTable);
  if (outputRegionTableFileName == "") {
    writeRegionTable = false;
  }

  if (writeRegionTable != false) {
    if (regionTableWriter.Initialize(outputRegionTableFileName) == 0) {
      cout <<"ERROR, could not initialize the modified region table " << outputRegionTableFileName << endl;
      exit(1);
    }
  }

  IDSScoreFunction<DNASequence, FASTQSequence> idsScoreFn;
  //
  // Initialize backup scores for alignment.
  //
  idsScoreFn.ins = insertion;
  idsScoreFn.del = deletion;
  idsScoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);


  SMRTSequence read;
  RegionType  adapterAndInsert[2] = {Adapter, Insert};
  int         numRegionTypes = 2;

  

  int m;
  set<string> availableMetrics, usedMetrics;
  availableMetrics.insert("alignment");
  availableMetrics.insert("adapters");
  availableMetrics.insert("coordinates");
  
  for (m = 0; m < metricsList.size(); m++ ) {
    if (availableMetrics.find(metricsList[m]) == availableMetrics.end()) {
      cout << "ERROR, the metric " << metricsList[m] << " is not valid." << endl;
      exit(1);
    }
    usedMetrics.insert(metricsList[m]);
  }

  //
  // Prepare the summary table if present.
  //

  if (printAnnotations) {
    annotationsFile << "readName\t";
    if (findFlushWithEnd) {
      annotationsFile << "\tfrontMaskIndex\tfrontTopScore\tbackMaskIndex\tbackTopScore";
    }
    if (usedMetrics.find("adapters") != usedMetrics.end()) {
      annotationsFile << "\tforwardAdapter\treverseAdapter";
    }
    if (usedMetrics.find("alignment") != usedMetrics.end()) {
      annotationsFile << "\tfrontAlignment\tbackAlignment" << endl;
    }
    if (usedMetrics.find("coordinates") != usedMetrics.end()) {
      annotationsFile << "\tstart\tend";
    }
    annotationsFile << endl;
  }
  if (maskBest) {
    annotationsFile << "\tread_type" << endl;
  }
  stringstream annotationsStream;
  int readIndex = 0;
  SMRTSequence ccsRead;

  while(reader.GetNext(read)) {

    bool bestIsCCS = false;
    if (maskBest) {
      ccsReader.GetNext(ccsRead);
      if (ccsRead.length > 0) {
        read = ccsRead;
        ccsRead.Free();
        bestIsCCS = true;
      }
    }
    PrintDot(readIndex, cout);
    ++readIndex;
    if (!maskCCS and 
        (readIndices.size() > 0 and find(readIndices.begin(), readIndices.end(), read.zmwData.holeNumber) == readIndices.end())) {
      continue;
    }
    vector<SMRTSequence> seqToMask;

    vector<int> frontAdapterIndices, backAdapterIndices, regionIndices;
    
    
    if (maskCCS or (maskBest and bestIsCCS )) {
      read.subreadStart = 0;
      read.subreadEnd = read.length;
      seqToMask.resize(seqToMask.size()+1);
      int last = seqToMask.size()-1;
      seqToMask[last] = read;
      seqToMask[last].subreadStart = 0;
      seqToMask[last].subreadEnd   = read.length;
      
      //
      // Add some faked adapter indices to signal that ther were
      // adapters on both sides of one of the inserts.
      //
      if (read.length > 0) {
        frontAdapterIndices.push_back(0);
        backAdapterIndices.push_back(0);
      }
      else {
        frontAdapterIndices.push_back(-1);
        backAdapterIndices.push_back(-1);
      }
    }
    else {
      //    cout << read.title << endl;
      int regionLowIndex, regionHighIndex;
      FindRegionIndices(read.zmwData.holeNumber, &regionTable, regionLowIndex, regionHighIndex);
      vector<int> adapterAndInsertIndices;
      CollectRegionIndices(read, regionTable, adapterAndInsertIndices, adapterAndInsert, numRegionTypes);
      SortRegionIndicesByStart(regionTable, adapterAndInsertIndices);
      int regionTableIndex;
      for (regionTableIndex = regionLowIndex; regionTableIndex < regionHighIndex; ++regionTableIndex) {
        int r;
        if (regionTable.GetType(regionTableIndex) != Insert) {
          // 
          // This is not an insert, therefore it shouldn't be modified,
          // so just write it out.
          //
          continue;
        }
      
        //
        // Check to see if this is an adapter/insert that is being aligned.
        //
        for (r = 0; r < adapterAndInsertIndices.size(); r++) {
          if (adapterAndInsertIndices[r] == regionTableIndex) {
            break;
          }
        }
      
        // Error check -- the adaper/insert list is built from the
        // region table, so it should be found.  
        if (r == adapterAndInsertIndices.size()) {
          cout << "ERROR! A table of regions was created that does not match the original." << endl;
          cout << "This is a bug." << endl;
          exit(1);
        }

        //
        // If this is an insert, look to see if it is flanked by
        // adapters.
        //
        int regionIndex = adapterAndInsertIndices[r];
        if (regionTable.GetType(regionIndex) == Insert) {

          int frontAdapterIndex = -1, backAdapterIndex = -1;
          if (findFlushWithEnd == true) {
            FindAdjacentAdapterIndices(regionTable, adapterAndInsertIndices,
                                       r, frontAdapterIndex, backAdapterIndex);
          }

          //
          // If there are no adapters, and you are looking for match adjacent to an adapter, continue.
          //
          ReadInterval  readInterval(regionTable.GetStart(regionIndex),
                                     regionTable.GetEnd(regionIndex)); 
       
          SMRTSequence seq;
          if (readInterval.end <= readInterval.start) {
            continue;
          }
          read.MakeSubreadAsReference(seq, readInterval.start, readInterval.end);
          read.SetSubreadTitle(seq, readInterval.start, readInterval.end);
          seqToMask.resize(seqToMask.size() + 1);
          int last = seqToMask.size()-1;
          seqToMask[last] = seq;
          seqToMask[last].subreadStart = readInterval.start;
          seqToMask[last].subreadEnd   = readInterval.end;
          frontAdapterIndices.push_back(frontAdapterIndex);
          backAdapterIndices.push_back(backAdapterIndex);
          regionIndices.push_back(regionIndex);
        }
      }
    }

    if (maskBest and ccsRead.length == 0) {
      //
      // Need to only mask the best sequence for this read.  The best is:
      //  1. The longest full-length sequence (if it exists)
      //  2. The longest sequence otherwise.
      //
      vector<bool> hasBothAdapters;
      hasBothAdapters.resize(seqToMask.size());
      fill(hasBothAdapters.begin(), hasBothAdapters.end(), false);
      int s;
      int longestWithAdaptersIndex = 0;
      int longestWithAdaptersLength = 0;
      int longestIndex = 0;
      int longestLength = 0;
      int numWithBothAdapters = 0;
      for (s = 0; s < seqToMask.size(); s++) {
        if (frontAdapterIndices[s] != -1 and
            backAdapterIndices[s]  != -1) {
          hasBothAdapters[s] = true;
          numWithBothAdapters++;
          if (longestWithAdaptersLength < seqToMask[s].length) {
            longestWithAdaptersLength = seqToMask[s].length;
            longestWithAdaptersIndex  = s;
          }
        }
        if (longestLength < seqToMask[s].length) {
          longestLength = seqToMask[s].length;
          longestIndex  = s;
        }
      }
      int bestIndex = 0;
      if (numWithBothAdapters > 0) {
        bestIndex = longestWithAdaptersIndex;
      }
      else {
        bestIndex = longestIndex;
      }
      
      if (seqToMask.size() > 0) {
        if (bestIndex < seqToMask.size()-1) {
          seqToMask.erase(seqToMask.begin() + bestIndex + 1, seqToMask.end());
          frontAdapterIndices.erase(frontAdapterIndices.begin() + bestIndex + 1, frontAdapterIndices.end());
          backAdapterIndices.erase(backAdapterIndices.begin() + bestIndex + 1, backAdapterIndices.end());
        }
        if (bestIndex > 0) {
          seqToMask.erase(seqToMask.begin(), seqToMask.begin() + bestIndex);
          frontAdapterIndices.erase(frontAdapterIndices.begin(), frontAdapterIndices.begin() + bestIndex);
          backAdapterIndices.erase(backAdapterIndices.begin(), backAdapterIndices.begin() + bestIndex);
        }
      }
    }

    int seqIndex;
    for (seqIndex = 0; seqIndex < seqToMask.size(); seqIndex++) {

      //
      // If the alignment mode is masking adjacent to adapters, and
      // there are no adjacent adapters, nothing can be done, so
      // continue.
      //
      int frontAdapterIndex = frontAdapterIndices[seqIndex];
      int backAdapterIndex  = backAdapterIndices[seqIndex];
      if (findFlushWithEnd == true) {
        //
        // Simple case, just align the adapter to the prefix and suffix of the read.
        //
        vector<int> frontAlignScores, backAlignScores;
        vector<int> frontAlignEnd, backAlignStart;
        vector<string> frontAlignStrings, backAlignStrings;
        frontAlignScores.resize(adapters.size());
        backAlignScores.resize(adapters.size());
        frontAlignEnd.resize(adapters.size());
        backAlignStart.resize(adapters.size());
        frontAlignStrings.resize(adapters.size());
        backAlignStrings.resize(adapters.size());
        fill(frontAlignScores.begin(), frontAlignScores.end(), 0);
        fill(frontAlignEnd.begin(), frontAlignEnd.end(), 0);
        fill(backAlignScores.begin(), backAlignScores.end(), 0);
        int a;
        for (a = 0; a < adapters.size(); a++) {

          if (frontAdapterIndex != -1) {
            // Align to theprefix.
            SMRTSequence prefix;
            seqToMask[seqIndex].MakeSubreadAsReference(prefix, 0, min(adapters[a].length * delta, (float)seq.length));
            Alignment prefixAlignment;
            int qvAwareScore;
            
            qvAwareScore = KBandAlign(prefix, adapters[a], SMRTDistanceMatrix, 
                                      insertion, deletion, k,
                                      scoreMat, pathMat,
                                      prefixAlignment, idsScoreFn, Global);            
            ComputeAlignmentStats(prefixAlignment, prefix.seq, adapters[a].seq, idsScoreFn);

            frontAlignScores[a] = prefixAlignment.score;
            frontAlignEnd[a]    = prefixAlignment.GenomicTEnd();
            if (usedMetrics.find("alignment") != usedMetrics.end()) {
              string qStr, alnStr, tStr;
              CreateAlignmentStrings(prefixAlignment, prefix.seq, adapters[a].seq, tStr, alnStr, qStr);
              frontAlignStrings[a] = "'" + qStr + " " + alnStr + " " + tStr + "'";
            }

              
            if (prefixAlignment.score < maxScore) {
              prefixAlignment.qName = seqToMask[seqIndex].title;
              prefixAlignment.tName = adapters[a].title;
              prefixAlignment.qLength = read.length;
              prefixAlignment.tLength = adapters[a].length;
              if (verbosity > 1) {
                cout << "prefix adapter index " << a << endl;
                StickPrintAlignment(prefixAlignment, prefix, adapters[a], cout, seqToMask[seqIndex].subreadStart, 0);
              }
            }
          }

          if (backAdapterIndex != -1) {
            SMRTSequence suffix;
            seqToMask[seqIndex].MakeSubreadAsReference(suffix, seqToMask[seqIndex].length - min(adapters[a].length * delta, (float)seqToMask[seqIndex].length), seqToMask[seqIndex].length);
            Alignment suffixAlignment;
            int qvAwareScore;
            qvAwareScore = KBandAlign(suffix, adapters[a], SMRTDistanceMatrix, 
                                      insertion, deletion, k,
                                      scoreMat, pathMat,
                                      suffixAlignment, idsScoreFn, Global);            
            ComputeAlignmentStats(suffixAlignment, suffix.seq, adapters[a].seq, idsScoreFn);

            backAlignScores[a] = suffixAlignment.score;
            backAlignStart[a]  = seqToMask[seqIndex].length - (suffix.length - suffixAlignment.GenomicTBegin());
            if (usedMetrics.find("alignment") != usedMetrics.end()) {
              string qStr, alnStr, tStr;
              CreateAlignmentStrings(suffixAlignment, suffix.seq, adapters[a].seq, tStr, alnStr, qStr);
              backAlignStrings[a] = "'" + qStr + " " + alnStr + " " + tStr + "'";
            }
            if (suffixAlignment.score < maxScore) {
              suffixAlignment.qName = seqToMask[seqIndex].title;
              suffixAlignment.tName = adapters[a].title;
              
              suffixAlignment.qLength = read.length;
              suffixAlignment.tLength = adapters[a].length;
              if (verbosity > 1) {
                cout << "suffix adapter index " << a << endl;
                StickPrintAlignment(suffixAlignment, suffix, adapters[a], cout, seqToMask[seqIndex].subreadStart + suffix.subreadStart, 0);
              }
            }
          }
        }

        vector<int>::iterator f;
        vector<int>::iterator r;
        int fi = -1, ri = -1;
        int optFScore=0, optRScore=0;
        if (frontAdapterIndex != -1) {
          int p;
          f = min_element(frontAlignScores.begin(), frontAlignScores.end());
          if (*f < maxScore) {
            fi = f - frontAlignScores.begin();
            if (writeRegionTable) {
              regionTable.SetStart(regionIndices[seqIndex], frontAlignEnd[fi]);
            }
            optFScore = *f;
          }
        }
        if (backAdapterIndex != -1) {
          int p;
          r = min_element(backAlignScores.begin(), backAlignScores.end());
          if (*r < maxScore) {
            ri = r - backAlignScores.begin();
            if (writeRegionTable) {
              regionTable.SetEnd(regionIndices[seqIndex], backAlignStart[ri]);
            }
            optRScore = *r;
          }
        }
        if (printAnnotations) {
          annotationsStream << seqToMask[seqIndex].title;
          annotationsStream << "\t" << fi << "\t" << optFScore << "\t" << ri << "\t" << optRScore;
          if (usedMetrics.find("adapters") != usedMetrics.end()) {
            annotationsStream << "\t" << (int) (frontAdapterIndex != -1) << "\t" <<(int) (backAdapterIndex != -1) <<  "\t";
          }
          if (usedMetrics.find("alignment") != usedMetrics.end()) {
            if (fi != -1) {
              annotationsStream << "\t" << frontAlignStrings[fi];
            }
            else {
              annotationsStream << "\tNA";
            }
            if (ri != -1) {
              annotationsStream << "\t" << backAlignStrings[ri];
            }
            else {
              annotationsStream << "\tNA";
            }
          }
          if (usedMetrics.find("coordinates") != usedMetrics.end()) {
            if (fi != -1) {
              annotationsStream << "\t" << frontAlignEnd[fi];
            }
            else {
              annotationsStream << "\t0";
            }
            if (ri != -1) {
              annotationsStream << "\t" << backAlignStart[ri];
            }
            else {
              annotationsStream << "\t" << seqToMask[seqIndex].length;
            }
          }
          if (maskBest) {
            if (bestIsCCS) {
              annotationsStream << "\tccs";
            }
            else {
              annotationsStream << "\traw";
            }
          }
        }

        if (fastaOutputFileName != "" and seqToMask[seqIndex].length > 0 ) {
          int readStart, readEnd;
          if (fi == -1) {
            readStart = 0;
          }
          else {
            readStart = frontAlignEnd[fi];
          }
          if (ri == -1) {
            readEnd = seqToMask[seqIndex].length;
          }
          else {
            readEnd = backAlignStart[ri];
          }
          FASTASequence subseq;
          if (readStart < readEnd) {
            subseq.CopySubsequence(seqToMask[seqIndex], readStart, readEnd);
            subseq.CopyTitle(seqToMask[seqIndex].title);
            subseq.PrintSeq(fastaOutFile);
            subseq.Free();
          }
        }
      }
      else {
        //
        // Search for a match anywhere in the read.  This is much
        // more slow than solely at ends, since each adapter must be
        // aligned to each position in the read.  This is sped up
        // using keyword matching.
        //
        cout << "This isn't finshed yet!" << endl;
        exit(0);
        assert(findFlushWithEnd == false); //check my logic

      }
      if (printAnnotations) {
        annotationsFile << annotationsStream.str() << endl;
        annotationsStream.str("");
      }
    }
  }
  cout << endl;
  
  if (findFlushWithEnd == true) {
    if (writeRegionTable) {
      regionTableWriter.WriteRows(regionTable.table);
      regionTableWriter.Finalize(regionTable.columnNames, 
                                 regionTable.regionTypes, 
                                 regionTable.regionDescriptions, 
                                 regionTable.regionSources);
      regionTableWriter.Close();
    }
  }

  //
  // Tidy up.
  //
  regionTableReader.Close();
  reader.Close();

}
