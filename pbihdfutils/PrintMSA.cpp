#include "data/hdf/HDFCmpFile.h"
#include "datastructures/alignment/CmpFile.h"
#include "CommandLineParser.h"
#include "datastructures/alignment/ByteAlignment.h"
#include "algorithms/alignment/printers/StickAlignmentPrinter.h"
#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFBasWriter.h"
#include "data/hdf/HDFRegionTableReader.h"
#include "data/hdf/HDFRegionTableWriter.h"
#include "utils/RegionUtils.h"
#include "datastructures/reads/RegionTable.h"
#include "files/ReaderAgglomerate.h"
#include "utils/FileOfFileNames.h"
#include "datastructures/alignment/ByteAlignment.h"
#include "Types.h"
#include "datastructures/matrix/Matrix.h"
#include <string>
#include <ext/hash_map>
#include <map>


#include "InsertedString.h"
#include "Profile.h"
#include "MSAUtils.h"

typedef map<string, int> MovieNameToArrayIndex;
using namespace std;

typedef vector<unsigned char> ByteAlignment;


void FormLabelMargin(int marginWidth, char label, string &margin) {
  if (marginWidth < 2) {
    marginWidth = 2;
  }
  margin.resize(marginWidth);
  margin[0] = ' ';
  margin[1] = label;
  fill(margin.begin() + 2, margin.end(), ' ');
}

void CondenseCharacterVector(vector<UChar> &fv, vector<char> &cv, int scale=5) {
  cv.resize(fv.size());
  int i;
  static const char charArray[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
  for (i = 0; i < fv.size(); i++) {
    if (!isnan(fv[i])) {
      int value  = (int)(fv[i] / scale);
      if (value > 9) { value = 9; }
      cv[i] = charArray[value];
    }
    else {
      cv[i] = 'N';
    }
  }
}


void StoreInsertionQVs(InsertedStringList &insertions, vector<char> &valueStr) {
  int i;
  for (i = 0; i < insertions.size(); i++) {
    int substrLen = insertions[i].insSeq.size();
    int substrPos = insertions[i].alnPos;
    // Modify the string to be from the valuestr.
    insertions[i].insQVSeq = "";
    insertions[i].insQVSeq.insert(0, &valueStr[substrPos], substrLen);
  }
}


void PrintInsertions(InsertedStringList &insertions,
                     int refStart, int readStart, int lineWidth, 
                     string &insMarginStr, string &insQVMarginStr,
                     bool printQVs,
                     ostream &out) {
  //
  // Assume insertions are already sorted by pos.
  //
  bool insertionsRemain = true;
  int curPos = refStart;
  int insIndex = 0;
  int nPrinted = 0;
  int curOffset;
  bool onNewLine = true;

  stringstream insertionQVStream;

  while (nPrinted < insertions.size()) {
    //
    // Do a cyclic search through the insertions to 
    //
    while ((readStart + insertions[insIndex].pos) < curPos or
           insertions[insIndex].printed == true) { 
      insIndex++;
      if (insIndex == insertions.size()) {
        out << endl;
        if (printQVs == true) {
          out << insQVMarginStr << insertionQVStream.str() << endl;
          insertionQVStream.str("");
        }
        insIndex = 0;
        curPos = refStart;
        onNewLine = true;
      }
    }
    curOffset = (readStart + insertions[insIndex].pos) - curPos;
    int i;
    if (onNewLine) { cout << insMarginStr; }
    onNewLine = false;
    //
    //  Print the insertion.
    for (i = 0; i < curOffset; i++) { out << " "; };
    out << insertions[insIndex].insSeq << " ";
    
    // Print the quality of the insertion.
    if (printQVs) {
      for (i = 0; i < curOffset; i++) { insertionQVStream << " "; };
      insertionQVStream << insertions[insIndex].insQVSeq << " ";
    }      

    insertions[insIndex].printed = true;
    nPrinted++;
    curPos = curPos + insertions[insIndex].insSeq.size() + 1 + curOffset;
  }
  if (!onNewLine) {
    out << endl;
  }
  if (printQVs == true and !onNewLine) {
    out << insQVMarginStr << insertionQVStream.str() << endl;
    insertionQVStream.str("");
  }

}


  
void MakeMatchesDot(string &ref, string &query) {
  int i;
  if (ref.size() != query.size()) return;
  for (i = 0; i < ref.size(); i++) {
    if (ref[i] == query[i]) {
      query[i] = '.';
    }
  }
}

void MakeMismatchesLowerCase(string &ref, string &query) {
  int i;
  if (ref.size() != query.size()) return;
  for (i = 0; i < ref.size(); i++) {
    if (ref[i] != query[i]) {
      query[i] = tolower(query[i]);
    }
  }
}

void SetRefString(DNASequence &genome, int pos, int width, string &refStr) {
  refStr.resize(width*2+1);
  fill(refStr.begin(), refStr.end(), ' ');
  int p;
  int begin, end;
  begin = (pos >= width) ? pos - width : 0;
  end   = (pos + width + 1 <= genome.length) ? pos + width + 1 : genome.length;
  
  for (p = pos; p >= begin; p-- ) {
    refStr[p - pos + width] = genome.seq[p];
    //
    // Check to see if the next iteration goes past the beginning of the genome.
    //
  }

  for (p = pos+1; p < end; p++) {
    refStr[p - pos + width] = genome.seq[p];
  }
}

int FindPosition(vector<int> &posMap, int pos, int &index) {
  index = 0;
  // 
  // Handle some boundary conditions.
  //
  if (posMap.size() == 0) {
    return 0;
  }
  if (pos < 0 or pos > posMap[posMap.size()-1]) {
    return 0;
  }
  for (index = 0; index < posMap.size(); index++) {
    if (posMap[index] == pos) return 1;
  }
  return 0;
}

void CopyFieldString(ByteAlignment &aln, vector<char> &fieldChars, int begin, int end, string &fieldStr) {
  int i;
  fieldStr = "";
  for (i = begin; i < end; i++) {
    if (RefChar[aln[i]] != ' ') {
      if (QueryChar[aln[i]] != ' ') 
        fieldStr.append(&fieldChars[i], 1);
      else 
        fieldStr.push_back(' ');        
    }
  }
}

/*
 * This function does two things, which is bad software engineering
 * but fine for code that should not last too long.  It fills the
 * query string with characters from the byte alignment, and finds the
 * positions in the byte alignment that allow the position in the
 * alignemnt at centerIndex to be printed at the center of the
 * alignment.
 * It now does three things.  The third result of this function is to
 * store the padding printed before the read characters start.
 */

void SetAlignString(ByteAlignment &aln, vector<int> &pos, int centerIndex, int width, string &alnStr, int &beginIndex, int &endIndex, string &padding) {
  if (width <= 0) return;
  if (aln.size() == 0) return;
  if (pos.size() == 0) return;
  // If this doesn't overlap, bail
  alnStr.resize(width*2+1);
  fill(alnStr.begin(), alnStr.end(), ' ');
  if (centerIndex + width < pos[0]) return;
  // pos is ths pos in the ref, index is the index in the aln string
  int centerPos = pos[centerIndex];
  int beginPos, endPos;
  beginPos = -1;
  endPos   = -1;
  int p;
  int bp;
  p = width;
  int nRefBases = 0;
  bp = centerIndex;
  //
  // Store aligment string from center to left.
  //
  while (p >= 0 and bp >= 0 and nRefBases < width + 1) {
    if (RefChar[aln[bp]] != ' ') {
      alnStr[p] = QueryChar[aln[bp]];
      p--;
      nRefBases++;
    }
    bp--;
  }
  //
  // Store the padding for use in printing other fields.
  //
  if (p > 0) {
    padding.resize(p);
    std::fill(padding.begin(), padding.end(), ' ');
  }

  //
  // Store the alignment string from the center to the right.
  //
  beginIndex = bp + 1;
  p = 0;
  nRefBases = 0;
  bp = centerIndex + 1;
  while (p < width and bp < aln.size() and nRefBases < width) {
    if (RefChar[aln[bp]] != ' ') {
      alnStr[p + width + 1] = QueryChar[aln[bp]];
      p ++;
    }
    bp++;
  }
  endIndex = bp;
}

int ComputeOffsetToRefPosition(ByteAlignment &aln, DNALength alnStartPos, DNALength refPos, DNALength &alnIndex, DNALength &refOffset) {
  alnIndex = 0;
  refOffset = 0;
  while(alnIndex < aln.size() and refOffset + alnStartPos != refPos) {
    // move past the current position.
    if (RefChar[aln[alnIndex]] != ' ') {
      refOffset++;
    }
    alnIndex++;
  }
  if (alnIndex < aln.size()) {
    return 1;
  }
  else {
    return 0;
  }
}


int main(int argc, char* argv[]) {
  CommandLineParser clp;
  string cmpH5FileName, genomeFileName;
  string basH5OutFileName, readsFileName;
  string posFileName = "";
  basH5OutFileName = "";
  readsFileName    = "";
  vector<int> positions;
  int width = 20;
  bool printMSA = false;
  bool showDeletionQV = false;
  bool showInsertionQV= false;
  bool showSubstitutionQV=false;
  bool showMergeQV    = false;
  bool insertionsOnly = false;
  bool hideInsertions = false;
  bool miscallView    = false;
  bool printProfile   = false;
  int  scale = 5;
  int  marginSize = 5;
  string alnStrMargin, insMargin, delQVMargin, subQVMargin, insQVMargin, mergeQVMargin;
  stringstream verboseHelpStream, helpStream;
  verboseHelpStream << "printmsa is a utility to either view the msa centered at a position" << endl
                    << "in a base.h5 file, or print all reads overlapping a position." << endl;
  clp.SetVerboseHelp(verboseHelpStream.str());
  //  clp.SetHelp(verboseHelpStream.str());
  clp.RegisterStringOption("cmpH5FileName", &cmpH5FileName, "The input hdf file");
  clp.RegisterStringOption("genome", &genomeFileName, "Name of the reference reads were aligned to (1 contig for now)");
  clp.RegisterPreviousFlagsAsHidden();
  clp.RegisterIntListOption("positions",   &positions, "Positions to add");
  clp.RegisterStringOption("posFile",      &posFileName, "File of positions to add");
  clp.RegisterIntOption("width",           &width, "Width of region to print", CommandLineParser::PositiveInteger);
  clp.RegisterFlagOption("msa",            &printMSA, "Print the MSA at this position");
  clp.RegisterFlagOption("profile",        &printProfile, "Print a profile of nucleotides at each position");
  clp.RegisterIntOption("scale",           &scale, "Scale quality values by 'scale'", CommandLineParser::PositiveInteger);
  clp.RegisterIntOption("margin",          &marginSize, "The size of the margin to print.", CommandLineParser::PositiveInteger);
  clp.RegisterFlagOption("showDelQV",      &showDeletionQV, "Show deletion qvs below each read.");  
  clp.RegisterFlagOption("showInsQV",      &showInsertionQV, "Show insertion qvs below each read.");  
  clp.RegisterFlagOption("showSubQV",      &showSubstitutionQV, "Show substitution qvs below each read.");  
  clp.RegisterFlagOption("showMergeQV",      &showMergeQV, "Show merge qvs below each read.");  
  clp.RegisterFlagOption("insertionsOnly", &insertionsOnly, "Only print inserted strings.");
  clp.RegisterFlagOption("hideInsertions", &hideInsertions, "Hide inserted strings.");
  clp.RegisterFlagOption("miscallView",    &miscallView, "Show matches as '.'");
  clp.RegisterStringOption("readsFile",    &readsFileName, "Input reads file or fofn");
  clp.RegisterStringOption("baseFile",     &basH5OutFileName, "Print reads to file.");

  clp.ParseCommandLine(argc, argv);
  FormLabelMargin(marginSize, ' ', alnStrMargin);
  FormLabelMargin(marginSize, 'I', insQVMargin);
  FormLabelMargin(marginSize, 'i', insMargin);
  FormLabelMargin(marginSize, 'd', delQVMargin);
  FormLabelMargin(marginSize, 's', subQVMargin);
  FormLabelMargin(marginSize, 'm', mergeQVMargin);
  vector<vector<int> > readsByPosition;
  
  vector<HDFBasReader > readers;
  vector<HDFRegionTableReader> regionTableReaders;
  vector<set<int> > printedRegionsByHoleNumber;

  FASTASequence genome;
  ReaderAgglomerate reader;
  reader.Initialize(genomeFileName);
  reader.GetNext(genome);

  readsByPosition.resize(genome.length);

  if (posFileName != "") {
    ifstream posIn;
    CrucialOpen(posFileName, posIn, std::ios::in);
    int pos;
    while(posIn >> pos) {
      positions.push_back(pos);
    }
    posIn.close();
  }

  CmpFile cmpFile;
	
	/*
	 * These readers pull information from the same pls file.
	 */
	HDFCmpFile<CmpAlignment> cmpReader;
  vector<RegionTable> regionTables;
  vector<string> readsFileNames;

  if (readsFileName != "") {
    FileOfFileNames::StoreFileOrFileList(readsFileName, readsFileNames);
  }
	MovieNameToArrayIndex movieNameToReaderIndex;
  int readsFileIndex;
  readers.resize(readsFileNames.size());
  regionTableReaders.resize(readsFileNames.size());
  regionTables.resize(readsFileNames.size());
  printedRegionsByHoleNumber.resize(readsFileNames.size());
  for (readsFileIndex = 0; readsFileIndex < readsFileNames.size(); readsFileIndex++) {
    readers[readsFileIndex].InitializeDefaultIncludedFields();
    readers[readsFileIndex].Initialize(readsFileNames[readsFileIndex]);
    readers[readsFileIndex].PrepareForRandomAccess();
		movieNameToReaderIndex[readers[readsFileIndex].GetMovieName()] = readsFileIndex;
    
    if (readers[readsFileIndex].hasRegionTable == false) {
      cout << "ERROR! printmsa must be ran on bas.h5 files that contain " << endl
           << "a region table." << endl;
      exit(1);
    }
    regionTableReaders[readsFileIndex].Initialize(readsFileNames[readsFileIndex]);
    regionTableReaders[readsFileIndex].ReadTable(regionTables[readsFileIndex]);
  }

  HDFBasWriter basWriter;
  HDFRegionTableWriter regionWriter;
  if (basH5OutFileName != "") {
    string movieName = "fake_movie";
    if (readers.size() > 0) {
      movieName = readers[0].GetMovieName();
      map<string,bool>::iterator includedFieldIt;
      for (includedFieldIt = readers[0].includedFields.begin();
           includedFieldIt != readers[0].includedFields.end();
           ++includedFieldIt) {
        if (includedFieldIt->second == true) {
          basWriter.IncludeField(includedFieldIt->first);
        }
      }
      basWriter.Initialize(basH5OutFileName, movieName);
      assert(regionTables.size() > 0);
      regionWriter.Initialize(basWriter.pulseDataGroup);
    }
  }

  if (showDeletionQV) {
    cmpReader.IncludeField("DeletionQV");
  }
  if (showInsertionQV) {
    cmpReader.IncludeField("InsertionQV");
  }
  if (showSubstitutionQV) {
    cmpReader.IncludeField("SubstitutionQV");
  }
  if (showMergeQV) {
    cmpReader.IncludeField("MergeQV");
  }
	if (cmpReader.Initialize(cmpH5FileName, H5F_ACC_RDWR) == 0) {
		cout << "ERROR, could not open the cmp file." << endl;
		exit(1);
	}
	cout << "Reading cmp file." << endl;
	cmpReader.Read(cmpFile);

  cout << "Building alignment table." << endl;
  int alignmentIndex;
  for (alignmentIndex = 0; alignmentIndex < cmpFile.alnInfo.alignments.size(); alignmentIndex++) {
    int refBegin = cmpFile.alnInfo.alignments[alignmentIndex].GetRefStart();
    int refEnd   = cmpFile.alnInfo.alignments[alignmentIndex].GetRefEnd();

    int genomePos;
    int p;
    for (p = 0; p < positions.size(); p++) {
      if (refBegin <= positions[p] and refEnd > positions[p]) {
        for (genomePos = refBegin; genomePos < refEnd; genomePos++) {
          readsByPosition[genomePos].push_back(alignmentIndex);
        }
        break;
      }
    }
  }

  SMRTSequence read;

  int p;
  for (p = 0; p < positions.size(); p++) {
    int a;
    int centerPos = positions[p];
    cout << "position " << centerPos << endl;
    string refStr;
    if (printMSA) {
      SetRefString(genome, centerPos, width, refStr);
    }
    cout << alnStrMargin << refStr << endl;
    Profile profile(width*2+1);
    for (a = 0; a < readsByPosition[positions[p]].size(); a++) {
      alignmentIndex = readsByPosition[positions[p]][a];
      int readStart, readEnd, refStart, refEnd;
      readStart = cmpFile.alnInfo.alignments[alignmentIndex].GetQueryStart();
      readEnd   = cmpFile.alnInfo.alignments[alignmentIndex].GetQueryEnd();
      refStart  = cmpFile.alnInfo.alignments[alignmentIndex].GetRefStart();
      refEnd    = cmpFile.alnInfo.alignments[alignmentIndex].GetRefEnd();
      int holeNumber = cmpFile.alnInfo.alignments[alignmentIndex].GetHoleNumber();
      int movieId    = cmpFile.alnInfo.alignments[alignmentIndex].GetMovieId();
			int refGroupId = cmpFile.alnInfo.alignments[alignmentIndex].GetRefGroupId();
			int refGroupIndex = cmpReader.refGroupIdToArrayIndex[refGroupId];
            int alnGroupId = cmpFile.alnInfo.alignments[alignmentIndex].GetAlnGroupId();
            string readGroupName = cmpReader.alnGroupIdToReadGroupName[alnGroupId];
            int readGroupIndex   = cmpReader.refAlignGroups[refGroupIndex]->experimentNameToIndex[readGroupName];
      int refStrand      = cmpFile.alnInfo.alignments[alignmentIndex].GetRCRefStrand();


      if (printMSA) {
			//
			// Read the alignment string.  All alignments 
			//
        int offsetBegin = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetBegin();
        int offsetEnd   = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetEnd();
		
        int alignedSequenceLength = offsetEnd - offsetBegin;
        string   alignedSequence;
        ByteAlignment byteAlignment, byteAlignmentRC;
        vector<char> insQVChars, delQVChars, subQVChars, mergeQVChars;

        if (alignedSequenceLength >= 0) {
          alignedSequence.resize(alignedSequenceLength);
          byteAlignment.resize(alignedSequenceLength);
        }
        byteAlignment = cmpFile.alnInfo.alignments[alignmentIndex].alignmentArray;

        CondenseCharacterVector(cmpFile.alnInfo.alignments[alignmentIndex].fields["InsertionQV"],    insQVChars, scale);
        CondenseCharacterVector(cmpFile.alnInfo.alignments[alignmentIndex].fields["DeletionQV"],     delQVChars, scale);
        CondenseCharacterVector(cmpFile.alnInfo.alignments[alignmentIndex].fields["SubstitutionQV"], subQVChars, scale);
        CondenseCharacterVector(cmpFile.alnInfo.alignments[alignmentIndex].fields["MergeQV"], mergeQVChars, scale);

        if (refStrand == 1) {
          byteAlignmentRC.resize(byteAlignment.size());
          MakeReverseComplementByteAlignment(&byteAlignment[0], byteAlignment.size(), &byteAlignmentRC[0]);
          byteAlignment = byteAlignmentRC;
          reverse(insQVChars.begin(), insQVChars.end());
          reverse(delQVChars.begin(), delQVChars.end());
          reverse(subQVChars.begin(), subQVChars.end());
          reverse(mergeQVChars.begin(), mergeQVChars.end());
        }

        string queryAlnString, refAlnString;
        queryAlnString.resize(byteAlignment.size());
        refAlnString.resize(byteAlignment.size());
        ByteAlignmentToRefString(&byteAlignment[0], byteAlignment.size(), &refAlnString[0]);
        ByteAlignmentToQueryString(&byteAlignment[0], byteAlignment.size(), &queryAlnString[0]); 
        
        vector<int> refPositions, queryPositions;
        ComputeQueryPositions(byteAlignment, queryPositions);
        ComputeRefPositions(byteAlignment,   refPositions);

        //
        // Look backwards to see when the alignment should begin.
        //
        int start = refStart;
        int nb = 0;

        
        int readPos;
        int readSubstringStart = 0, readSubstringEnd = 0;
        if (FindPosition(refPositions, centerPos - refStart, readPos)) {

          string alnStr;
          string alnPadding;
          int alnBeginIndex = 0, alnEndIndex = 0;
          SetAlignString(byteAlignment, refPositions, readPos, width, alnStr, alnBeginIndex, alnEndIndex, alnPadding);
          

          MakeMismatchesLowerCase(refStr, alnStr);
          if (miscallView) {
            MakeMatchesDot(refStr, alnStr);
          }
          if (!insertionsOnly) {
            cout << alnStrMargin << alnStr << " " << refStrand << " " << holeNumber << " " << movieId << endl;
          }
          InsertedStringList insertions, insertionQVs;
          if (!hideInsertions or printProfile) {
            StoreInsertedStrings(byteAlignment, refPositions, queryPositions,
                                 insertions, alnBeginIndex, alnEndIndex);

            if (showInsertionQV) {
              string insQVStr;
              CopyFieldString(byteAlignment, insQVChars, alnBeginIndex, alnEndIndex, insQVStr);
              cout << insQVMargin << alnPadding << insQVStr << endl;
            }

            if (showDeletionQV) {
              string delQVStr;
              CopyFieldString(byteAlignment, delQVChars, alnBeginIndex, alnEndIndex, delQVStr);
              cout << delQVMargin << alnPadding << delQVStr << endl;
            }
            if (showSubstitutionQV) {
              string subQVStr;
              CopyFieldString(byteAlignment, subQVChars, alnBeginIndex, alnEndIndex, subQVStr);
              cout << subQVMargin << alnPadding << subQVStr << endl;
            }
            if (showMergeQV) {
              string mergeQVStr;
              CopyFieldString(byteAlignment, mergeQVChars, alnBeginIndex, alnEndIndex, mergeQVStr);
              cout << mergeQVMargin << alnPadding << mergeQVStr << endl;
            }

            if (showInsertionQV) {
              StoreInsertionQVs(insertions, insQVChars);
            }

            if (!hideInsertions) {
              PrintInsertions(insertions, centerPos - width, 
                              refStart, refStr.size(), insMargin, insQVMargin, showInsertionQV, cout);
            }
          }
          if (printProfile) {
            profile.StoreProfile(alnStr, centerPos - width, refStart, insertions);
          }
        }
      }

      if (readsFileName != "") {
        // 
        // Look to see if the reads will be printed.
        //
        MovieNameToArrayIndex::iterator movieNameIt;
        string movieName;
        if (cmpFile.movieInfo.FindMovie(movieId, movieName) == 0) {
          cout << "ERROR in alignment string.  The movie index: " << movieId << " was specified in alignment "
               << alignmentIndex << " but does not exist in the movie info dataset." << endl;
          exit(1);
        }
        
        movieNameIt = movieNameToReaderIndex.find(movieName);
        if (movieNameIt == movieNameToReaderIndex.end()) {
          cout << "Error, an alignment from movie " << movieName << " was specified, but no movie with " << endl
               << " this name exists in " << readsFileName << endl;
          exit(1);
        }
        
        int readerIndex = movieNameIt->second;

        readers[readerIndex].GetReadAt(holeNumber, read);
        read.zmwData.holeNumber = holeNumber;
        
        cout << "alignment index: " << alignmentIndex << " has hole number " << holeNumber << " and movie " << movieId << " with reader " << readerIndex << " zmw id " <<
          read.zmwData.holeNumber << endl;
        if (printedRegionsByHoleNumber[readerIndex].find(read.zmwData.holeNumber) ==
            printedRegionsByHoleNumber[readerIndex].end()) {
          int uniqueHoleNumber = holeNumber + movieId*100000;
          read.zmwData.holeNumber = uniqueHoleNumber;
          basWriter.Write(read);
          read.zmwData.holeNumber = holeNumber;
          int regionLowIndex, regionHighIndex, regionIndex;
          FindRegionIndices(read, &regionTables[readerIndex], regionLowIndex, regionHighIndex);
          for (regionIndex = regionLowIndex; regionIndex < regionHighIndex; regionIndex++) {
            RegionAnnotation region;
            region = regionTables[readerIndex].table[regionIndex];
            region.row[0] = uniqueHoleNumber;
            cout << "Writing row for " << region.row[0] << endl;
            regionWriter.Write(region);
          }
          printedRegionsByHoleNumber[readerIndex].insert(read.zmwData.holeNumber);
        }
      }
    }
    if (printProfile) {
      profile.Print(alnStrMargin.size());
    }
  }
  
  if (readsFileName != "") {
    basWriter.Close();
    regionWriter.Finalize(regionTables[0].columnNames,
                          regionTables[0].regionTypes,
                          regionTables[0].regionDescriptions,
                          regionTables[0].regionSources);
    regionWriter.Close();
  }
  cout << "Done." << endl;
  return 0;
}
