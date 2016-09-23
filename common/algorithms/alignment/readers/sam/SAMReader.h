#ifndef SAM_READER_H_
#define SAM_READER_H_

#include "datastructures/alignmentset/SAMAlignment.h"
#include "datastructures/alignmentset/SAMKeywordValuePair.h"
#include "datastructures/alignmentset/ReadGroup.h"
#include "datastructures/alignmentset/ReferenceSequence.h"
#include "datastructures/alignmentset/AlignmentSet.h"
#include "utils/StringUtils.h"
#include "utils.h"

#include <sstream>

template<typename T_ReferenceSequence=SAMReferenceSequence, typename T_ReadGroup=SAMReadGroup, typename T_SAMAlignment=SAMAlignment>
class SAMReader {
  public:
  int lineNumber;
  ifstream samFile;
  bool Initialize(string samFileName) {
    CrucialOpen(samFileName, samFile, std::ios::in);
		return true;
  }

  void Close() {
    samFile.close();
  }

  enum LineType {Blank, HSHeader, HSSequence, HSReadGroup, HSProgram, HSComment, Alignment, Error};

  int GetLine(istream &in, string &line) {
    return (getline(in, line));
  }

  bool LineTypeIsHeader(LineType lineType) {
    return (lineType == HSHeader or
            lineType == HSSequence or
            lineType == HSReadGroup or 
            lineType == HSProgram or 
            lineType == HSComment);
  }

  bool PeekLineIsHeader(istream &in) {
    if (in and in.peek() == '@') {
      return true;
    }
    else {
      return false;
    }
  }

  LineType GetLineType(string &line) {
    if (line.length() == 0) {
      return Blank;
    }
    else if (line[0] == '@') {
      stringstream strm(line);
      string tag;
      strm >> tag;
      if (tag == "@HD") { return HSHeader; }
      else if (tag == "@SQ") { return HSSequence; }
      else if (tag == "@RG") { return HSReadGroup; }
      else if (tag == "@PG") { return HSProgram; }
      else if (tag == "@CO") { return HSComment; }
      else { return Error; }
    }
    else {
      return Alignment;
    }
  }
 
  void StoreKVPairs(string line, vector<SAMKeywordValuePair> &kvPairs, string token="\t") {
    //
    // Split on tab delineated line.
    //
    vector<string> kvPairStrings;
    Tokenize(line, token, kvPairStrings);
    KeywordValueStringsToPairs(kvPairStrings, kvPairs);
  }
 
  int StoreHeader(vector<SAMKeywordValuePair> &kvPairs, AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments) {
    alignments.header.StoreValues(kvPairs, lineNumber);
  }

  void StoreReferenceSequence(vector<SAMKeywordValuePair> &kvPairs,
                              AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments) {
    alignments.references.push_back(T_ReferenceSequence());
    int lastRefIndex = alignments.references.size() - 1;
    alignments.references[lastRefIndex].StoreValues(kvPairs, lineNumber);
		alignments.refNameToIndex[alignments.references[lastRefIndex].sequenceName] = lastRefIndex;
  }

  void StoreReadGroup(vector<SAMKeywordValuePair> &kvPairs,
                      AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments) {
    alignments.readGroups.push_back(T_ReadGroup());
    int lastReadGroupIndex = alignments.readGroups.size() - 1;
    alignments.readGroups[lastReadGroupIndex].StoreValues(kvPairs, lineNumber);
  }

  void StoreAlignment(string & line,
                      AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments) {
    alignments.alignments.push_back(T_SAMAlignment());   
    int lastAlignmentIndex = alignments.alignments.size() - 1;
    alignments.alignments[lastAlignmentIndex].StoreValues(line, lineNumber);
  }

  void StoreProgram(vector<SAMKeywordValuePair> &kvPairs,
                    AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments ) {
    //
    // Do this later
    //
  }

  void Read(string samFileName, AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments) {
    Initialize(samFileName);
    Read(alignments);
  }

  vector<string> ReadHeader(AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments) {
    vector<string> allHeaders;
    string line;
    LineType lineType;
    lineNumber = 0;
    while (samFile and PeekLineIsHeader(samFile)) {
      getline(samFile, line);
      lineType = GetLineType(line);
      if (LineTypeIsHeader(lineType)) {
        allHeaders.push_back(line);
        stringstream strm(line);
        string tag;
        strm >> tag;
        string remainder;
        getline(strm, remainder);
        vector<SAMKeywordValuePair> kvPairs;
        StoreKVPairs(remainder, kvPairs);
        if (lineType == HSHeader) {
          StoreHeader(kvPairs, alignments);
        }
        else if (lineType == HSSequence) {
          StoreReferenceSequence(kvPairs, alignments);
        }
        else if (lineType == HSReadGroup) {
          StoreReadGroup(kvPairs, alignments);
        }
        else if (lineType == HSProgram) {
          StoreProgram(kvPairs, alignments);
        }
        else if (lineType == HSComment) {
          // do nothing with comments for now.
        }
      }
      ++lineNumber;
    }
    return allHeaders;
  }

  void Read(AlignmentSet<T_ReferenceSequence, T_ReadGroup, T_SAMAlignment> &alignments) {
    string line;
    LineType lineType;
    lineNumber = 0;
    ReadHeader(alignments);
    while (getline(samFile, line)) {
      lineType = GetLineType(line);
      if (LineTypeIsHeader(lineType)) {
        cout << "ERROR! Header line found outside of the header at " << lineNumber << endl;
        exit(1);
      }
      else if (lineType == Alignment) {
        StoreAlignment(line, alignments);
      }
      else {
        cout << "Error, line type unknown at " << lineNumber << endl;
        cout << line << endl;
        exit(1);
      }
      ++lineNumber;
    }
  }

  bool GetNextAlignment( SAMAlignment& alignment, bool allowUnaligned=false) {
    if (samFile) {
      string line;
      if (getline(samFile, line)) {
        alignment.StoreValues(line, lineNumber, allowUnaligned);
        ++lineNumber;
        return true;
      }
      else {
        return false;
      }
    }
    else {
      return false;
    }
  }
};

#endif
