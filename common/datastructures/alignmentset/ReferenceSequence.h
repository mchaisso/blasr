#ifndef ALIGNMENT_SET_REFERENCE_SEQUENCE_H_
#define ALIGNMENT_SET_REFERENCE_SEQUENCE_H_

#include "SAMKeywordValuePair.h"
#include <string>

class SAMReferenceSequence {
 public:
  string sequenceName;
  unsigned int length;
  string GetSequenceName() {
    return sequenceName;
  }
  //
  // By definining accessor functions here but no data, we can
  // economize on the amount of space used for each element.   This is
  // no big deal for references, but for pairwise alignments, it is
  // big.
  //
  string GetMD5() {
    return "";
  }
  unsigned int GetLength() {
    return length;
  }
  string GetGenomeAssembly() {
    return "";
  }
  string GetSpecies() {
    return "";
  }
  string GetURI() {
    return "";
  }

  enum SAMReferenceSequenceRequiredFields { SQ_SN, SQ_LN };
  static const char* SAMReferenceSequenceFieldNames[];
  void StoreValues(vector<SAMKeywordValuePair> &kvPairs, int lineNumber=0) {
    int i;
    vector<bool> usedFields;
    usedFields.resize(SQ_LN);
    fill(usedFields.begin(), usedFields.end(), false);
    for (i = 0; i < kvPairs.size(); i++) {
      if (kvPairs[i].key == "SN") {
        sequenceName = kvPairs[i].value;
        usedFields[SQ_SN] = true;
      }
      else if (kvPairs[i].key == "LN") {
        StoreValue(kvPairs[i].value, length);
        usedFields[SQ_SN] = true;
      }
    }
    for (i = 0; i < usedFields.size(); i++) {
      if (usedFields[i] == false) {
        cout << "SQ specifier missing " << SAMReferenceSequenceFieldNames[i] << endl;
      }
    }
  }
};
const char* SAMReferenceSequence::SAMReferenceSequenceFieldNames[] = {"SN", "LN"};

class SAMFullReferenceSequence : public SAMReferenceSequence {
 public:
  string md5;
  string species;
  string uri;
  string genomeAssembly;
  string GetMD5() {
    return md5;
  }
  string GetSpecies() {
    return species;
  }
  string GetURI() {
    return uri;
  }
  string GetGenomeAssembly() {
    return genomeAssembly;
  }

  enum FullReferenceSequenceRequiredFields { SQ_AS, SQ_M5, SQ_SP, SQ_UR };
  static const char* SAMFullReferenceSequenceFieldNames[];
  void StoreValues(vector<SAMKeywordValuePair> &kvPairs, int lineNumber=0) {
    SAMReferenceSequence::StoreValues(kvPairs, lineNumber);
    int i;
    for (i = 0; i < kvPairs.size(); i++ ){
      if (kvPairs[i].key == "AS") {
        genomeAssembly = kvPairs[i].value;
      }
      else if (kvPairs[i].key == "M5") {
        md5 = kvPairs[i].value;
      }
      else if (kvPairs[i].key == "SP") {
        species = kvPairs[i].value;
      }
      else if (kvPairs[i].key == "UR") {
        uri = kvPairs[i].value;
      }
    }
  }
};

const char* SAMFullReferenceSequence::SAMFullReferenceSequenceFieldNames[] = {"AS", "M5", "SP", "UR"};

#endif
