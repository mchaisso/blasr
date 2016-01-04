#ifndef ALIGNMENT_SET_READ_GROUP_H_
#define ALIGNMENT_SET_READ_GROUP_H_

#include <string>

/*
 * Minimal components of read group.  Define only required items. Use
 * if memory is scarce.  
 */
#include "SAMKeywordValuePair.h"

class SAMReadGroup {
 public:
  string id;
  void StoreValues(vector<SAMKeywordValuePair> &kvPairs, int lineNumber = 0) {
    int i;
    bool idIsStored = false;
    for (i = 0; i < kvPairs.size(); i++ ){
      if (kvPairs[i].key == "ID") {
        id = kvPairs[i].value;
        idIsStored = true;
      }
    }
    if (idIsStored == false) {
      cout << "ReadGroup missing id at " << lineNumber << endl;
      exit(1);
    }
  }
};

/*
 * Full read group. Use when all data are required.
 */
class SAMFullReadGroup : public SAMReadGroup {
 public:
  string centerName;
  string description;
  string date;
  string flowOrder;
  string keySequence;
  string library;
  string processingProgram;
  int    insertSize;
  string platform;
  string platformUnit;
  string sample;
	string bindingKit;
	string sequencingKit;
	string basecallerVersion;
	string movieName;
	string changelistId;
	SAMFullReadGroup() {
		movieName=basecallerVersion=sequencingKit="";
	}

  void StoreValues(vector<SAMKeywordValuePair> &kvPairs, int lineNumber = 0) {
    SAMReadGroup::StoreValues(kvPairs, lineNumber);
    string kwPair;
    string key, valueStr;
    int i;

    for (i = 0; i < kvPairs.size(); i++) {
      if (kvPairs[i].key == "CN") {
        centerName = kvPairs[i].value;
      }
      else if (kvPairs[i].key == "DS") {
        description = kvPairs[i].value;
				vector<SAMKeywordValuePair> descKvPairs;
				vector<string> kvPairStrings;
				Tokenize(description, ";", kvPairStrings);
				KeywordValueStringsToPairs(kvPairStrings, descKvPairs, "=");
				//				StoreValues(descKvPairs, lineNumber);
				kvPairs.insert(kvPairs.end(), descKvPairs.begin(), descKvPairs.end());
      }
      else if (kvPairs[i].key == "DT") {
        date = kvPairs[i].value;
      }
      else if (kvPairs[i].key == "FO") {
        flowOrder = kvPairs[i].value;
      }
      else if (kvPairs[i].key == "LB") {
        library = kvPairs[i].value;
      }
      else if (kvPairs[i].key == "PG") {
        processingProgram = kvPairs[i].value;
      }
      else if (kvPairs[i].key == "PI") {
        StoreValue(kvPairs[i].value, insertSize);
      }
      else if (kvPairs[i].key == "SM") {
        sample = kvPairs[i].value;
      }
			else if (kvPairs[i].key == "CHANGELISTID") {
				changelistId = kvPairs[i].value;
			}
			else if (kvPairs[i].key == "BINDINGKIT") {
				bindingKit = kvPairs[i].value;
			}
			else if (kvPairs[i].key == "SEQUENCINGKIT") {
				sequencingKit = kvPairs[i].value;
			}
			else if (kvPairs[i].key == "BASECALLERVERSION") {
				basecallerVersion = kvPairs[i].value;
			}
			else if (kvPairs[i].key == "PU") {
				movieName = kvPairs[i].value;
			}
    }
  }
};


#endif
