#ifndef SAM_KEYWORD_VALUE_PAIR_H_
#define SAM_KEYWORD_VALUE_PAIR_H_

#include <string>
#include <iostream>
#include <vector>

#include <stdlib.h>
#include "utils/StringUtils.h"

using namespace std;

class SAMKeywordValuePair {
 public:
  string key;
  string value;
};

bool SplitSAMKeyValuePair(string &kvPair, string &key, string &value, string token=":") {
  int sepIndex = kvPair.find_first_of(token);
  if (sepIndex == kvPair.npos) {
    return false;
  }
  else {
    key = kvPair.substr(0, sepIndex);
    value = kvPair.substr(sepIndex+1);
    return true;
  }
}

bool SplitSAMTypedKeyValuePair(string kvPair, string &key, string &kvType, string &value, char token=':') {
  vector<string> strValues;
  ParseSeparatedList(kvPair, strValues, token, 3);
  if (strValues.size() != 3) {
    return false;
  }
  else {
    key = strValues[0];
    kvType = strValues[1];
    value  = strValues[2];
    return true;
  }
}

template<typename T>
void StoreValue(string &valueStr, T&value) {
  stringstream strm(valueStr);
  if (!(strm >> value)) {
    cout <<"Error parsing " << valueStr << endl;
    exit(1);
  }
}

void KeywordValueStringsToPairs(vector<string> &kvStrings, vector<SAMKeywordValuePair> &kvPairs, string token=":") {

  kvPairs.resize(kvStrings.size());

  if (kvStrings.size() == 0) {
    return;
  }

  int i;
  for (i = 0; i < kvStrings.size(); i++ ) {
    SplitSAMKeyValuePair(kvStrings[i], kvPairs[i].key, kvPairs[i].value, token);
  }
}

class TypedKeywordValuePair {
 public:
  static bool Separate(string &kvPair, string &kvKey, string &kvType, string &kvValue) {
    if (SplitSAMTypedKeyValuePair(kvPair, kvKey, kvType, kvValue) == false) {
      return false;
    }
    return true;
  }
};


template<typename T_Value>
class KeywordValuePair {
 public:

  static bool Parse(string &kvPair, const char *key, T_Value &result) {
    string kvKey, kvValue;
    SplitSAMKeyValuePair(kvPair, kvKey, kvValue);
    if (kvKey != key) {
      return false;
    }
    else {
      stringstream strm(kvValue);
      if ( ! (kvValue >> result) ) {
        return false;
      }
      else {
        return true;
      }
    }
  }
  
  static bool Store(string &valueStr, T_Value &value) {
    return (stringstream(valueStr) >> value);
  }


};

#endif 
