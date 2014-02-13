#ifndef CDNA_GTFDB_H_
#define CDNA_GTFDB_H_

#include <map>
#include <set>
#include <string>
#include "GencodeGFFReader.h"
#include "GencodeGFFEntry.h"
#include "utils.h"

using namespace std;


class GTFDBEntry {
 public:
  int start, end;
  string type;
  string line;
  string gene;
  string transcriptID;
  int operator<(const GTFDBEntry &rhs) const {
    return start < rhs.start;
  }
 GTFDBEntry(int s, int e, string t, string g, string tid, string l) : 
  start(s), 
    end(e), 
    type(t), 
    gene(g),
    transcriptID(tid),
    line(l) {}
};

class GTFChromosome {
 public:
  multiset<GTFDBEntry> entries;
};

class GTFDB : public map<string, GTFChromosome*> {
 public:
  void InsertEntry(string line) {
    GencodeGFFEntry fullEntry;
    stringstream strm(line);
    ReadGencodeGFFLine(strm, fullEntry);
    if (find(fullEntry.chr) == end()) {
      GTFChromosome *newChr = new GTFChromosome;
      insert(pair<string, GTFChromosome*>(fullEntry.chr, newChr));
    }
    (*this)[fullEntry.chr]->entries.insert(GTFDBEntry(fullEntry.start, fullEntry.end, fullEntry.genomicLocusType, fullEntry.geneName, fullEntry.transcriptId, line));
  }

  GTFChromosome *GetChromosome(string chr) {
    if (find(chr) == end()) {
      return NULL;
    }
    else {
      return (*this)[chr];
    }
  };
  
  void ReadGTFFile(string gtfFileName) {
    ifstream in;
    CrucialOpen(gtfFileName, in, std::ios::in);
    string line;
    while(getline(in, line)) {
      if (line.size() == 0 or line[0] == '#') {
        continue;
      }
      InsertEntry(line);
    }
  }

  int FindEntries(string chr, int start, int end,
                  multiset<GTFDBEntry>::iterator & entryBegin, 
                  multiset<GTFDBEntry>::iterator & entryEnd) {
    GTFChromosome* chrPtr = GetChromosome(chr);

    if (chrPtr == NULL) {
      return 0;
    }
    multiset<GTFDBEntry>::iterator entryIt, prevIt;
    GTFDBEntry query(start, 0, "", "", "", "");
    entryBegin = chrPtr->entries.lower_bound(query);
    entryEnd = chrPtr->entries.end();
    if (entryBegin == chrPtr->entries.end()) {
      return 0;
    }
    
    entryIt = entryBegin;
    int numEntries = 0;
    prevIt = entryIt;
    if (prevIt != chrPtr->entries.begin()) {
      --prevIt;
      if (prevIt != chrPtr->entries.end() and 
          prevIt != chrPtr->entries.begin() and 
          ((prevIt->start <= start and 
            prevIt->end >= start) or
           (prevIt->start <= start and
            prevIt->end >= start) or
           (prevIt->start <= start and
            prevIt->end >= end) or
           (prevIt->start >= start and
            prevIt->end <= end))) {
        --entryBegin;
        --prevIt;
      }
      entryIt = entryBegin;
    }
    while(entryIt != chrPtr->entries.end() and 
          ((entryIt->start <= start and 
            entryIt->end >= start) or 
           (entryIt->start <= end and 
            entryIt->end >= end) or
           (entryIt->start <= start and
            entryIt->end >= end) or
           (entryIt->start >= start and
            entryIt->end <= end))) {
      numEntries++;
      ++entryIt;
    }
    entryEnd = entryIt;
    return numEntries;
  }

  int FindEntries(string chr, int start, int end, vector<GTFDBEntry> &dbEntries) {
    int numEntries;
    multiset<GTFDBEntry>::iterator entryBegin;
    multiset<GTFDBEntry>::iterator entryEnd;
    numEntries = FindEntries(chr, start, end, entryBegin, entryEnd);
    if (numEntries == 0) {
      return 0;
    }
    for (; entryBegin != entryEnd; ++entryBegin) {
      dbEntries.push_back(*entryBegin);
    }
    return dbEntries.size();
  }
  int FindEntries(string chr, int start, int end, vector<string> &entries) {
    int numEntries;
    multiset<GTFDBEntry>::iterator entryIt;
    multiset<GTFDBEntry>::iterator entryEnd;
    numEntries = FindEntries(chr, start, end, entryIt, entryEnd);
    if (numEntries == 0) {
      return 0;
    }
    for (; entryIt != entryEnd; ++entryIt) {
      entries.push_back(entryIt->line);
    }
    return entries.size();
  }
};

#endif
