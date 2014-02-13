#include "cdna/GTFDB.h"
#include "utils.h"
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

class Alignment {
public:
  string type;
  string chr;
  string gene;
  string transcriptId;
  int start, end, score;
  float identity;
  int strand;
  int tLength;
};

void ReadGFF3Line(istream &in, Alignment &aln) {
  string tmp;
  in >> aln.chr >> tmp >> aln.type >> aln.start >> aln.end >> tmp >> aln.strand;
}

void ReadM4Line(istream &in, Alignment &aln) {
  string name, tmp;
  string qs, qst, qe;
  in >> name >> aln.chr >> aln.score >> aln.identity >> qs >> qst >> qe >> tmp >> aln.strand >> aln.start >> aln.end >> aln.tLength;
}

int main(int argc, char* argv[]) {

  string gtfFileName, inputFileName;
  if (argc < 3) {
    cout << "usage: matchGTFExons align_file gtfFile [-format (gff3)]" << endl;
    exit(1);
  }

  inputFileName = argv[1];
  gtfFileName = argv[2];
  int argi = 3;
  string format = "gff3";
  while (argi < argc) {
    if (strcmp(argv[argi], "-format") == 0) {
      ++argi;
      format = argv[argi];
    }
    else {
      cout << "bad option " << argv[argi] << endl;
      exit(0);
    }
    ++argi;
  }
  ifstream inFile;
  CrucialOpen(inputFileName, inFile, std::ios::in);

  GTFDB gtfdb;
  gtfdb.ReadGTFFile(gtfFileName);
  
  while(inFile) {
    string line;
    getline(inFile, line);
    stringstream strm(line);
    Alignment aln;
    if (format == "gff3") {
      if (line[0] == '#') {
        continue;
      }
      ReadGFF3Line(strm, aln);
      if (aln.type == "gene" or aln.type == "mRNA") {
        continue;
      }
    }
    else if (format == "m4") {
      ReadM4Line(strm, aln);
      if (aln.strand == 1) {
        swap(aln.start, aln.end);
        aln.start = aln.tLength - aln.start;
        aln.end   = aln.tLength - aln.end + 1;
      }
    }
    vector<GTFDBEntry> entryLines;
    gtfdb.FindEntries(aln.chr, aln.start, aln.end, entryLines);
    int e;
    cout << "searching for " << aln.chr << " " << aln.start << " " << aln.end << endl;
    bool foundMatch = false;
    int maxEntryLine = 0;
    int maxEntryOverlap = 0;
    if (entryLines.size() > 0) {
      for(e = 0; e < entryLines.size(); e++) {
        float overlap;
        int overlapStart = max(aln.start, entryLines[e].start);
        int overlapEnd   = min(aln.end, entryLines[e].end);
        int intvStart    = min(aln.start, entryLines[e].start);
        int intvEnd      = max(aln.end, entryLines[e].end);

        if (overlapEnd < overlapStart) {
          overlap = 0;
        }
        else {
          if (aln.end <= aln.start or intvEnd <= intvStart) {
            overlap = 0;
          }
          else {
            overlap = 100 * (overlapEnd - overlapStart) / (intvEnd - intvStart*1.0);
          }
        }
        if (overlap > maxEntryOverlap) {
          maxEntryLine = e;
          maxEntryOverlap = overlap;
        }
      }
    }
    if (entryLines.size() == 0 or maxEntryOverlap  == 0) {
      cout << line << endl;
      cout << "   *** NO MATCHES ***" << endl;
    }
    else {
      e = maxEntryLine;
      cout << " " << maxEntryOverlap << " match " << aln.chr << " " << entryLines[e].start << " " << entryLines[e].end << " " << entryLines[e].type << " " << entryLines[e].gene << " " << entryLines[e].transcriptID << endl;
    }
  }
}
