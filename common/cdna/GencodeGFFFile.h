#ifndef CDNA_GENCODE_GFF_FILE_H_
#define CDNA_GENCODE_GFF_FILE_H_

#include "GencodeGFFEntry.h"
#include "GencodeGFFReader.h"
#include "utils/FileUtils.h"
#include "utils.h"

#include <vector>
using namespace std;

class GencodeGFFFile {
 public:
  vector<GencodeGFFEntry> entries;
  /*
   * Read all gff entries into a file. 
   * TODO: limit to subsets of gene_type .
   */
  void ReadAll(string gffFileName) {

    ifstream in;
    CrucialOpen(gffFileName, in, std::ios::in);

    while(in) {
      string line;
      getline(in, line);
      if (line.size() == 0 or line[0] == '#') {
        continue;
      }
      else {
        stringstream strm(line);
        GencodeGFFEntry entry;
        ReadGencodeGFFLine(strm, entry);
        entries.push_back(entry);
      }
    }
  }
};



#endif
