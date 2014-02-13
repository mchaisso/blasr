#ifndef REPMASK_TABLE_H_
#define REPMASK_TABLE_H_


#include <iostream>
#include <string>
#include <vector>

#include "RepmaskRow.h"

#include "../../utils.h"

using namespace std;

class RepmaskTable: public vector<RepmaskRow> {
 public:
  void Read(string &tableFileName) {
    ifstream in;
    CrucialOpen(tableFileName, in, std::ios::in);
    Read(in);
  }

  void Read(istream &in) {
    //
    // The first 3 lines of the table are a header.
    //
    string throwaway;
    getline(in,throwaway);
    getline(in, throwaway);
    getline(in, throwaway);
    while(in) {
      RepmaskRow tmpRow;
      in >> tmpRow;
      this->push_back(tmpRow);
    }
  }
};



#endif
