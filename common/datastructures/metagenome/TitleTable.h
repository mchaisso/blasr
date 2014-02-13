#ifndef DATASTRUCTURES_METAGENOMICS_TITLE_TABLE_H_
#define DATASTRUCTURES_METAGENOMICS_TITLE_TABLE_H_

#include <string.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "../../utils.h"
using namespace std;
class TitleTable {
 public:
	char **table;
	int tableLength;
  TitleTable() {
    table = NULL;
    tableLength = 0;
  }
	void Copy(char **src, int nSrc) {
		table = new char*[nSrc];
		tableLength = nSrc;
		int i;
		for (i = 0; i < nSrc; i++ ){
			int lenStrI = strlen(src[i]);
			table[i] = new char[lenStrI+1];
			memcpy(table[i], src[i], lenStrI);
			table[i][lenStrI] = '\0';
		}
	}

  void Write(string &name) {
    ofstream out;
    CrucialOpen(name, out, std::ios::out);
    Write(out);
  }

	void Write(ofstream &out) {
		int i;
		for (i = 0; i < tableLength;i++) {
			out << table[i] << endl;
		}
	}
  void Read(string &inFileName) {
    ifstream in;
    CrucialOpen(inFileName, in, std::ios::in);
    Read(in);
  }
  void CopyFromVector(vector<string> &titles) {
    tableLength = titles.size();
    table = new char*[tableLength];
    int i;
    for (i = 0; i < tableLength; i++) {
      table[i] = new char[titles[i].size() + 1];
      memcpy(table[i], titles[i].c_str(), titles[i].size());
      table[i][titles[i].size()] = '\0';
    }
  }

  void Read(ifstream &in) {
    vector<string> titles;
    while(in) {
      string title;
      if (in >> title) {
        titles.push_back(title);
      }
    }
    if (titles.size() > 0) {
      CopyFromVector(titles);
    }
    else {
      tableLength = 0;
      table = NULL;
    }
  }

  void Free() {
    int i;
    for (i = 0; i < tableLength; i++) {
      delete[] table[i];
    }
    delete[] table;
  }

  bool Lookup(string title, int &index) {
    int i;
    for (i = 0; i < tableLength; i++) {
      if (table[i] == title) { 
        index = i; 
        return true;
      }
    }
    return false;
  }
    
	static void ResetTableToIntegers(char **table, int *tableLengths, int nTable) {
		int i;
		for (i = 0; i < nTable; i++ ) {
			delete[] table[i];
			 stringstream namestrm;
			namestrm << i;
			string name;
			name = namestrm.str();
			table[i] = new char[name.size()+1];
			memcpy( table[i], name.c_str(), name.size());
			table[i][name.size()] = '\0';
			tableLengths[i] = (int) name.size() + 1;
		}
  }
	
};


#endif
