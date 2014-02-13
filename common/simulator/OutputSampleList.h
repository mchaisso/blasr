#ifndef SIMULATOR_OUTPUT_SAMPLE_LIST_H_
#define SIMULATOR_OUTPUT_SAMPLE_LIST_H_

#include "OutputSample.h"
#include <vector>

using namespace std;
class OutputSampleList: public vector<OutputSample> {
 public:

  void Write(ofstream &out) {
    int nElem = this->size();
    out.write((char*)&nElem, sizeof(int));
    if (nElem > 0) {
      int i;
      for (i = 0; i < nElem; i++) {
        (*this)[i].Write(out);
      }
    }
  }

  void Read(ifstream &in) {
    int nElem;
    in.read((char*) &nElem, sizeof(int));
    if (nElem > 0) {
      this->resize(nElem);
      int i;
      for (i = 0; i < nElem; i++) {
        (*this)[i].Read(in);
      }
    }
  }
};

    

#endif
