#ifndef DATASTRUCTURES_ALIGNMENT_CIGAR_STRING_H_
#define DATASTRUCTURES_ALIGNMENT_CIGAR_STRING_H_

#include <string>
#include <sstream>

using namespace std;
class CigarString : public string {
 public:
  int Vectorize(vector<int> &lengths, vector<char> &operations) {
    stringstream strm;
    strm.str(*this);
    while (strm) {
      int l;
      char o;
      if ((strm >> l >> o)) {
        lengths.push_back(l);
        operations.push_back(o);
      }
    }
  }
};


#endif
