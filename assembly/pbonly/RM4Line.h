#ifndef RM4Line_H_
#define RM4Line_H_

#include <istream>

class RM4Line {
 public:
  //m120308_205818_42142_cunknown0_s1_p0/94/0_419 m120308_205818_42142_cunknown0_s1_p0/94/0_419 -2095 100 0 0 419 491 0 0 419 419 0 13726 -377.391 -2095 1

  string qName, tName;
  int    alignScore;
  float  pctIdentity;
  unsigned int qStrand, tStrand;
  unsigned int qBegin, qEnd, qLength;
  unsigned int tBegin, tEnd, tLength;
  
};

std::istream& operator >>(std::istream &is, RM4Line &rm4Line) {
  if (! (is >> rm4Line.qName) ) return is;
  if (! (is >> rm4Line.tName) ) return is;
  if (! (is >> rm4Line.alignScore) ) return is;
  if (! (is >> rm4Line.pctIdentity) ) return is;
  if (! (is >> rm4Line.qStrand ) ) return is;
  if (! (is >> rm4Line.qBegin) ) return is;
  if (! (is >> rm4Line.qEnd) ) return is;
  if (! (is >> rm4Line.qLength) ) return is;
  if (! (is >> rm4Line.tStrand) ) return is;
  if (! (is >> rm4Line.tBegin) ) return is;
  if (! (is >> rm4Line.tEnd) ) return is;
  if (! (is >> rm4Line.tLength) ) return is;
  string remainder;
  getline(is, remainder);
  return is;
}

  
  


#endif
