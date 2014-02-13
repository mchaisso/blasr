#ifndef PBIHDFUTILS_INSERTED_STRING_H_
#define PBIHDFUTILS_INSERTED_STRING_H_
class InsertedString {
public:
  string insSeq;
  string insQVSeq;
  int    pos;
  int    queryPos;
  int    alnPos;
  bool   printed;
  InsertedString() {
    pos = -1;
    printed = false;
  }
  InsertedString(string &s, int p, int q, int a) {
    insSeq = s;
    pos = p;
    queryPos = q;
    alnPos   = a;
    printed = false;
  }
};
typedef vector<InsertedString> InsertedStringList;


#endif
