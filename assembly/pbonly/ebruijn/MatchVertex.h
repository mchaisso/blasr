#ifndef EBRUIJN_MATCH_VERTEX_H_
#define EBRUIJN_MATCH_VERTEX_H_

#include <vector>
#include <set>
using namespace std;
class MatchVertex {
 public:
  int pos;
  int seqIndex, seqPos;
  char strand;
  set<MatchVertex*> next;
  set<MatchVertex*> prev;
  
  MatchVertex() {
    pos = 0;
  }
  
  MatchVertex(int p) : pos(p) {}

  bool operator<(const MatchVertex &rhs) const {
    return pos < rhs.pos;
  }
};

typedef vector<MatchVertex> MatchVertexList;

#endif
