#ifndef PBONLY_OVERLAP_H_
#define PBONLY_OVERLAP_H_

#include <set>

#include "Types.h"
#include "ReadNode.h"

using namespace std;
class Overlap {
 public:
  Read *targetRead;
  ReadNode  *srcBeginNode, *srcEndNode, *targetBeginNode, *targetEndNode;

  Overlap( Read* targetReadP,
           ReadNode *srcBeginNodeP,  ReadNode *srcEndNodeP,  ReadNode *targetBeginNodeP,  ReadNode *targetEndNodeP) : targetRead(targetReadP),
    srcBeginNode(srcBeginNodeP), 
    srcEndNode(srcEndNodeP), 
    targetBeginNode(targetBeginNodeP), 
    targetEndNode(targetEndNodeP) {
      assert(targetRead != NULL);
      assert(srcBeginNode != NULL);
      assert(srcEndNode   != NULL);
      assert(targetBeginNode != NULL);
      assert(targetEndNode   != NULL);
    }


  int operator<(const Overlap &rhs) const { 
    assert(rhs.srcBeginNode != NULL);
    assert(srcBeginNode != NULL);
    return srcBeginNode->pos < rhs.srcBeginNode->pos;
  }
};


class CompareOverlapPointers {
 public:
  int operator()(const Overlap *lhs, const Overlap *rhs) const { 
    assert(lhs != NULL);
    assert(rhs != NULL);
    return lhs->srcBeginNode->pos < rhs->srcBeginNode->pos;
  }
};

typedef set<Overlap*, CompareOverlapPointers> OverlapSet;


#endif
