#ifndef READ_NODE_H_
#define READ_NODE_H_
#include <assert.h>
#include <set>
#include "Types.h"
#include <vector>

using namespace std;

class ReadNode;
class Read;
class  ReadNodeSuperset : public set<ReadNode*> {
 public:
  static UInt counter;
  static UInt countDeleted;
  UInt id;
  bool traversed;
  void ClearFlags() {
    traversed = false;
  }
  ReadNodeSuperset() {
    // Make a counter so we can keep track of things.
    id = ++counter;
    ClearFlags();
  }

  bool ContainsNode(ReadNode* readNode) {
    return (find(readNode) != end());
  }


};

UInt ReadNodeSuperset::counter = 0;
UInt ReadNodeSuperset::countDeleted = 0;

class OverlapReadNodes {
 public:
  ReadNode *thisEnd, *overlapBegin, *overlapEnd;
 OverlapReadNodes(ReadNode *te, ReadNode *ob, ReadNode *oe) : thisEnd(te), overlapBegin(ob), overlapEnd(oe) {}
};

class ReadNode {
 public:
  UInt pos;
  UInt coverage;
  static UInt counter;
  UInt id;
  ReadNodeSuperset *super;
  Read *read;
  typedef vector<OverlapReadNodes> Overlaps;
  Overlaps overlaps;
   
  void AddOverlap(ReadNode *endNode, ReadNode *overlapBeginNode, ReadNode *overlapEndNode) {
    // Make sure all are assigned
    assert(endNode != NULL);
    assert(overlapBeginNode != NULL);
    assert(overlapEndNode != NULL);
    overlaps.push_back(OverlapReadNodes(endNode, overlapBeginNode, overlapEndNode));
  }

  void AddNodeToSuperset(ReadNode *node) {
    if (super == NULL) {
      // 
      // If the superset does not exist yet, make it, and add self to
      // it. 
      //
      super = new ReadNodeSuperset;
      super->insert(this);
    }
    super->insert(node);
    node->super = super;
  }

  void AbsorbSuperset(ReadNode* node) {
    if (super == NULL and node->super == NULL) {
      //
      // No supersets exist, make one out of these two nodes.
      //
      AddNodeToSuperset(node);
    }
    else if (super == NULL) {
      node->AddNodeToSuperset(this);
      super = node->super;
    }
    else if (node->super == NULL) {
      AddNodeToSuperset(node);
    }
    else {
      //
      // Both this super and the other contain nodes, merge them, and
      // delete the other superset.
      //
      ReadNodeSuperset::iterator ssIt, ssEnd;
      ReadNodeSuperset *supersetToBeMerged = node->super;
      for (ssIt = node->super->begin(); ssIt != node->super->end(); ++ssIt) {
        (*ssIt)->super = super;
        super->insert(*ssIt);
      }
      delete supersetToBeMerged;
    }
  }


  int operator<(const UInt rhsPos) {
    return pos < rhsPos;
  }
  
  int operator<(const ReadNode &rhs) {
    return pos < rhs.pos;
  }

  ReadNode() {
    super = NULL;
    read  = NULL;
    coverage = 0;
    id = ++counter;
  }

 ReadNode(UInt posP) : pos(posP) {
    super = NULL;
    read  = NULL;
    coverage = 0;
  }

  ReadNode(UInt posP, Read* readPtr)  {
   pos      = posP;
   super    = NULL;
   read     = readPtr;
   coverage = 0;
   id = ++counter;
  }
  
  int operator<(const ReadNode &rhs) const {
    return pos < rhs.pos;
  }

  bool SupersetContainsRead(Read *read) {
    ReadNodeSuperset::iterator it, end;
    if (super == NULL) {
      return false;
    }
    else {
      end = super->end();
      it  = super->begin();
      while (it != end) {
        if ((*it)->read == read) {
          return true;
        }
      }
      return false;
    }
  }

  void RemoveFromSuperset() {
    if (super != NULL) {
      ReadNodeSuperset::iterator it;
      it = super->find(this);
      assert(it != super->end());
      super->erase(it);
      //
      // Removed the last element of the superset, erase this so that
      // dangling pointers are not laying around.
      //
      if (super->size() == 0) {
        ReadNodeSuperset::countDeleted++;
        delete super;
        super = NULL;
      }
    }
  }
};


UInt ReadNode::counter = 0;

class CompareReadNodePointers {
 public:
  int operator()(const ReadNode *lhs, const ReadNode *rhs) {
    return lhs->pos < rhs->pos;
  }
};

typedef set<ReadNode*, CompareReadNodePointers> ReadNodeSet;


#endif
