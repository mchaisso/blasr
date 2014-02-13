#ifndef READ_H_
#define READ_H_

#include "SameMolecule.h"
#include "ReadNode.h"

#include <assert.h>
#include "Overlap.h"
#include <algorithm>
// Each read is a subread from the same molecule.  Define that here.
class SameMolecule;

class Read {
 public:
  ReadNodeSet orderedReadNodes;
  ReadNodeSet orderedSameMoleculeNodes;
  SameMolecule *myMolecule;
  UInt subreadBegin, subreadEnd;
  UInt readLength;
  UInt fullReadLength;
  // Overlaps ordered by read
  set<Read*> overlapsWith;
  // Overlaps ordered by position
  OverlapSet overlaps;
  static UInt counter;
  UInt id;
  Read() {
    subreadBegin = subreadEnd = 0;
    myMolecule = NULL;
    id = counter;
    counter++;
  }

  bool OverlapsWith(Read* targetPtr) {
    return overlapsWith.find(targetPtr) != overlapsWith.end();
  }

  void InitializeEndNodes() {
    if (readLength == 0) {
      //
      // This operation is undefined for length 0 reads.
      // Bail.
      //
      return;
    }
    ReadNode *firstNode = new ReadNode(0, this);
    ReadNode *lastNode  = new ReadNode(readLength, this);
    orderedReadNodes.insert(firstNode);
    firstNode->coverage = 1;
    orderedReadNodes.insert(lastNode);
    lastNode->coverage = 1;
  }

  void AddOverlap(Read* destRead, ReadNode *beginNode, ReadNode *endNode, ReadNode *targetBegin, ReadNode *targetEnd) {
    overlaps.insert(new Overlap(destRead, beginNode, endNode, targetBegin, targetEnd));
    overlapsWith.insert(destRead);
  }


  void RemoveInterval(ReadNode *beginNode, ReadNode *endNode) {
    assert(beginNode != NULL);
    assert(endNode!= NULL);

    //
    // First decrement the coverage along the path corresponding to
    // this read interval.
    //
    ReadNodeSet::iterator it = orderedReadNodes.find(beginNode);
    ReadNodeSet::iterator end = orderedReadNodes.find(endNode);

    assert(it != orderedReadNodes.end());
    assert(end != orderedReadNodes.end());
    end++;
    for (; it != end; it++) {
      (*it)->coverage--;
    }

    //
    // Next, remove these intervals from any superset.
    //
    beginNode->RemoveFromSuperset(); 
    endNode->RemoveFromSuperset();
  }

  bool RemoveOverlapsWith(Read *targetRead) {
    set<Read*>::iterator ovpIt;
    ovpIt = overlapsWith.find(targetRead);
    if (ovpIt == overlapsWith.end()) {
      return false;
    }
    else {
      overlapsWith.erase(ovpIt);
      return true;
    }
  }
  
  void RemoveOverlaps() {
    //
    // Remove all intervals from this read, and decrement any overlaps
    // it has elsewhere.
    //
    OverlapSet::iterator ovpIt, ovpEnd;
    ovpEnd = overlaps.end();
    ovpIt =  overlaps.begin();
    
    //
    // This creates a bunch of dangling pointers, however this is ok,
    // because the entire set is erased after this loop.
    //
    for (; ovpIt != ovpEnd; ovpIt++) {
      Overlap* overlap = (*ovpIt);
      assert(overlap->targetRead != NULL);
      overlap->targetRead->RemoveInterval(overlap->targetBeginNode, overlap->targetEndNode);
      bool result;
      result = overlap->targetRead->RemoveOverlapsWith(this);
      assert(result == true); 
      delete overlap;
    }
    overlaps.clear();
    overlapsWith.clear();
  }

  void RemoveIntervals() {
    //
    // Removes all intervals.
    // 
    ReadNodeSet::iterator intvIt, intvEnd;
    intvIt = orderedReadNodes.begin();
    intvEnd = orderedReadNodes.end();
    for (; intvIt != intvEnd; ++intvIt) {
      ReadNode *nodePtr;
      nodePtr=  (*intvIt);
      assert(nodePtr != NULL);
      nodePtr->RemoveFromSuperset();
      delete nodePtr;
    }
    orderedReadNodes.clear();
  }

  void ClearReadSupsersetFlags() {
    ReadNodeSet::iterator intvIt, intvEnd;
    intvIt = orderedReadNodes.begin();
    intvEnd = orderedReadNodes.end();
    for (; intvIt != intvEnd; ++intvIt) {
      if ((*intvIt)->super != NULL) {
        (*intvIt)->super->ClearFlags();
      }
    }
  }

  void MergeOverlaps() {
    
    ReadNodeSet::iterator intvIt, intvEnd;
    intvEnd = orderedReadNodes.end();
    for ( intvIt = orderedReadNodes.begin(); intvIt != intvEnd; ++intvIt) {
      ReadNode *readNode = *intvIt;

      int overlapIndex;
      //
      // For every read that overlaps with this one, starting at the
      // node from intvIt, make sure they have corresponding ReadNodes
      // and that the corresponding read nodes share the same super node.
      //
      // If the full pairwise alignment between the two reads existed,
      // that should be used to map from this onto the overlap read.
      // The full overlaps take a lot of space, so the initial try is
      // to not use the full overlaps and instead just the coordinates
      // of the overlaps.
      //
      for (overlapIndex = 0; overlapIndex < readNode->overlaps.size(); overlapIndex++) {
        
        //
        // These variables define the endpoints of the overlap.
        //
        ReadNode *thisEndNode, *overlapBeginNode, *overlapEndNode;
        ReadNodeSet::iterator thisReadNodeIt, thisReadNodeItEnd, 
          overlapReadNodeIt, overlapReadNodeItEnd,
          thisNextIt, overlapNextIt;
        //
        // The overlap ends in this read at 'thisEndNode'
        //
        thisEndNode = readNode->overlaps[overlapIndex].thisEnd;
        
        //
        // Assign the node boundaries of the overlapping read.
        //
        overlapBeginNode = readNode->overlaps[overlapIndex].overlapBegin;
        overlapEndNode   = readNode->overlaps[overlapIndex].overlapEnd;
        cout << id << " overlaps with " << overlapBeginNode->read->id << endl;
        Read *overlappingRead = overlapBeginNode->read;
        assert(overlappingRead != NULL);
        

        bool result;
        // 
        // Set up the iterators to pass through this read.
        //
        result = FindNodeItAfter(readNode->pos, thisReadNodeIt);  assert(result);
        result = FindNodeItAt(thisEndNode->pos, thisReadNodeItEnd);     assert(result);
        assert(overlapBeginNode->read != NULL);
        assert(overlapEndNode->read  == overlapBeginNode->read);
        
        cout << "read " << overlappingRead->id << " nodes: " << endl;
        overlappingRead->DumpNodes();

        cout << "This read " << id << " nodes: " << endl;
        DumpNodes();

        result = overlappingRead->FindNodeItAt(overlapBeginNode->pos, overlapReadNodeIt);  assert(result);
        result = overlappingRead->FindNodeItAt(overlapEndNode->pos, overlapReadNodeItEnd);  assert(result);

        //
        // Build a reverse map from super nodes to readnodes.  This
        // allows us to make sure that any read nodes that are being
        // added do not already connect to this read later on
        //

        map<ReadNodeSuperset*, ReadNode*>  thisReadSuperMap, overlapReadSuperMap;
        ReadNodeSet::iterator thisReadIt, overlapReadIt;
        thisReadIt = intvIt;
        if (intvIt != orderedReadNodes.end() and intvIt != thisReadNodeItEnd) {
          //
          // The first node (corresponding to intvIt) is created due
          // to an overlap with overlapRead, and so by definition the
          // two must share a supernode.  Because of that, skip past
          // the first node.
          // 
          ++thisReadIt;
          for (; thisReadIt != thisReadNodeItEnd; ++thisReadIt) {
            if ((*thisReadIt)->super != NULL) {
              thisReadSuperMap[(*thisReadIt)->super] = (*thisReadIt);
            }
          }
        }
        overlapReadIt = overlapReadNodeIt;
        if (overlapReadNodeIt != overlapRead->orderedReadNodes.end() and
            overlapReadNodeIt != overlapReadNodeItEnd) {
          ++overlapReadIt;
          for (; overlapReadIt != overlapReadNodeItEnd; ++overlapReadIt) {
            if ((*overlapReadIt)->super != NULL) {
              overlapReadSuperMap[(*overlapReadIt)->super] = (*overlapReadIt);
            }
          }
        }

        cout << " this (" << this->id << ") super map has " << thisReadSuperMap.size() 
             << " overlap (" << overlapReadIt->id << ") has " << overlapReadSuperMap.size() << endl;
        
        vector<ReadNode*> thisNewReadNodes, overlapNewReadNodes;
        
        UInt thisStartPos = (*intvIt)->pos;
        UInt overlapStartPos = (*overlapReadNodeIt)->pos;
        map<ReadNodeSuperset*, ReadNode*>::iterator thisSuperMapIt, overlapSuperMapIt;

        while (thisReadNodeIt != thisReadNodeItEnd or
               overlapReadNodeIt != overlapReadNodeEnd) {

          // Set up iterators that point at the nodes after the current node.
          thisNextIt = thisReadNodeIt;  ++thisNextIt;
          overlapNextIt = overlapReadNodeIt;  ++overlapNextIt;
          
          thisSuperMapIt = thisReadSuperMap.find((*overlapNextIt)->super);

          ReadNodeSet::iterator thisTmpEnd, overlapTmpEnd;

          if (thisSuperMapIt != thisReadSuperMap.end()) {
            thisTmpEnd = FindReadNodeItAt(thisSuperMapIt->second->pos);
          }
          else {
            thisTmpEnd = thisReadNodeItEnd;
          }

          overlapSuperMapIt = overlapReadSuperMap.find((*thisNextIt)->super);
          if (overlapSuperMapIt != overlapReadSuperMap.end()) {
            overlapTmpEnd = overlapRead->FindReadNodeItAt(overlapSuperMapIt->second->pos);
          }
          else {
            overlapTmpEnd = overlapReadNodeEnd;
          }

          while (thisReadNodeIt    != joinedNodeIt or
                 overlapReadNodeIt != overlapReadNodeEnd) {




          ReadNodeSuperset *nextSuper;
          nextSuper = (*thisNextIt)->super;
          cout << "adding for " << (*thisReadNodeIt)->id << ", " << (*thisReadNodeIt)->pos << "    ovp: " 
               << (*overlapReadNodeIt)->id << ", " << (*overlapReadNodeIt)->pos << endl;


            

          }

          if (nextSuper != NULL and nextSuper->ContainsNode(*overlapReadNodeIt) == true) {
            //
            // The nodes are already connected, simply advance.
            //
            ++thisReadNodeIt;
            ++overlapReadNodeIt;
          }
          else {
            //
            // The next two nods do not belong to the same superset.
            // There are a few cases here: either the nodes are
            // equidistant, in which case the supernodes are simply
            // merged, and no new nodes are created.   Otherwise, the
            // one of the reads needs a corresponding node to be
            // created, and then merged into the superset of the other
            // read node.
            // 
            UInt thisIntvLength = (*thisNextIt)->pos - (*thisReadNodeIt)->pos;
            UInt overlapIntvLength = (*overlapNextIt)->pos - (*overlapReadNodeIt)->pos;

            if (thisIntvLength == overlapIntvLength) {
              //
              // The two intervals are of the same length, but for some reason
              // have not been merged into the same supernode, do that here.
              //
              (*intvIt)->AbsorbSuperset(*overlapReadNodeIt);
              ++thisReadNodeIt;
              ++overlapReadNodeIt;
            }
            else if (thisIntvLength < overlapIntvLength) {
              ReadNode *newOverlapNode = new ReadNode;
              newOverlapNode->pos = (*overlapReadNodeIt)->pos + thisIntvLength;
              (*overlapReadNodeIt)->AddNodeToSuperset(newOverlapNode);
              overlapNewReadNodes.push_back(newOverlapNode);
              ++thisReadNodeIt;
            }
            else {
              assert(overlapIntvLength < thisIntvLength);
              ReadNode *newNode = new ReadNode;
              newNode->pos = (*thisReadNodeIt)->pos + overlapIntvLength;
              (*overlapReadNodeIt)->AddNodeToSuperset(*intvIt);
              thisNewReadNodes.push_back(newNode);
              ++overlapReadNodeIt;
            }
          }
        } // end looping over nodes.
        //
        // At this point, have merged nodes on this read, and overlap
        // read. These need to be added back to the read for
        // bookkeeping and MSA's later on.
        //
        cout << "node " << (*intvIt)->id << " at " << (*intvIt)->pos << " inserting " << thisNewReadNodes.size() << endl;
        cout << "overlap node " << (*overlapReadNodeIt)->id << " at " << (*overlapReadNodeIt)->pos << " inserting " << overlapNewReadNodes.size() << endl;
        orderedReadNodes.insert(thisNewReadNodes.begin(), thisNewReadNodes.end());
        (*overlapReadNodeIt)->read->orderedReadNodes.insert(overlapNewReadNodes.begin(), overlapNewReadNodes.end());
      }
    }
  }
  
  int CollectNextSuperset(ReadNodeSuperset *super, set<ReadNodeSuperset*> &nextSuper) {

    ReadNodeSuperset::iterator readNodeIt, readNodeEnd;
    readNodeEnd = super->end();

    //
    // Find all of the superset nodes that reads go to after this one.
    // This is consistent if there is only one next supernode.
    //

    //      int CollectDestinationSupersets(
    for (readNodeIt = super->begin(); readNodeIt != readNodeEnd; ++readNodeIt) {
      ReadNode *readNodeNext;
      if ((*readNodeIt)->read->FindNodeAfter((*readNodeIt)->pos + 1, readNodeNext)) {
        //          cout << "fnn for " << (*readNodeIt)->pos << " found " << readNodeNext->pos<<endl;
        //          cout << "super for " << super->id << " is " << readNodeNext->id << endl; 
        assert(readNodeNext != NULL);
        if (readNodeNext->super != NULL) {
          nextSuper.insert(readNodeNext->super);
        }
      }
    }        
    return nextSuper.size();
  }

  
  int CountConsistentSupersets() {
    ReadNodeSet::iterator intvIt, intvEnd;
    intvIt = orderedReadNodes.begin();
    intvEnd = orderedReadNodes.end();
    int numConsistent = 0;

    //
    // Examine the superset for all read nodes in this read.
    //
    for (; intvIt != intvEnd; ++intvIt) {
      ReadNodeSuperset *super = (*intvIt)->super;
      if (super == NULL or super->size() <= 1 or super->traversed) {
        continue;
      }
      //
      // Mark this node so that it is not reused.
      //
      super->traversed = true;
      set<ReadNodeSuperset*> nextSuper;
      int numNextSuper = CollectNextSuperset(super, nextSuper);

      if (numNextSuper == 1) {
        numConsistent++;
      }
      else {

        ReadNodeSet::iterator nextReadNodeIt = intvIt;
        ++nextReadNodeIt;
        int readNodeLength = 0;
        if (nextReadNodeIt != intvEnd) {
          readNodeLength = (*nextReadNodeIt)->pos - (*intvIt)->pos;
        }
        
        cout << "super " << super->id << " of len " << readNodeLength << " has " << nextSuper.size() << " destinations. " << endl;
        set<ReadNodeSuperset*>::iterator nextSuperIt, nextSuperEnd;
        nextSuperEnd = nextSuper.end();
        for (nextSuperIt = nextSuper.begin(); nextSuperIt != nextSuperEnd; ++nextSuperIt) {
          set<ReadNodeSuperset*> nextSuperDest;
          CollectNextSuperset(*nextSuperIt, nextSuperDest);
          set<ReadNodeSuperset*>::iterator nextSuperDestIt, nextSuperDestEnd;
          nextSuperDestEnd = nextSuperDest.end();
          for ( nextSuperDestIt = nextSuperDest.begin();
                nextSuperDestIt != nextSuperDestEnd;
                ++nextSuperDestIt ) {
            cout << "  " << super->id << " -> " << (*nextSuperIt)->id << " -> " << (*nextSuperDestIt)->id << endl;
          }
        }
      }
    }
    return numConsistent;
  }
    
  void RemoveIntervalsAndOverlaps() {
    //    cout << "removing intervals and overlaps for " << id << endl;
    RemoveOverlaps();
    RemoveIntervals();
  }
  
  bool FindNextSupersetContainingRead(Read *read, ReadNodeSet::iterator &it, ReadNodeSet::iterator &end) {
    while (it != end) {
      if ((*it)->super != NULL and (*it)->SupersetContainsRead(read)) {
        return true;
      }
      ++it;
    }
    return false;
  }


  void AddInterval(UInt begin, UInt end, ReadNode* &beginNode, ReadNode *&endNode) {
    //
    // Locate either an existing node point on the read, or the one
    // before where this will be inserted. 
    //
    ReadNode *nodeBeforeBegin, *nodeAfterBegin;
    FindNodeAtOrBefore(begin, nodeBeforeBegin);
    FindNodeAfter(begin, nodeAfterBegin);
    
    assert(nodeBeforeBegin != NULL);
    assert(nodeAfterBegin != NULL);

    ReadNode *nodeBeforeEnd, *nodeAfterEnd;
    FindNodeAtOrBefore(end, nodeBeforeEnd);
    FindNodeAfter(end, nodeAfterEnd);

    assert(nodeBeforeBegin != NULL);
    assert(nodeAfterEnd != NULL);

    // 
    // An existing node has been found, will possibly shorten it. 
    // 
    if (nodeBeforeBegin->pos == begin) { 
      beginNode = nodeBeforeBegin; 
    } 
    else { 
      // 
      // Will need to divide the interval, and fix the lengths. 
      // 
      beginNode = new ReadNode(begin, this); 
      orderedReadNodes.insert(beginNode); 
    }

    if (nodeAfterEnd->pos == end) { 
      endNode = nodeAfterEnd; 
    } 
    else { 
      endNode = new ReadNode(end, this); 
      orderedReadNodes.insert(endNode); 
    } 

    // 
    // Now increment the coverage of this interval. 
    // 
    ReadNodeSet::iterator beginIt, endIt; 
    ReadNode findBeginQuery(begin); 
    ReadNode findEndQuery(end); 
    beginIt = orderedReadNodes.find(&findBeginQuery); 
    endIt   = orderedReadNodes.find(&findEndQuery); 
  
    (*endIt)->coverage = (*beginIt)->coverage = nodeBeforeBegin->coverage; 
  
    while(beginIt != endIt) { 
      (*beginIt)->coverage++; 
      ++beginIt; 
    } 
  } 

  ReadNodeSet::iterator GetNodesEndIt() {
    return orderedReadNodes.end();
  }

  bool FindNodeItAfter(UInt pos, ReadNodeSet::iterator &it) {
    ReadNode query(pos);
    it = orderedReadNodes.lower_bound(&query); 
    if (it == orderedReadNodes.end()) { 
      return false; 
    } 
    else { 
      return true; 
    } 
  }

  bool FindNodeAfter(UInt pos, ReadNode *&readNode) { 
    ReadNodeSet::iterator it; 
    if (FindNodeItAfter(pos, it)) {
      readNode = *it;
      return true;
    }
    else {
      readNode = NULL;
      return false;
    }
  } 
   
  bool FindNodeItAt(UInt pos, ReadNodeSet::iterator &it) {
    ReadNode query(pos);
    it = orderedReadNodes.find(&query); 
    if (it != orderedReadNodes.end()) { 
      return true; 
    } 
    else { 
      return false; 
    } 

  }
  bool FindNodeAt(UInt pos, ReadNode *&readNode) { 
    ReadNodeSet::iterator it; 
    if ( FindNodeItAt(pos, it) ) {
      readNode = *it;
      return true;
    }
    else {
      readNode = NULL;
      return false;
    }
  } 
  void DumpNodes() { 
    ReadNodeSet::iterator it, end; 
    cout << "read of length " << readLength << endl; 
    if (orderedReadNodes.size() == 0) { 
      cout << "  No nodes." << endl; 
      return; 
    } 
    for (it = orderedReadNodes.begin(), 
           end = orderedReadNodes.end();  
         it != end;  
         ++it) { 
      cout << (*it)->pos << " "; 
    } 
    cout << endl; 
  } 

  bool FindNodeAtOrBefore(UInt pos, ReadNode *&readNode) { 
    ReadNodeSet::iterator it; 
    ReadNode query(pos); 
    it = orderedReadNodes.lower_bound(&query); 
    if (it == orderedReadNodes.end()) { 
      readNode = NULL; 
      return false; 
    } 
    if ((*it)->pos == pos) { 
      readNode = (*it); 
      return true; 
    } 
    else { 
      --it; 
      readNode = *it; 
      assert((*it)->pos != pos); 
      return true; 
    } 
  } 
  
  void PrintLowCoverageIntervals(int coverageCutoff) { 
    ReadNodeSet::iterator it, end, next; 
    end = orderedReadNodes.end(); 
    for (it = orderedReadNodes.begin(); it != end; ++it) { 
      next = it; 
      next++; 
      if ((*it)->coverage < coverageCutoff) { 
        UInt length; 
        if (next != end) { 
          length = (*next)->pos - (*it)->pos; 
        } 
        else { 
          length = 0; 
        } 
        cout << (*it)->pos << "\t" << length << "\t" << (*it)->coverage << "\t" << (*it)->id; 
        if ((*it)->super != NULL) { 
          cout << "\t" << (*it)->super->size(); 
        } 
        cout << endl; 
      } 
    } 
  } 

  bool IsLikelyChimeric(int coverageCutoff) { 
  
    if (orderedReadNodes.size() <= 1) { 
      // 
      // Empty, can't say anything about this. 
      // 
      return false; 
    } 
    ReadNodeSet::iterator it, end, next; 
    end = orderedReadNodes.end(); 
    ReadNodeSet::iterator maxFromBegin, maxFromEnd; 
    int totalCoverage = 0; 
    for (it = orderedReadNodes.begin(); it != end; ++it) { 
      ++(next = it); 
      if (next != end) { 
        int length =  (*next)->pos - (*it)->pos; 
        totalCoverage += length * (*it)->coverage; 
      } 
    } 
    // 
    // Move up in coverage from left to right. 
    // 
    ReadNodeSet::iterator prev, begin; 
    prev = it = orderedReadNodes.begin(); 
    it++; 
    while (it != end and (*it)->coverage >= (*prev)->coverage) { 
      ++it; 
      ++prev; 
    } 
    maxFromBegin = prev; 

    // 
    // Move up in coverage from right to left. 
    //
    it = orderedReadNodes.end();
    --it;
    prev = it;

    maxFromEnd = it;
    --maxFromEnd;
    if (prev == maxFromBegin) {
      maxFromEnd = prev;
    }
    else {
      while (maxFromEnd != maxFromBegin and (*maxFromEnd)->coverage >= (*prev)->coverage) {
        prev = maxFromEnd;
        --maxFromEnd;
      }
      maxFromEnd = prev;
    }
    
    // 
    // Single peak profiles are ok.
    //
    if ( maxFromBegin == maxFromEnd ) {
      return false;
    }
    //
    // Now look for minimum coverage between.
    //
    int minCoverage = (*maxFromEnd)->coverage;
    if ((*maxFromBegin)->coverage < (*maxFromEnd)->coverage ) {
      minCoverage = (*maxFromBegin)->coverage;
    }
    for (it = maxFromBegin; it != maxFromEnd; ++it) {
      if ( (*it)->coverage < minCoverage ) {
        minCoverage = (*it)->coverage;
      }
    }
    if (minCoverage < coverageCutoff) {
      cout << "min coverage " << minCoverage << endl;
      cout << "max from begin " << (*maxFromBegin)->coverage << " " << (*maxFromBegin)->id << endl;
      cout << "max from end   " << (*maxFromEnd)->coverage << " " << (*maxFromEnd)->id << endl;
      return true;
    }
    else {
      return false;
    }
  }
  
  bool RemoveIfLowCoverage(int minCoverage) {
    ReadNodeSet::iterator it, end;
    end = orderedReadNodes.end();
    if (orderedReadNodes.size() <= 1) {
      //
      // Do not do anything with singleton reads.
      //
      return false;
    }
    bool allLowCoverage = true;
    UInt maxCoverage = 0;
    for (it = orderedReadNodes.begin(); it != end; it++) {
      UInt coverage = (*it)->coverage;
      maxCoverage = max(maxCoverage, coverage);
      if (coverage >= minCoverage) {
        allLowCoverage = false;
        break;
      }
    }

    if (allLowCoverage == true) {
      //      cout << "removing read with max coverage " << maxCoverage << endl;
      RemoveIntervalsAndOverlaps();
      return true;
    }
    return false;
  }

  bool RemoveLowCoverage(int minCoverage) {
    
    ReadNodeSet::iterator it, end;
    end = orderedReadNodes.end();
    if (orderedReadNodes.size() <= 1) {
      //
      // Do not do anything with singleton reads.
      //
      return false;
    }
    for (it = orderedReadNodes.begin(); it != end; it++) {
      int coverage = (*it)->coverage;
      if (coverage < minCoverage) {
        //
        // Delete this node, and any overlaps it is part of.
        //
        RemoveIntervalsAndOverlaps();
        return true;
      }
    }
  }
};

UInt Read::counter = 0;

#endif
