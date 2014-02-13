#ifndef DATASTRUCTURES_ALIGNMENT_SET_ALIGNMENTSET_H_
#define DATASTRUCTURES_ALIGNMENT_SET_ALIGNMENTSET_H_

#include <vector>
#include <map>
#include "datastructures/alignmentset/ReadGroup.h"
#include "datastructures/alignmentset/ReferenceSequence.h"
#include "datastructures/alignmentset/SAMAlignment.h"
#include "SAMHeader.h"
#include "FASTASequence.h"

template<typename T_ReferenceSequence=SAMReferenceSequence, typename T_ReadGroup=SAMReadGroup, typename T_Alignment=SAMAlignment>
class AlignmentSet {
 public:
  SAMHeader header;
  vector<T_ReferenceSequence> references;
  vector<T_ReadGroup> readGroups;
  vector<T_Alignment> alignments;

  //
  //  Rearrange references such that they are placed in the same order
  //  as fastaReferences
  //  
  void RearrangeReferences(vector<FASTASequence> & fastaReferences) {
      int i = 0;
      map<string, int> fastaRefToIndex;
      map<string, int>::iterator it;
      for (i = 0; i<fastaReferences.size(); i++) {
          it = fastaRefToIndex.find(fastaReferences[i].GetName());
          if (it != fastaRefToIndex.end()) {
              cout<<"Error, reference with name \""<<fastaReferences[i].GetName()
				  <<"\" in the reference genome is not unique"<<endl;
              exit(1);
          }
          fastaRefToIndex[fastaReferences[i].GetName()] = i;
      }
      vector<T_ReferenceSequence> newreferences;
      for (i = 0; i < references.size(); i++) {
          newreferences.push_back(T_ReferenceSequence());
      }
      for (i = 0; i < references.size(); i++) {
          it = fastaRefToIndex.find(references[i].sequenceName);
          if (it == fastaRefToIndex.end()) {
              cout<<"Error, can not find reference name "<<references[i].sequenceName
                  <<" in the reference genome."<<endl;
              exit(1);
          }
          newreferences[(*it).second] = references[i];
      }
      references = newreferences;
  }

};


#endif
