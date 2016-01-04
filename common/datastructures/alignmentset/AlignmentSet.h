#ifndef DATASTRUCTURES_ALIGNMENT_SET_ALIGNMENTSET_H_
#define DATASTRUCTURES_ALIGNMENT_SET_ALIGNMENTSET_H_

#include <vector>
#include <map>
#include <set>
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
  map<string, int> refNameToIndex;
  //
  //  Rearrange references such that they are placed in the same order
  //  as fastaReferences
  //  
  void GetRepresentativeChangelistId(string &changelistId) {
		if (readGroups.size() > 0) {
		  changelistId = readGroups[0].changelistId;
	  } 
		else {
			changelistId = "0";
		}
  }
  void GetRepresentativeSequencingKit(string &sequencingKit) {
		if (readGroups.size() > 0) {
		  sequencingKit = readGroups[0].sequencingKit;
	  } 
		else {
			sequencingKit = "0";
		}
		
	}
  void GetRepresentativeBindingKit(string &bindingKit) {
		if (readGroups.size() > 0) {
		  bindingKit = readGroups[0].bindingKit;
	  } 
		else {
			bindingKit = "0";
		}
	}

  void RearrangeReferences(vector<FASTASequence> & fastaReferences) {
      int i = 0;
      map<string, int> fastaRefToIndex;
      map<string, int>::iterator it;
			set<string> refNames;
			for (i = 0; i < references.size(); i++) {
				refNames.insert(references[i].GetSequenceName());
			}
			int rank = 0;
      for (i = 0; i<fastaReferences.size(); i++) {
          it = fastaRefToIndex.find(fastaReferences[i].GetName());
          if (it != fastaRefToIndex.end()) {
              cout<<"Error, reference with name \""<<fastaReferences[i].GetName()
				  <<"\" in the reference genome is not unique"<<endl;
              exit(1);
          }
					if (refNames.find(fastaReferences[i].GetName()) != refNames.end()) {
						fastaRefToIndex[fastaReferences[i].GetName()] = rank;
						++rank;
					}
      }
      vector<T_ReferenceSequence> newreferences(references.size());
      for (i = 0; i < references.size(); i++) {
          it = fastaRefToIndex.find(references[i].sequenceName);
          if (it == fastaRefToIndex.end()) {
              cout<<"Error, can not find reference name "<<references[i].sequenceName
                  <<" in the reference genome."<<endl;
              exit(1);
          }
          newreferences[(*it).second] = references[i];
					refNameToIndex[references[i].sequenceName] = (*it).second;
      }
      references = newreferences;

  }

};


#endif
