#ifndef ALIGNMENT_SET_REF_GROUP_H_
#define ALIGNMENT_SET_REF_GROUP_H_


#include <map>

#include "ReferenceSequence.h"

template<typename T_ReferenceSequence=FullRefererenceSequence>
class RefGroup {
 public:
 map<string, T_ReferenceSequence*>  referenceMap;

 vector<T_ReferenceSequence> references;
 
 void AddReferenceSequence(const T_ReferenceSequence &refSeq) {
   references.push_back(refSeq);
 }

 void BuildMap() {
   int i;
   for (i = 0; i < references.size(); i++) {
     referenceMap[references[i].sequenceName] = i;
   }
 }
};


#endif
