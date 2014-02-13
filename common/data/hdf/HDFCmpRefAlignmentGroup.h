#ifndef DATA_HDF_HDF_CMP_REF_ALIGNMENT_GROUP_H_
#define DATA_HDF_HDF_CMP_REF_ALIGNMENT_GROUP_H_

#include "H5Cpp.h"
#include "data/hdf/HDFData.h"
#include "data/hdf/HDFGroup.h"
#include "data/hdf/HDFCmpExperimentGroup.h"
#include <map>
#include <string>
using namespace std;

class HDFCmpRefAlignmentGroup {
 public:
	HDFGroup  refGroup;
	string refGroupName;
	vector<HDFCmpExperimentGroup*> readGroups;
	HDFAtom<string> annotationStringAtom;
	map<string,int> experimentNameToIndex;
    // A RefAlignmentGroup may contain one or more 
    // ExperimentGroups. The following shows a 
    // RefAlignmentGroup containing two ExperimentGroups.
    // /ref00001/m121219_103658_42194_c000447032559900001500000112311426_s1_p0
    // /ref00001/m121219_103658_42194_c000447032559900001500000112311426_s2_p0
    //
    // But these ExperimentGroups are not necessarily
    // grouped by movie. For example, the following groups
    // can be seen in deep-sorted cmp.h5
    // /ref00001/rg8953-0
    // /ref00001/rg2453-1

	int Initialize(CommonFG &group, string _refGroupName) {
		refGroupName = _refGroupName;
		refGroup.Initialize(group, _refGroupName);
		//		annotationStringAtom.Initialize(refGroup.group, "annotationString");
	}

  void Create(HDFGroup parent, string refGroupNameP) {
    refGroupName = refGroupNameP;
    parent.AddGroup(refGroupName);
    refGroup.Initialize(parent, refGroupName);
  }

  HDFCmpExperimentGroup *GetExperimentGroup(string readGroupName) {
    //
    // In contrast to initialization, only create one group.
    //
    map<string,int>::iterator it = experimentNameToIndex.find(readGroupName);
		if (it != experimentNameToIndex.end()) {
      assert(it->second < readGroups.size());
      return readGroups[it->second];
    }

    //
    // Allocate the new group structure
    //
    int newReadGroupIndex = readGroups.size();
    HDFCmpExperimentGroup* readGroupPtr = new HDFCmpExperimentGroup;
    readGroups.push_back(readGroupPtr);
    experimentNameToIndex[readGroupName] = newReadGroupIndex;

    //
    // Now add it to the cmp.h5 file.
    //
    if (readGroupPtr->Create(refGroup, readGroupName) == 0) {
      delete readGroupPtr;
      experimentNameToIndex[readGroupName] = -1;
      return NULL;
    }
    return readGroupPtr;
  }

  bool ContainsExperimentGroup(string readGroupName) {
    return experimentNameToIndex.find(readGroupName) != experimentNameToIndex.end();
  }

	HDFCmpExperimentGroup* InitializeExperimentGroup(string experimentGroupName, set<string> &includedFields) {
		if (refGroup.ContainsObject(experimentGroupName)) {
			HDFCmpExperimentGroup* newGroup = new HDFCmpExperimentGroup;

			if (newGroup->Initialize(refGroup, experimentGroupName, includedFields) == 0) {
				cout << "ERROR, could not initialize the exp group." << endl;
				exit(1);
			}
			experimentNameToIndex[experimentGroupName] = readGroups.size();
			readGroups.push_back(newGroup);
			return newGroup;
		}
		else {
			return NULL;
		}
	}

  
	HDFCmpExperimentGroup* InitializeExperimentGroup(string experimentGroupName) {
    set<string> EMPTYIncludedFields;
    return InitializeExperimentGroup(experimentGroupName, EMPTYIncludedFields);
	}
};


#endif
