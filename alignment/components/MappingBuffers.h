#ifndef ALIGNMENT_COMPONENTS_MAPPING_BUFFERS_H_
#define ALIGNMENT_COMPONENTS_MAPPING_BUFFERS_H_

#include "Enumerations.h"
#include "tuples/DNATuple.h"
#include "datastructures/anchoring/WeightedInterval.h"


//
// Define a list of buffers that are meant to grow to high-water
// marks, and not shrink down past that.   The memory is reused rather
// than having multiple calls to new.
//
class MappingBuffers {
public:
  vector<int> hpInsScoreMat, insScoreMat;
  vector<int> kbandScoreMat;
  vector<Arrow> hpInsPathMat, insPathMat;
  vector<Arrow> kbandPathMat;
  vector<int>   scoreMat;
  vector<Arrow> pathMat;
  vector<ChainedMatchPos> matchPosList;
  vector<ChainedMatchPos> rcMatchPosList;
  vector<BasicEndpoint<ChainedMatchPos> > globalChainEndpointBuffer;
  vector<Fragment> sdpFragmentSet, sdpPrefixFragmentSet, sdpSuffixFragmentSet;
  TupleList<PositionDNATuple> sdpCachedTargetTupleList;
  TupleList<PositionDNATuple> sdpCachedTargetPrefixTupleList;
  TupleList<PositionDNATuple> sdpCachedTargetSuffixTupleList;
  std::vector<int> sdpCachedMaxFragmentChain;
  vector<double> probMat;
  vector<double> optPathProbMat;
  vector<float>  lnSubPValueMat;
  vector<float>  lnInsPValueMat;
  vector<float>  lnDelPValueMat;
  vector<float>  lnMatchPValueMat;
  vector<int>    clusterNumBases;
  ClusterList    clusterList;
  ClusterList    revStrandClusterList;
};

#endif
