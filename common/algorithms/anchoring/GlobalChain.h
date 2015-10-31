#ifndef GLOBAL_CHAIN_H_
#define GLOBAL_CHAIN_H_
#include <vector>
#include "PrioritySearchTree.h"


using namespace std;
template<typename T_Fragment,typename T_Endpoint>
	void FragmentSetToEndpoints(T_Fragment *fragments, int nFragments, vector<T_Endpoint> &endpoints) {
	
	endpoints.resize(nFragments*2);
	
	int i;
	int ep = 0;
	for (i = 0; i < nFragments; i++) {
		endpoints[ep].FragmentPtrToStart(&fragments[i]);
		ep++;
		endpoints[ep].FragmentPtrToEnd(&fragments[i]);
		ep++;
	}
}


template<typename T_Fragment>
UInt RestrictedGlobalChain( T_Fragment *fragments, 
													 DNALength nFragments, 
													 float maxIndelRate,
														vector<VectorIndex> &optFragmentChainIndices,
														vector<UInt> &scores,
														vector<UInt> &prevOpt) {
	// assume fragments are sorted by t
	
	UInt f1, f2;
	scores.resize(nFragments);
	prevOpt.resize(nFragments);
	std::fill(scores.begin(), scores.end(), 0);

	UInt globalOptScore = 0;
	UInt globalOptIndex = 0;
	for (f1 = 0; f1 < nFragments; f1++) {
		prevOpt[f1] = f1;
		scores[f1] = 1;
	}
	
	for (f1 = 0; f1 < nFragments - 1; f1++) {
		UInt maxF1Score = 0;
		UInt maxF1Prev  = f1;
		for (f2 = f1+1; f2 < nFragments; f2++ ){ 
			//
			//  Check to see if the fragments may be connected within the
			//  expected indel rate.
			//
			if (fragments[f2].GetQ() > fragments[f1].GetQ() + fragments[f1].GetQW() and
					fragments[f2].GetT() > fragments[f1].GetT() + fragments[f1].GetQW()) {
				//
				// Compute drift from diagonal.
				//
				UInt tDiff, qDiff;
				tDiff = fragments[f2].GetT() - (fragments[f1].GetT() + fragments[f1].GetQW());
				qDiff = fragments[f2].GetQ() - (fragments[f1].GetQ() + fragments[f1].GetQW());
				UInt tIns, qIns;
				tIns = qIns = 0;
				UInt maxDiff = max(tDiff, qDiff);
				UInt minDiff = min(tDiff, qDiff);
				if (maxDiff - minDiff < minDiff*maxIndelRate) {
					//
					// The fragment is sufficiently close to the diagonal to
					// consider it as a chain.  
					//
					if (scores[f2] < scores[f1] + 1) {
						scores[f2] = scores[f1] + 1;
						prevOpt[f2] = f1;
						if (scores[f2] > globalOptScore) {
							globalOptScore = scores[f2];
							globalOptIndex = f2;
						}
					}
				}
			}
		}
	}
	UInt index = globalOptIndex;
	UInt prevIndex;
	while(index != prevOpt[index]) {
		optFragmentChainIndices.push_back(index);
		assert(optFragmentChainIndices.size() < nFragments);
		prevIndex = index;
		index = prevOpt[index];
		// Make sure there was no problem with backtracking.
		assert(index < nFragments);
		assert(index <= prevOpt[prevIndex]);
	}
	optFragmentChainIndices.push_back(index);
		std::reverse(optFragmentChainIndices.begin(), optFragmentChainIndices.end());
	return optFragmentChainIndices.size();
}
	

template<typename T_Fragment, typename T_Endpoint>
	int GlobalChain( T_Fragment *fragments, 
									 DNALength nFragments, 
									 vector<VectorIndex> &optFragmentChainIndices,
									 vector<T_Endpoint> *bufEndpointsPtr = NULL) {
	

	//
	// Initialize the fragment score to be the length of each fragment.
	//
	if (nFragments == 0) {
		return 0;
	}

	DNALength f;
	for (f = 0; f < nFragments; f++) { 
		fragments[f].SetScore(fragments[f].GetLength());
	}

	//
	// Add the start/end points of each fragment. This allows separate scoring
	// of start points and activation of endpoints.
	//
	vector<T_Endpoint> endpoints;
	vector<T_Endpoint> *endpointsPtr;
	if (bufEndpointsPtr != NULL) {
		endpointsPtr = bufEndpointsPtr;
	}
	else {
		endpointsPtr = &endpoints;
	}
				
	
	FragmentSetToEndpoints<T_Fragment, T_Endpoint>(fragments, nFragments, *endpointsPtr);

	//
	// The Starting points of all the fragmements are in order, 
	// but not necessarily all of the end endpoints, so
	// the list must be resorted.
	//
	std::sort(endpointsPtr->begin(), endpointsPtr->end(), typename T_Endpoint::LessThan());
	
	PrioritySearchTree<T_Endpoint> pst;

	pst.CreateTree(*endpointsPtr);

	VectorIndex p;
	VectorIndex maxScoringEndpoint = 0;
	bool maxScoringEndpointFound = false;
	for (p = 0; p < endpointsPtr->size(); p++) {
		int x = (*endpointsPtr)[p].p.GetX();
		int y = (*endpointsPtr)[p].p.GetY();
		if ((*endpointsPtr)[p].GetSide() == T_Endpoint::Start) {
			int maxPointIndex;
			if (pst.FindIndexOfMaxPoint((*endpointsPtr), (*endpointsPtr)[p].GetKey(), maxPointIndex)) {
				(*endpointsPtr)[p].SetChainPrev((*endpointsPtr)[maxPointIndex].GetFragmentPtr());
				assert((*endpointsPtr)[maxPointIndex] < (*endpointsPtr)[p]);
				//				assert((*endpointsPtr)[p].GetY() <= (*endpointsPtr)[p].GetY());
				int score = (*endpointsPtr)[maxPointIndex].GetScore() + (*endpointsPtr)[p].GetScore();
				(*endpointsPtr)[p].SetScore(score);
			}
			else {
				(*endpointsPtr)[p].SetChainPrev(NULL);
			}
		}	else {
			assert((*endpointsPtr)[p].GetSide() == T_Endpoint::End);
			// 
			// The score of the fragment should be already set.  So simply activate
			// it here (make the point be visible in a search).
			//
			pst.Activate((*endpointsPtr), p);
			if (maxScoringEndpointFound == false or
					(*endpointsPtr)[maxScoringEndpoint].GetScore() < (*endpointsPtr)[p].GetScore()) {
				maxScoringEndpoint = p;
				maxScoringEndpointFound = true;
			}
		}
	}
	

	
	// 
	// Now compute the chain of optimum fragments
	//
	T_Fragment *optFragmentPtr;
	if (maxScoringEndpointFound == false) {
		// 
		// Null case, no endpoints have been processed.
		//
		return 0;
	}
	
	optFragmentPtr = (*endpointsPtr)[maxScoringEndpoint].GetFragmentPtr();

	unsigned int numIter = 0;
	while (optFragmentPtr != NULL) {
		optFragmentChainIndices.push_back((int) (optFragmentPtr - &fragments[0]));
		optFragmentPtr = (T_Fragment*) optFragmentPtr->GetChainPrev();
		// 
		// Do a sanity check to make sure this loop is finite -- the optimal
		// fragment chain should never contain more fragments than what are
		// input.
		//
		assert(numIter < nFragments);
		++numIter;
	}
  reverse(optFragmentChainIndices.begin(), optFragmentChainIndices.end());
	return optFragmentChainIndices.size();

}


template<typename T_Fragment, typename T_Endpoint>
int GlobalChain(vector<T_Fragment> &fragments, vector<VectorIndex> &optFragmentChainIndices) {
	return GlobalChain<T_Fragment, T_Endpoint>(&fragments[0], fragments.size(), optFragmentChainIndices);
}

template<typename T_Fragment, typename T_Endpoint>
	int GlobalChain(vector<T_Fragment> &fragments, DNALength start, DNALength end, vector<VectorIndex> &optFragmentChainIndices,vector<T_Endpoint> *bufEndpointsPtr = NULL) {
	return GlobalChain<T_Fragment, T_Endpoint>(&fragments[start], end - start, optFragmentChainIndices, bufEndpointsPtr);
}



#endif
