#ifndef SDP_GLOBAL_CHAIN_H_
#define SDP_GLOBAL_CHAIN_H_

#include "algorithms/alignment/sdp/SparseDynamicProgramming.h"
#include "algorithms/alignment/sdp/SDPFragment.h"



template<typename T_Fragment>
	int SDPGlobalChain( T_Fragment *fragments, 
											DNALength nFragments, 
											vector<VectorIndex> &optFragmentChainIndices, int tupleSize, vector<Fragment> &sdpFragments) {


	//	vector<Fragment> sdpFragments;
	sdpFragments.clear();

	//
	// First count all tuples.
	//
	int i, j, s;
	int nSdpFragments = 0;
	for (i = 0; i < nFragments; i++) {
		assert(fragments[i].GetLength() >= tupleSize);
		nSdpFragments += 3; //fragments[i].GetLength() - tupleSize + 1;
	}
	sdpFragments.resize(nSdpFragments);
	s = 0;
	int maxX = 0;
	for (i = 0; i < nFragments; i++) {
		for (j = 0; j < 3; j++) { //fragments[i].GetLength() - tupleSize + 1; j++) {
			sdpFragments[s].x = fragments[i].GetX() + j;
			sdpFragments[s].y = fragments[i].GetY() + j;
			sdpFragments[s].length = tupleSize;
			//
			// Use the weight as a proxy for wehre
			if (j > 0) {
				sdpFragments[s].weight = 1;
			}
			else {
				sdpFragments[s].weight = 0;
			}
			if ((int)sdpFragments[s].x >= maxX) {
				maxX = sdpFragments[s].x;
			}
			s++;
		}
	}
	std::sort(sdpFragments.begin(), sdpFragments.end(), LexicographicFragmentSort<Fragment>());
	
	s = 0;
  int fCur = 0;
  while (s + 1 <= sdpFragments.size()) {
    sdpFragments[fCur] = sdpFragments[s];
    while (s < sdpFragments.size() and sdpFragments[fCur].x == sdpFragments[s].x and sdpFragments[fCur].y == sdpFragments[s].y) {
			if (sdpFragments[s].weight == 0) {
				sdpFragments[fCur].weight = 0;
			}
      s++;
    }

    fCur++;

  }
  sdpFragments.resize(fCur);

	for (s = 0; s < sdpFragments.size(); s++) {
		sdpFragments[s].index = s;
		sdpFragments[s].chainPrev = 0;
		sdpFragments[s].above = -1;
	}


	//
	// Assume inversion will be in rc max frament chain set.
	//
	std::vector<int> maxFragmentChain;
	SDPLongestCommonSubsequence(maxX, sdpFragments, tupleSize, 2, 2, -5, maxFragmentChain, Local);

	s = 0;
	int f = 0;
	int sdpi = 0;
	int chainLength = 0;
	i = 0;
	while (f < nFragments and i < maxFragmentChain.size()) {

		while (i < maxFragmentChain.size() and 
					 sdpFragments[maxFragmentChain[i]].weight != 0) {
			i++;
		}
		sdpi = maxFragmentChain[i];
		if (i >= maxFragmentChain.size()) {
			break;
		}
		while ( f < nFragments and 
					 (fragments[f].GetX() != sdpFragments[sdpi].x or fragments[f].GetY() != sdpFragments[sdpi].y)) {
			f++;
		}

		//		cout << "i: " << i << " " << maxFragmentChain.size() << endl;
		if (f < nFragments) {
			optFragmentChainIndices.push_back(f);
			chainLength += fragments[f].GetLength();
		}
		i++;
	}
	//	cout << "done with " << chainLength << endl;
	return chainLength;
}



#endif
