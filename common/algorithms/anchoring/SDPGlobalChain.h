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
			sdpFragments[s].weight = tupleSize;
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
		sdpi = maxFragmentChain[i];
		while (f < nFragments and 
					 (fragments[f].GetX() != sdpFragments[sdpi].x or fragments[f].GetY() != sdpFragments[sdpi].y)) {
			f++;
		}
		/*
		cout << "f: "<< f << " " << nFragments << endl;
		cout << "match: " << fragments[f].GetX() << " " << fragments[f].GetY() << " " << sdpFragments[sdpi].x << " " << sdpFragments[sdpi].y << " " << fragments[f].GetLength() << " " << sdpFragments[sdpi].length << endl;
		*/
		while (f < nFragments and 
					 i < maxFragmentChain.size() and 
					 sdpFragments[maxFragmentChain[i]].x + tupleSize <= fragments[f].GetX() + fragments[f].GetLength()) {
			i++;
		}
		//		cout << "i: " << i << " " << maxFragmentChain.size() << endl;
		if (f < nFragments) {
			optFragmentChainIndices.push_back(f);
			chainLength += fragments[f].GetLength();
		}
	}
	//	cout << "done with " << chainLength << endl;
	return chainLength;
}



#endif
