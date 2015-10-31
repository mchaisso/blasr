#ifndef SDP_GLOBAL_CHAIN_H_
#define SDP_GLOBAL_CHAIN_H_

#include "algorithms/alignment/sdp/SparseDynamicProgramming.h"
#include "algorithms/alignment/sdp/SDPFragment.h"

class IndexedFragment: public Fragment {
 public:
	int i;
};


template<typename T_Fragment>
	int SDPGlobalChain( T_Fragment *fragments, 
											DNALength nFragments, 
											vector<VectorIndex> &optFragmentChainIndices, int tupleSize, vector<Fragment> &sdpFragments2) {



	//
	// First count all tuples, making sure they are bigger than the sdp
	//
	int i, j, s;
	int nSdpFragments = 0;
	vector<IndexedFragment> sdpFragments;
	sdpFragments.resize(nFragments);
	s = 0;
	int maxX = 0;
	for (i = 0; i < nFragments; i++) {
		sdpFragments[i].i = i;
		sdpFragments[i].x = fragments[i].GetX();
		sdpFragments[i].y = fragments[i].GetY();
		sdpFragments[i].length = tupleSize;
		sdpFragments[i].weight = 0;

		if ((int)sdpFragments[i].x >= maxX) {
			maxX = sdpFragments[i].x;
		}
	}
	std::sort(sdpFragments.begin(), sdpFragments.end(), LexicographicFragmentSort<Fragment>());
	for (s = 0; s < sdpFragments.size(); s++) {
		sdpFragments[s].index = s;
		sdpFragments[s].chainPrev = 0;
		sdpFragments[s].above = -1;
	}


	//
	// Assume inversion will be in rc max frament chain set.
	//
	std::vector<int> maxFragmentChain;
	SDPLongestCommonSubsequence(maxX, sdpFragments, tupleSize, 2, 2, -5, maxFragmentChain, Global);


	s = 0;
	int f = 0;
	int sdpi = 0;
	int chainLength = 0;
	i = 0;
	bool print = maxFragmentChain.size() > 20000;
	for (i = 0; i < maxFragmentChain.size(); i++) {
		optFragmentChainIndices.push_back(sdpFragments[maxFragmentChain[i]].i);
	}
	chainLength = maxFragmentChain.size()*tupleSize;

	/*	while (f < nFragments and i < maxFragmentChain.size()) {

		sdpi = maxFragmentChain[i];

		while ( f < nFragments and 
					 (fragments[f].GetX() != sdpFragments[sdpi].x or fragments[f].GetY() != sdpFragments[sdpi].y)) {
			f++;
		}

		if (print) {
			cout << f << "\t" << sdpi << "\t" << fragments[f].GetX() << "\t" << sdpFragments[sdpi].x << "\t" << fragments[f].GetY() << "\t" << sdpFragments[sdpi].y << endl;
		}
			
		if (f < nFragments) {
			optFragmentChainIndices.push_back(f);
			fragments[f].l = tupleSize;
			chainLength += fragments[f].l;
		}
		i++;
	}
	*/
	cout << "mfc: "<< maxFragmentChain.size() << " cl: " << chainLength << endl;
	//	cout << "done with " << chainLength << endl;
	return chainLength;
}



#endif
