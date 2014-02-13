#ifndef BWT_H_
#define BWT_H_

// 
// Implement fast quicksorting on suffixes according to Burrows and Wheeler.
//

namespace pbbwt {
	template<typename T_KEY>
	int RadixSort(vector<T_KEY> &V, int keyLength, int alphSize) {
		// 
		// 
		//
		vector<int> bucketCount, bucketCur, prev;
		bucketCount.resize(alphSize);
		prev.resize(alphSize);
		bucketCur.resize(alphSize);

		int c = keyLength - 1;
		//
		// 
		//
		while (c >= 0) {
			std::fill(bucketCount.begin(), bucketCount.end(), 0);
			std::fill(prev.begin(), prev.end(), -1);
			int i;

			//
			// Count the number of elements in every bucket, and link
			// each value to the previous value in every bucket.
			//
			for (i = 0; i < l; i++ ){ 
				V[i] = prev[s[V[i]+c]];
				prev[s[V[i]+c]] = i;
				bucketCount[s[V[i]+c]]++;
			}
			
			// 
			// Restore the elements in V.
			//
					
			std::fill(bucketCur.begin(), bucketCur.end(), 0);
			for (i = 0; i < l; i++ ) {
				

			}
		}
		// 
	}
}


#endif
