#ifndef DATASTRUCTURES_SAF_HDF_MOVIE_INFO_H_
#define DATASTRUCTURES_SAF_HDF_MOVIE_INFO_H_

#include "Types.h"
#include <string>

using namespace std;

class MovieInfo {
 public:
	vector<string> name;
	vector<UInt>   run;
	vector<UInt>   experiment;
	vector<UInt>   id;
	vector<string> bindingKit;
	vector<string> sequencingKit;
	vector<string> softwareVersion;
	int FindMovie(int idKey, string &nameVal) {
		int i;
		for (i = 0; i < id.size(); i++) {
			if (id[i] == idKey) {
				nameVal = name[i];
				return 1;
			}
		}
		return 0;
	}
};

#endif
