#ifndef DATASTRUCTURES_SAF_ALN_GROUP_H_
#define DATASTRUCTURES_SAF_ALN_GROUP_H_

class AlnGroup {
 public:
	vector<unsigned int> id;
	vector<string>       path;
	int FindPath(int idKey, string &val) {
		int i;
		for (i = 0; i < id.size(); i++) {
			if (idKey == id[i]) {
				val = path[i];
				return 1;
			}
		}
		return 0;
	}
};

#endif
