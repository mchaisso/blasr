#include "../common/FASTAReader.h"
#include "../common/FASTASequence.h"
#include "../common/utils.h"

#include <string>
#include <map>
#include <vector>
using namespace std;

void FilterReadName(string in, string &out) {
	int slashPos;
	slashPos = in.find_first_of('/');
	out.assign(in, 0, slashPos);
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: checkfilt reads.fa alignfile" << endl;
		exit(1);
	}

	string readsFileName = argv[1];
	string alignFileName = argv[2];

	FASTAReader reader;
	reader.Initialize(readsFileName);
	
	vector<FASTASequence*> reads;
	FASTASequence read;
	map<string,int> nameIndexMap;
	int readIndex = 0;
	while(reader.GetNext(read)) {
		FASTASequence *readPtr = new FASTASequence;
		*readPtr = read;
		reads.push_back(readPtr);
		string newName;
		FilterReadName(read.title, newName);
		nameIndexMap[newName] = readIndex;
		++readIndex;
	}

	ifstream hitsIn;
	CrucialOpen(alignFileName, hitsIn);
	int nFilt = 0;
	int nHits = 0;
	while(hitsIn) {
		float score;
		int  length;
		float acc;
		string name;
		string alignstr;
		if (!(hitsIn >> score >> length >> acc >> name >> alignstr)) { break;}
		++nHits;
		string filtName;
		FilterReadName(name, filtName);
		if (nameIndexMap.find(filtName) != nameIndexMap.end()) {
			int readIndex = nameIndexMap[filtName];
			bool isFiltered = false;
			int p;
			int nN = 0;
			for (p = 0; p < reads[readIndex]->length; p++) {
				if (reads[readIndex]->seq[p] == 'N') { nN++; }
			}
			cout << nN << " " << reads[readIndex]->length << endl;
			if (nN == reads[readIndex]->length) {
				nFilt++;
			}
		}
	}
	cout << nFilt << " / " << nHits << endl;
}

