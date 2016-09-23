#ifndef SIMULATOR_OUTPUT_SAMPLE_LIST_SET_H_
#define SIMULATOR_OUTPUT_SAMPLE_LIST_SET_H_

#include <map>
#include <string>
#include <iostream>
#include "OutputSampleList.h"


using namespace std;

typedef map<string, OutputSampleList> OutputSampleListMap;
class OutputSampleListSet {
 public:
  OutputSampleListMap listMap;
  vector<string> keys;
  int keyLength;
  int nSufficient;
  int sampleSpaceSize;
  int keySize;
	int minSamples;
	int maxSamples;
  vector<int> lengths;
  OutputSampleListSet(int keySizeP) {
    minSamples = 500;
		maxSamples = 2000;
    nSufficient = 0;
    keySize = keySizeP;
    sampleSpaceSize = 1 << (2*keySize);
  }

  void Write(ofstream &out) {
    // Say how many elements to write.
    OutputSampleListMap::iterator mapIt;
    int setSize = listMap.size();
    out.write((char*) &setSize, sizeof(int));
    int keySize = 0;
    // Say how large each element is.
    if (listMap.size() > 0) {
      keySize = listMap.begin()->first.size();
    }
    out.write((char*)&keySize, sizeof(int));

    for (mapIt = listMap.begin(); mapIt != listMap.end(); ++mapIt) {
      string mapItKey = mapIt->first;
      out.write((char*) mapItKey.c_str(), sizeof(char) * mapItKey.size());
      mapIt->second.Write(out);
    }
    int numLengths = lengths.size();
    out.write((char*) &numLengths, sizeof(int));
    int i;
    for ( i = 0; i < lengths.size(); i++) {
      out.write((char*) &lengths[i], sizeof(int));
    }
  }

  void Read(string &inName) {
    ifstream in;
    CrucialOpen(inName, in, std::ios::in|std::ios::binary);
    Read(in);
    in.close();
  }

  void Read(ifstream &in) {
    int setSize;
    in.read((char*) &setSize, sizeof(int));
    in.read((char*) &keyLength, sizeof(int));

    if (keyLength == 0 or setSize == 0) { return; }
    char *key = new char[keyLength+1];
    key[keyLength] = '\0';
    int i;
    for (i = 0; i < setSize; i++) {
      in.read(key, sizeof(char)*keyLength);
      listMap[key].Read(in);
    }
    int numLengths;
    in.read((char*) &numLengths, sizeof(int));
    if (numLengths > 0) {
      lengths.resize(numLengths);
    }
    for (i = 0; i < numLengths; i++) {
      in.read((char*) &lengths[i], sizeof(int));
    }
  }
  
  void AppendOutputSample(string key, OutputSample &sample) {

    if (listMap[key].size() < minSamples) {
			if (listMap[key].size() < maxSamples) {
				listMap[key].push_back(sample);
			}
      if (listMap[key].size() == minSamples) {
        nSufficient++;
        cout << nSufficient << " / " << sampleSpaceSize << endl;
				if (sampleSpaceSize - nSufficient == 10) {
					int i;
					const char nucs[] = {'A','C','G','T'};
					for (i = 0; i < sampleSpaceSize; i++) {
						string s;
						int j;
						int v = i;
						for (j = 0; j < 5; j++) {
							s.push_back(nucs[v & 3]);
							v >>=2;
						}
						if (listMap.find(s) == listMap.end()) {
							cout << "Last 10: " << s << endl;
						}
					}
				}
      }
    }
  }

  bool Sufficient() {
    return nSufficient + 1 >= sampleSpaceSize;
  }

  void SampleRandomSample(string key, OutputSample &sample) {
    if (listMap.find(key) == listMap.end()) {
      cout << listMap.size() << endl;
      cout <<"ERROR, " << key << " is not a sampled context." << endl;
      int i;
      for (i = 0; i < key.size(); i++) {
        char c = toupper(key[i]);
        if (c != 'A' and c != 'C' and c != 'G' and c != 'T') {
          cout << "The nucleotide " << c << " is not supported." << endl;
        }
      }
      exit(1);
    }
    sample = listMap[key][RandomInt(listMap[key].size())];
  }


};

#endif
