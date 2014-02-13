#include "data/hdf/HDFArray.h"
#include "data/hdf/HDF2DArray.h"
#include "H5Cpp.h"
#include <string>
#include <iostream>
#include <stdlib.h>
#include <assert.h>

using namespace std;
using namespace H5;

/*
 * Agglomerate tester for hdf reading.  Right now this only reads arrays and a 
 * 2D array.
 */

int main(int argc, char* argv[]) {
	if (argc < 1) {
		cout << "usage: testHDFReading hdfFile" << endl;
		exit(0);
	}

	string hdfFileName = argv[1];
	
	HDFFile hdfFile;
	hdfFile.Open(hdfFileName, H5F_ACC_RDONLY);

	
	HDFArray<int> nElemArray;
	HDFArray<char> baseArray;
	HDFArray<unsigned char> qualArray;
	HDF2DArray<uint16_t> xyArray;

	nElemArray.Initialize(hdfFile.rootGroup, "PulseData/BaseCalls/ZMW/NumEvent");
	baseArray.Initialize(hdfFile.rootGroup, "PulseData/BaseCalls/Basecall");
	qualArray.Initialize(hdfFile.rootGroup, "PulseData/BaseCalls/QualityValue");
	xyArray.Initialize(hdfFile.rootGroup, "PulseData/BaseCalls/ZMW/HoleXY", 2);
	int size = nElemArray.size();

	cout << "got size: " << size << endl;


	int i;
	int cur = 0;
	int curXY = 0;
	for (i = 0; i < size; i++ ){
		int nElem;
		nElemArray.Read(i, i+1, &nElem);
		if (nElem > 0) {
		char *seq;
		seq = new char[nElem+1];
		seq[nElem] = '\0';
		baseArray.Read(cur, cur + nElem, seq);
		unsigned char *qual = new unsigned char[nElem];
		qualArray.Read(cur, cur + nElem, qual);

		cout << seq << endl;
		int j;
		for (j = 0; j < nElem; j++ ){
			cout << (int) qual[j] << " ";
		}
		cout << endl;
		}

		uint16_t *holeXY = new uint16_t[2];
		xyArray.Read(curXY, curXY + 1, holeXY);
		curXY += 1;
		cout <<" hole xy: " << holeXY[0] << " " << holeXY[1] << endl;
		cur+= nElem;
	}

	return 0;
}
