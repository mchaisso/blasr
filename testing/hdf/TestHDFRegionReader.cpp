#include "data/hdf/HDFRegionTableReader.h"
#include "datastructures/reads/RegionTable.h"
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
	if (argc < 2) {
		cout << "usage: testHDFRegionReader hdfFile" << endl;
		exit(0);
	}
	string inFileName = argv[1];

	HDFRegionTableReader reader;
	reader.Initialize(inFileName);

	RegionTable table;
	reader.ReadTableAttributes(table);
	RegionAnnotation annotation;
	int rowNumber = 0;
	while(reader.GetNext(annotation)) {
		int i;
		cout << rowNumber << ": ";
		for (i = 0;i<5;i++) {cout << annotation.row[i] << " ";};
		cout <<endl;
		++rowNumber;
	}
	return 0;
}
	
