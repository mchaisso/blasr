#ifndef DATA_HDF_HDF_2DARRAY_H_
#define DATA_HDF_HDF_2DARRAY_H_

#include <assert.h>
#include <iostream>
#include "H5Cpp.h"
#include "HDFConfig.h"
#include "HDFData.h"
#include "BufferedHDF2DArray.h"

using namespace std;
using namespace H5;


/*
 *
 * Implementation of a 2-D array for IO from an HDF array.
 * This is templated, but specialized for a few data types, so that 
 * the HDF data types do not need to be specified by somebody when reading.
 *
 * Currently no support exists for reading non-contiguous blocks of data, and
 * the main intended use is to read in increments of rows.

 int main(int argc, char* argv[]) {
	if (argc < 1) {
		cout << "usage: testHDFReading hdfFile" << endl;
		exit(0);
	}

	string hdfFileName = argv[1];
	
	H5File hdfFile;
	hdfFile.openFile(hdfFileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	HDF2DArray<uint16_t> xyArray;
	xyArray.Initialize(hdfFile, "PulseData/BaseCalls/ZMW/HoleXY");
	int curX = 0;
	xyArray.Read(curX, curX + 1, 0, 2, holeXY);

	or, to read a row:
	xyArray.Read(curX, curX+1, holeXY);

 *
 */
template<typename T>
class HDF2DArray : public BufferedHDF2DArray<T> {
 public:
	void WriteRow(const T*data, int dataLength, int destRow=-1) {
		this->writeBuffer = (T*) data;
		this->bufferIndex = dataLength;
		this->bufferSize  = dataLength;
		this->Flush(destRow);
		//
		// Reset status of buffer so that no methods are tricked into
		// thinking this is a valid pointer.
		//
		this->writeBuffer = NULL;
		this->bufferIndex = 0;
		this->bufferSize  = 0;
	}

};

#endif
