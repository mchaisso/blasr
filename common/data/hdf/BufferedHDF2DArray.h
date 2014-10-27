#ifndef DATA_BUFFERED_HDF_HDF_2DARRAY_H_
#define DATA_BUFFERED_HDF_HDF_2DARRAY_H_

#include "H5Cpp.h"
#include "HDFConfig.h"
#include "HDFData.h"
#include "HDFGroup.h"
#include "HDFWriteBuffer.h"
#include <assert.h>
#include <iostream>

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
	BufferedHDF2DArray<uint16_t> xyArray;
	xyArray.Initialize(hdfFile, "PulseData/BaseCalls/ZMW/HoleXY");
	int curX = 0;
	xyArray.Read(curX, curX + 1, 0, 2, holeXY);

	or, to read a row:
	xyArray.Read(curX, curX+1, holeXY);

 *
 */
template<typename T>
class BufferedHDF2DArray : public HDFData, public HDFWriteBuffer<T> {
	hsize_t   nDims;
	hsize_t   *dimSize;
	int       maxDims;
	int       rowLength, colLength;
 public:
	unsigned int GetNRows() {
		return rowLength;
	}
	unsigned int GetNCols() {
		return colLength;
	}

 BufferedHDF2DArray(CommonFG *_container, string _datasetName) : HDFData(_container, _datasetName) {
	}

 BufferedHDF2DArray() : HDFData() {
    maxDims = 0;
		nDims = 2;
		dimSize =NULL;
		rowLength = -1;
		colLength = -1;
	}
  
  void Close() {

		//
		// Clean up the write buffer.
		//
		//		Flush();
		if (dimSize != NULL) {
			delete[] dimSize;
      dimSize = NULL;
		}
		this->HDFWriteBuffer<T>::Free();
  }

	~BufferedHDF2DArray() {
    Close();
	}
	/*
	 * Initialize HDF2D for reading.  No write buffer initialization is
	 * required.  The assumption is that the dataspace is in two
	 * dimensions, and this exits without grace if it is not. 
	 */
	int Initialize(HDFGroup &group, string datasetName, int _rowLength=0, int _bufferSize=0, bool createIfMissing=true) {
    bool groupContainsDataset = group.ContainsObject(datasetName);
    if (groupContainsDataset == false) {
      //
      // Do some error checking.
      //
      if (createIfMissing == false) {
        cout << "ERROR! Could not open dataset " << datasetName << endl;
        exit(1);
      }
      if (_rowLength == 0) {
        cout << "ERROR!  Improper usage of BufferedHDF2DArray::Initialize.  The 2D Array "<<endl
             << "is being created but is given a number of columns of 0." << endl;
        exit(1);
      }
      Create(&group.group, datasetName, _rowLength);
    }
    else {
      InitializeDataset(group.group, datasetName);
      try {
        dataspace = dataset.getSpace();
      }
      catch(H5::DataSetIException &e) { 
        cout << e.getDetailMsg() << endl;
        exit(1);
      }

      maxDims   = MAX_DIMS;
      try {
        nDims     = dataspace.getSimpleExtentNdims();
        /*
         * Prevent abuse of this class for multidimensional IO.
         */
        if (nDims != 2) {
          cout << "ERROR in HDF format: dataset: " << datasetName << " should be 1-D, but it is not." << endl;
          exit(1);
        }

        /*
         * Load in the size of this dataset, and make a map to the whole thing.
         */
        dimSize = new hsize_t[nDims];
        dataspace.getSimpleExtentDims(dimSize);
        rowLength = dimSize[0];
        colLength = dimSize[1];
        if (rowLength == 0) {
            dataspace.close();
            // DONT create a real dataspace if the size is 0.
            // cout << "WARNING, trying to create a zero sized dataspace." << endl;
            return 1;
        }
        fullSourceSpace = DataSpace(2, dimSize);
        dataspace.close();
      }
      catch(Exception &e) {
        cout << e.getDetailMsg() << endl;
        exit(1);
      }
    }
      return 1;
	}

	int size() {
        // Why assert nDims == 1 for 2D Array?
		assert(nDims == 1);
		dataspace.getSimpleExtentDims(dimSize);
		return dimSize[0];
	}

	/*
	 * Read rows in the range (startX, endX] in to dest.
	 */

	void Read(int startX, int endX, DataType typeID, T*dest) {
		Read(startX, endX, 0, dimSize[1], typeID, dest);
	}
	
	void Read(int startX, int endX, T*dest) {
		Read(startX, endX, 0, dimSize[1], dest);
	}
	/*
	 * This is the non-specialized definition.  Since this should only
	 * operate on specialized types, report an error and bail.
	 */
	void Read(int startX, int endX, int startY, int endY, T* dest) {
		assert("ERROR, calling Read with an unsupported type. Use Read(startx,endx, starty,endy,datatype, dest) instead." == 0);
		exit(1);
	}

	void Read(int startX, int endX, int startY, int endY, DataType typeID, T *dest) {
		hsize_t memSpaceSize[2] = {0, 0};
		memSpaceSize[0] = endX - startX;
		memSpaceSize[1] = endY - startY;
		hsize_t sourceSpaceOffset[2] = {0, 0};
		sourceSpaceOffset[0] = startX;
		sourceSpaceOffset[1] = startY;
		
		DataSpace destSpace(2, memSpaceSize);		
		fullSourceSpace.selectHyperslab(H5S_SELECT_SET, memSpaceSize, sourceSpaceOffset);
		dataset.read(dest, typeID, destSpace, fullSourceSpace);
		destSpace.close();
	}
	
	void Create(CommonFG *_container, string _datasetName, int _rowLength) {
    container   = _container;
    datasetName = _datasetName;
    rowLength   = _rowLength;
		//
		// Make life easy if the buffer is too small to fit a row --
		// resize it so that rows may be copied and written out in an
		// atomic unit.
		//
		if (this->bufferSize < rowLength) {
			// When the buffer size is greater than 0, the write buffer
			// should exist.
      if (this->bufferSize > 0) {
        assert(this->writeBuffer != NULL);
        delete[] this->writeBuffer;
      }
			this->writeBuffer = new T[rowLength];
			this->bufferSize = rowLength;
		}

		hsize_t dataSize[2]    = {0, rowLength};
		hsize_t maxDataSize[2] = {H5S_UNLIMITED, rowLength};
		DataSpace fileSpace(2, dataSize, maxDataSize);
		DSetCreatPropList cparms;

		/*
		 * For some reason, chunking must be enabled when creating a dataset
		 * that  has an unlimited dimension.  Of course, this is not
		 * mentioned in the hdf5 c++ documentation, because that
		 * docuemntation was written for people who enjoy learning how to
		 * use an API by reading comments in source code.
		 */
		hsize_t chunkDims[2] = {16384, rowLength};
		cparms.setChunk( 2, chunkDims );
		TypedCreate(fileSpace, cparms);
    fileSpace.close();
    
    //
    // Set some flags that indicate this dataset is ready for writing.
    //
    fileDataSpaceInitialized = true;
    isInitialized = true;
	}

	void TypedCreate(DataSpace &fileSpace, DSetCreatPropList &cparms) {
		assert("Error, calling HDF2DArray<T>::TypedCreate on an unsupported type.  A specialization must be written in HDF2DArray.h" == 0);
	}
	// Append
	void TypedWriteRow(const T*, const DataSpace &memoryDataSpace, const DataSpace &fileDataSpace) {
		assert("Error, calling HDF2DArray<T>::TypedWriteRow on an unsupported type.  A specialization must be written in HDF2DArray.h" == 0);
	}


	/*
	 * This code is copied directly form BufferedHDFArray.  I'm not sure
	 * how to set up the objects nicely to share the code between the
	 * two since the Flush() function is different.  There probably is a
	 * design pattern or simply better way to engineer this, but for now
	 * it's 15 lines of code.
	 */
	 
	void WriteRow(const T *data, int dataLength, int destRow=-1) {
		// Fill the buffer with data. When there is overflow, write
		// that out to disk.
		//
		int dataIndex = 0;
		int bufferCapacity;
		int bufferFillSize = 0;
		bool flushBuffer;
		while(dataIndex < dataLength) {
			//
			// Compute the capacity of this buffer to fit an integral number
			// of rows into it.
			//
			bufferCapacity = (this->bufferSize / rowLength)*rowLength - this->bufferIndex;
			flushBuffer = false;
			if (bufferCapacity  > dataLength - dataIndex) {
				bufferFillSize = dataLength - dataIndex;
			}
			else {
				bufferFillSize = bufferCapacity;
				flushBuffer = true;
			}
			memcpy((void*) &this->writeBuffer[this->bufferIndex], (void*) &data[dataIndex], sizeof(T)*bufferFillSize);
			dataIndex   += bufferFillSize;
			this->bufferIndex += bufferFillSize;
			if (flushBuffer) {
				Flush(destRow);
			}
      //
      //  When not appending, increment the position of where the data
      //  is to be written.
      //
      if (destRow != -1) {
        destRow += this->bufferIndex / rowLength;
      }
		}
	}

	void Flush(int destRow = -1) {

    //
    // A default writeRow of -1 implies append
    //
    int numRowsToCreate, numDataRows;
    //
    // this->bufferIndex points after the end of the last data in the
    // buffer (full rows), so this->bufferIndex / rowLength is the
    // number of number of rows to create.
    //
    numDataRows = this->bufferIndex / rowLength;

    if (destRow < 0) {
      numRowsToCreate = this->bufferIndex / rowLength;  
    }
    else {
      numRowsToCreate = this->bufferIndex / rowLength + destRow;
    }
    if (numDataRows > 0) {
      assert(fileDataSpaceInitialized);
      
      DataSpace fileSpace;
      fileSpace = dataset.getSpace();
			
      //
      // Load the current size of the array on disk.
      //
      hsize_t fileArraySize[2], fileArrayMaxSize[2], blockStart[2];
      fileSpace.getSimpleExtentDims(fileArraySize, fileArrayMaxSize);

      // Save this for later to determine the offsets
      blockStart[0] = fileArraySize[0];
      blockStart[1] = fileArraySize[1];

      //
      // Calculate the number of rows to create.  This is dependent
      // on the current file size, the destination of where the data
      // will go, and how much to write.
      //

      if (destRow == -1) {
        fileArraySize[0] += numDataRows;
      }
      else {
        // If the data cannot fit in the current file size, extend
        // it,  otherwise, do not toch the file array size.
        if (destRow + numDataRows > fileArraySize[0]) {
          fileArraySize[0] = destRow + numDataRows;
        }
      }

      //
      // Make room in the file for the array.
      //
      dataset.extend(fileArraySize);
		
      DataSpace extendedSpace = dataset.getSpace();
      //
      // Store the newly dimensioned dataspaces.
      //
      fileSpace.getSimpleExtentDims(fileArraySize, fileArrayMaxSize);			
      int extendedSize = extendedSpace.getSimpleExtentNpoints();
      //
      // Configure the proper addressing to append to the array.
      //
      hsize_t dataSize[2];
      dataSize[0] = numDataRows;
      dataSize[1] = rowLength;
      hsize_t offset[2];
      //
      // Determine which row to write to.
      //
      if (destRow == -1) {
        offset[0] = blockStart[0];
      }
      else {
        offset[0] = destRow;
      }
      offset[1] = 0;
      extendedSpace.selectHyperslab(H5S_SELECT_SET, dataSize, offset);
      DataSpace memorySpace(2, dataSize);

      //
      // Finally, write out the data.  
      // This uses a generic function which is specialized with
      // templates later on to t
      // memorySpace addresses the entire array in linear format
      // fileSpace addresses the last dataLength blocks of dataset.
      //
      TypedWriteRow(this->writeBuffer, memorySpace, extendedSpace);
      memorySpace.close();
      extendedSpace.close();
      fileSpace.close();
    }
		this->ResetWriteBuffer();
	}
};

#define DEFINE_TYPED_WRITE_ROW(T, Pred) template<>\
void BufferedHDF2DArray<T>::TypedWriteRow(const T *data, const DataSpace &memorySpace, const DataSpace &fileSpace) {\
	dataset.write(data, Pred, memorySpace, fileSpace);\
}

DEFINE_TYPED_WRITE_ROW(int, PredType::NATIVE_INT);
DEFINE_TYPED_WRITE_ROW(unsigned int, PredType::NATIVE_UINT);
DEFINE_TYPED_WRITE_ROW(unsigned char, PredType::NATIVE_UINT8);
DEFINE_TYPED_WRITE_ROW(uint16_t, PredType::NATIVE_UINT16);
DEFINE_TYPED_WRITE_ROW(int16_t, PredType::NATIVE_INT16);


#define DEFINE_READ(T, Pred) template<>\
void BufferedHDF2DArray<T>::Read(int startX, int endX, int startY, int endY, T* dest) {\
	Read(startX, endX, startY, endY, Pred, dest);\
}


DEFINE_READ(int, PredType::NATIVE_INT);
DEFINE_READ(unsigned int, PredType::NATIVE_UINT);
DEFINE_READ(char, PredType::NATIVE_INT8);
DEFINE_READ(unsigned char, PredType::NATIVE_UINT8);
DEFINE_READ(uint16_t, PredType::NATIVE_UINT16);
DEFINE_READ(int16_t, PredType::NATIVE_INT16);

#define DEFINE_TYPED_CREATE(T, Pred)template<>\
void BufferedHDF2DArray<T>::TypedCreate(DataSpace &fileSpace, DSetCreatPropList &cparms) {\
	dataset = container->createDataSet(datasetName.c_str(), Pred, fileSpace, cparms);\
}

DEFINE_TYPED_CREATE(int, PredType::NATIVE_INT)
DEFINE_TYPED_CREATE(unsigned int, PredType::NATIVE_UINT)
DEFINE_TYPED_CREATE(char, PredType::NATIVE_INT8)
DEFINE_TYPED_CREATE(unsigned char, PredType::NATIVE_UINT8)
DEFINE_TYPED_CREATE(uint16_t, PredType::NATIVE_UINT16)
DEFINE_TYPED_CREATE(int16_t, PredType::NATIVE_INT16)


#endif
