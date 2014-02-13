#ifndef DATA_HDF_HDF_WRITE_BUFFER_H_
#define DATA_HDF_HDF_WRITE_BUFFER_H_


template<typename T>
class HDFWriteBuffer {
	public:
	T         *writeBuffer;
	int       bufferIndex;
	int       bufferSize;

	HDFWriteBuffer() {
		writeBuffer = NULL;
		bufferIndex = 0;
		bufferSize  = 0;
	}

	void InitializeBuffer(int pBufferSize) {
		bufferSize = pBufferSize;
		if (bufferSize > 0) {
			writeBuffer = new T[bufferSize];
		}
		else {
			writeBuffer = NULL;
		}
	}		

  void Free() {
		if (writeBuffer) {
			delete[] writeBuffer;
			writeBuffer = NULL;
		}
  }

	~HDFWriteBuffer() {
    Free();
	}

	void ResetWriteBuffer() {
		bufferIndex = 0;
	}

	bool WriteBufferEmpty() {
		return (bufferIndex == 0);
	}

};


#endif
