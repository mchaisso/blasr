#ifndef FLAT_MATRIX_H_
#define FLAT_MATRIX_H_

#include <vector>
#include <iostream>
#include <assert.h>
#include "Types.h"
#define rc2index(r, c, rowSize) ( (r) * (rowSize) + c)

using namespace std;
template<typename T>
void CreateFlatMatrix(int rows, int cols, vector<T> &matrix) {
	matrix.resize(rows*cols);
}

template<typename T>
void PrintFlatMatrix(vector<T> &matrix, int rows, int cols, std::ostream &out, int width=6) {
	PrintFlatMatrix((const T*) &matrix[0], rows, cols, out, width);
}

template<typename T>
void PrintFlatMatrix(const T* matrix, int rows, int cols, std::ostream &out, int width=6) { 
	int r, c, i;
	i = 0;
	for (r = 0; r < rows; r++ ) {
		for (c = 0; c < cols; c++ ) {
			out.width(width);
			out << matrix[i] << " ";
			i++;
		}
		out << endl;
	}
}

template<typename T>
class FlatMatrix2D {
 public:
	T* matrix;
	int nRows, nCols;
	int totalSize;
	int RC2Index(int row, int col) {
		return row * nCols + col;
	}

	T* operator[](int row) {
		return &matrix[row*nCols];
	}

	T operator()(int row, int col) {
		int index = RC2Index(row,col);
		return matrix[index];
	}

	FlatMatrix2D() {
		matrix = NULL;
		nRows = nCols = totalSize = 0;
	}
	unsigned int Size() {
		return nRows * nCols;
	}
  
  void Fill(T value) {
    fill(matrix, &matrix[totalSize], value);
  }
  //
  // Pseudonym for grow.
  //
  void Resize(int _nRows, int _nCols) {
    Grow(_nRows, _nCols);
  }

	void Resize(unsigned int totalSize) {
		if (matrix != NULL) {
			delete[] matrix;
		}
		matrix = new T[totalSize];
	}

	FlatMatrix2D(int _nRows, int _nCols) {
		totalSize = 0;
		matrix = NULL;
		//		assert(_nRows > 0);
		//		assert(_nCols > 0);
		if (_nRows == 0 or _nCols == 0) {
			nRows = _nRows; nCols = _nCols;
			return;
		}
		else {
			Grow(_nRows, _nCols);
		}
	}

	void Clear() {
		delete[] matrix;
		matrix = NULL;
		nRows = nCols = 0;
		totalSize = 0;
	}

	void Grow(int _nRows, int _nCols) {
		nRows = _nRows;
		nCols = _nCols;
		if (nRows * nCols > totalSize) {
			if (totalSize != 0)
				delete[] matrix;
			totalSize = nRows * nCols;
			matrix = new T[totalSize];
		}
	}

	T Get(int r, int c) {
		return matrix[r * nCols + c];
	}
	T Set(int r, int c, T v) {
		return matrix[r*nCols+c] = v;
	}
	int Index(int r, int c) {
		return r*nCols + c;
	}
	void Print(std::ostream &out) {
		PrintFlatMatrix(matrix, nRows, nCols, out);
	}

	void Allocate(UInt _nRows, UInt _nCols) {
		nRows = _nRows;
		nCols = _nCols;
		matrix = new T[nRows * nCols];
	}

	void Initialize(T value) {
		T* matPtr, *matEnd;
		matEnd = &matrix[nRows*nCols];
		for (matPtr = &matrix[0]; matPtr != matEnd; ++matPtr) {
			*matPtr = value;
		}
	}
};


template<typename T>
class FlatMatrix3D {
 public:
	T* matrix;
	// 
	// For some reason it makes sense to go from rows,cols to x,y,z 
	// for referencing coordinates.
	//
	int nx, ny, nz;
	int xy;
	int totalSize;
	FlatMatrix3D() {
		nx = 0, ny = 0, nz = 0;
    xy = 0;
		totalSize = 0;
		matrix = NULL;
	}

	FlatMatrix3D(int _nx, int _ny, int _nz) {
		totalSize = 0;
		matrix =NULL;
		assert(_nx > 0);
		assert(_ny > 0);
		assert(_nz > 0);
		Grow(_nx, _ny, _nz);
	}

	void Grow(int _nx, int _ny, int _nz) {
		nx = _nx;
		ny = _ny;
		nz = _nz;
		if (nx * ny * nz > totalSize) {
			if (matrix != NULL) {
				delete[] matrix;
			}
			totalSize = nx*ny*nz;
			matrix = new T[totalSize];
		}
		xy = nx*ny;
	}
	int Index(int x, int y, int z) {
		return z * xy + y * nx + x;
	}
	T Get(int x, int y, int z) {
		return matrix[Index(x,y,z)];
	}
	T Set(int x, int y, int z, T v){ 
		return matrix[Index(x,y,z)] = v;
	}
};
		
	



#endif

#ifdef BLAH

#include <vector>
#include <iostream>
#include <assert.h>
#include "../../Types.h"
#define rc2index(r, c, rowSize) ( (r) * (rowSize) + c)

using namespace std;
template<typename T>
void CreateFlatMatrix(int rows, int cols, vector<T> &matrix) {
	matrix.resize(rows*cols);
}

template<typename T>
void PrintFlatMatrix(vector<T> &matrix, int rows, int cols, std::ostream &out, int width=6) {
	PrintFlatMatrix((const T*) &matrix[0], rows, cols, out, width);
}

template<typename T>
void PrintFlatMatrix(const T* matrix, int rows, int cols, std::ostream &out, int width=6) { 
	int r, c, i;
	i = 0;
	for (r = 0; r < rows; r++ ) {
		for (c = 0; c < cols; c++ ) {
			out.width(width);
			out << matrix[i] << " ";
			i++;
		}
		out << endl;
	}
}

template<typename T>
class FlatMatrix2D {
 public:
	T* matrix;
	int nRows, nCols;
	int totalSize;
	int RC2Index(int row, int col) {
		return row * nCols + col;
	}

	T* operator[](int row) {
		return &matrix[row*nCols];
	}

	T operator()(int row, int col) {
		int index = RC2Index(row,col);
		return matrix[index];
	}

	FlatMatrix2D() {
		matrix = NULL;
		nRows = nCols = totalSize = 0;
	}
	unsigned int Size() {
		return nRows * nCols;
	}
	void Resize(unsigned int totalSize) {
		if (matrix != NULL) {
			delete[] matrix;
		}
		matrix = new T[totalSize];
	}

	FlatMatrix2D(int _nRows, int _nCols) {
		totalSize = 0;
		matrix = NULL;
		//		assert(_nRows > 0);
		//		assert(_nCols > 0);
		if (_nRows == 0 or _nCols == 0) {
			nRows = _nRows; nCols = _nCols;
			return;
		}
		else {
			Grow(_nRows, _nCols);
		}
	}
	void Grow(int _nRows, int _nCols) {
		nRows = _nRows;
		nCols = _nCols;
		if (nRows * nCols > totalSize) {
			if (totalSize != 0)
				delete[] matrix;
			totalSize = nRows * nCols;
			matrix = new T[totalSize];
		}
	}

	T Get(int r, int c) {
		return matrix[r * nCols + c];
	}
	T Set(int r, int c, T v) {
		return matrix[r*nCols+c] = v;
	}
	int Index(int r, int c) {
		return r*nCols + c;
	}
	void Print(std::ostream &out) {
		PrintFlatMatrix(matrix, nRows, nCols, out);
	}

	void Allocate(UInt _nRows, UInt _nCols) {
		nRows = _nRows;
		nCols = _nCols;
		matrix = new T[nRows * nCols];
	}

	void Initialize(T value) {
		T* matPtr, *matEnd;
		matEnd = &matrix[nRows*nCols];
		for (matPtr = &matrix[0]; matPtr != matEnd; ++matPtr) {
			*matPtr = value;
		}
	}

    ~FlatMatrix2D(){
        if (matrix != NULL) {
            delete [] matrix;
            matrix = NULL;
        }
    }
};


template<typename T>
class FlatMatrix3D {
 public:
	T* matrix;
	// 
	// For some reason it makes sense to go from rows,cols to x,y,z 
	// for referencing coordinates.
	//
	int nx, ny, nz;
	int xy;
	int totalSize;
	FlatMatrix3D() {
		nx = 0, ny = 0, nz = 0;
    xy = 0;
		totalSize = 0;
		matrix = NULL;
	}

	FlatMatrix3D(int _nx, int _ny, int _nz) {
		totalSize = 0;
		matrix =NULL;
		assert(_nx > 0);
		assert(_ny > 0);
		assert(_nz > 0);
		Grow(_nx, _ny, _nz);
	}

	void Grow(int _nx, int _ny, int _nz) {
		nx = _nx;
		ny = _ny;
		nz = _nz;
		if (nx * ny * nz > totalSize) {
			if (matrix != NULL) {
				delete[] matrix;
			}
			totalSize = nx*ny*nz;
			matrix = new T[totalSize];
		}
		xy = nx*ny;
	}
	int Index(int x, int y, int z) {
		return z * xy + y * nx + x;
	}
	T Get(int x, int y, int z) {
		return matrix[Index(x,y,z)];
	}
	T Set(int x, int y, int z, T v){ 
		return matrix[Index(x,y,z)] = v;
	}
    ~FlatMatrix3D() {
        if (matrix != NULL) {
            delete [] matrix;
            matrix = NULL;
        }
    }
};
		
	



#endif
