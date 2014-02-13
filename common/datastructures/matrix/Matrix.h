#ifndef DATASTRUCTURES_MATRIX_MATRIX_H_
#define DATASTRUCTURES_MATRIX_MATRIX_H_

#include <vector>
#include <iostream>
#include <fstream>
#include "../../Types.h"

using namespace std;
template<typename T>
void CreateMatrix(int rows, int cols, vector<T*> matrix) {
	matrix.resize(rows);
	matrix[0] = new T[rows*cols];
	VectorIndex r = 1;
	for (r = 1; r < rows; r++) {
		matrix[r] = &matrix[cols * r];
	}
}

/*
 *	Implement a matrix as an array into a large allocated buffer of size
 *	nRows * nCols.
 */

template<typename T>
class Matrix {
	VectorIndex nRows;
	VectorIndex nCols;
	T** matrix;
	VectorIndex matrixBufferSize;
	VectorIndex matrixSize;
	VectorIndex rowsBufferSize;
 public:
	Matrix() {
		nRows = 0;
		rowsBufferSize = 0;
		nCols = 0;
		matrixBufferSize = 0;
		matrix = NULL;
	}

	unsigned int size() {
		return nRows * nCols;
	}

  unsigned int GetNCols() {
    return nCols;
  }

  unsigned int GetNRows() {
    return nRows;
  }
	void Resize(VectorIndex nRowsP, VectorIndex nColsP) {

		nRows = nRowsP;
		nCols = nColsP;
		matrixSize = nRows * nCols;
		if (nRows * nCols > matrixBufferSize) {
			matrixBufferSize = nRows * nCols;
			if (nRows > rowsBufferSize) {
				if (matrix != NULL) { 
					delete[] matrix; matrix = NULL; 
				}
			}
			if (matrix == NULL) {
				matrix = new T*[nRows];
			}
			else {
				if (matrix[0] != NULL) {
					delete[] matrix[0];
				}
				if (nRows > rowsBufferSize) {
					delete[] matrix;
				}
			}
			matrix[0] = new T[matrixBufferSize];
			VectorIndex rowIndex;
			for (rowIndex = 1; rowIndex < nRows; rowIndex++ ){
				matrix[rowIndex] = &matrix[0][nCols * rowIndex];
			}
		}
	}

	Matrix(VectorIndex nRowsP, VectorIndex nColsP) {
		Resize(nRowsP, nColsP);
	}

	void Reference(Matrix<T> &rhs) {
		matrix = rhs.matrix;
	}
	void Initialize(T value) {
		std::fill(&matrix[0][0], &matrix[0][matrixSize], value);
	}

	T* operator[](VectorIndex rowIndex) {
		assert(rowIndex < nRows);
		return matrix[rowIndex];
	}

	void Free() {
		if (matrix != NULL) {
			if (matrix[0] != NULL) {
				delete[] matrix[0];
			}
			delete[] matrix;
		}
	}
	
	void Print(ofstream &out) {
		VectorIndex i;
		VectorIndex j;
		for (i = 0; i < nRows; i++) {
			for (j = 0; j < nCols; j++ ){ 
				out << matrix[i][j] << " ";
			}
			out << endl;
		}
	}
};
		
#endif
