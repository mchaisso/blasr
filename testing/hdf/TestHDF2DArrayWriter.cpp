#include "../..//common/data/hdf/HDFFile.h"
#include "../../common/data/hdf/HDF2DArray.h"


int main() {
	HDFFile hdfFile;
  hdfFile.Open("test_write_2d_int_array.hdf");

  HDF2DArray<int> intArray;
  intArray.Initialize(hdfFile.rootGroup, "IntData", 5);
  int value[15] = {1,2,3,4,5, 11,12,13,14,15, 21,22,23,24,25};
  //  intArray.WriteRow(value, 15);
  // This should appear appended.
  intArray.WriteRow(value, 15, 8);
  intArray.WriteRow(value, 15, 1);


  HDFFile appendedFile;
  string name2("test_2d_append.hdf");
  /*  appendedFile.Open(name2);
  //  HDF2DArray<int> intArray;
  intArray.Initialize(appendedFile.rootGroup, "IntData", 5);
  intArray.WriteRow(value, 15);    
  intArray.WriteRow(value, 15);    
  intArray.WriteRow(value, 15);    
  intArray.WriteRow(value, 15);    
  appendedFile.Close();
  */
  HDFFile ap2;
  ap2.Open(name2);
  intArray.Initialize(ap2.rootGroup, "IntData");
  intArray.Read(0, 1, value);

};
