#ifndef DATA_HDF_HDF_ATOM_H_
#define DATA_HDF_HDF_ATOM_H_

#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include <stdint.h>
#include "H5Cpp.h"

#include "HDFConfig.h"
#include "HDFGroup.h"
#include "HDFData.h"

using namespace std;
using namespace H5;

template<typename T>
class HDFAtom : public HDFData {
 public:
	Attribute attribute;

	bool initialized;
	HDFAtom() {
		initialized = false;
	}
	~HDFAtom() {
		if (initialized) {
			attribute.close();
		}
	}

  H5Object *GetObject() {
    return NULL;
  }

  int Initialize(H5Object &object, string attributeName, bool createIfMissing=false) {
		attribute = object.openAttribute(attributeName.c_str());
		initialized = true;
		return 1;
	}
	
	int Initialize(HDFGroup &group, string attributeName, bool createIfMissing=false) {
		return Initialize(group.group, attributeName);
	}
	
	int Initialize(HDFData &data, string attributeName, bool createIfMissing=false) {
		return Initialize(data.dataset, attributeName);
	}

  int Initialize(Group &object, string attributeName, bool createIfMissing=false) {
    try {
      attribute = object.openAttribute(attributeName.c_str());
    }
    catch (Exception e) {
      cout << "ERROR. Could not open attribute " << attributeName << endl;
      exit(1);
    }
		initialized  = true;
		return 1;
	}

	int Initialize(H5File &hdfFile, string groupName, string attributeName) {
		HDFGroup group;
		group.Initialize(hdfFile, groupName);
		attribute = group.group.openAttribute(attributeName.c_str());
		initialized = true;
		return 1;
	}

	//
	// This handles creation of all non-string types.  A specialization
	// for strings is provided below.
	//
	void Create(H5Object &object, string atomName) {
		hsize_t defaultDims[] = {1};
		DataSpace defaultDataSpace(1, defaultDims);
		TypedCreate(object, atomName, defaultDataSpace);
	}
	

	void Create(H5Object &object, string name, string value) {
		StrType strType(0, value.size());
		hsize_t defaultDims[] = {};
		//		DataSpace defaultDataSpace(1, defaultDims);
		attribute = object.createAttribute(name.c_str(), strType, DataSpace(0,NULL));
		initialized = true;
		attribute.write(strType, value.c_str());
	}
	void Write(string &value) {
		StrType strType(0, value.size());
		hsize_t defaultDims[] = {};
		//		DataSpace defaultDataSpace(1, defaultDims);
		attribute.write(strType, value.c_str());
	}

  void Create(H5Object &object, string name, vector<int> &vect) {
    hsize_t length = vect.size();
    ArrayType arrayDataType(PredType::NATIVE_INT, 1, &length);
    hsize_t one = 1;
    attribute = object.createAttribute(name.c_str(), PredType::NATIVE_INT, DataSpace(1, &length));
    attribute.write(PredType::NATIVE_INT, &((vect)[0]));    
  }
	

	
	
  void Create(H5Object &object, string name, vector<string> &vect) {
    hsize_t length = vect.size();
    StrType strType(0,H5T_VARIABLE);
    ArrayType arrayDataType(strType, 1, &length);
    hsize_t one = 1;
    attribute = object.createAttribute(name.c_str(), strType, DataSpace(1, &length));
    attribute.write(strType, &((vect)[0]));    
  }

	void TypedCreate(H5Object &object, string &atomName, DataSpace &dataSpace) {
		assert("Calling HDFAtom<T>::typedCreate on an unsupported type" == 0);
	}
	
	void Write(T value) {
		assert("Calling HDFAtom<T>::Write on an unsupported type" == 0);
	}

	void Read(T& value) {
		assert("Calling read on an unsupported type!" == 0);
	}

};

//
// Special create for strings.  Since this uses a StrType for the
// typename rather than specifying a PredType, it mertis its own
// function.
//


template<>
void HDFAtom<string>::Create(H5Object &object, string atomName) {
	StrType strType(0, H5T_VARIABLE);
	hsize_t defaultDims[] = {1};
	DataSpace defaultDataSpace(1, defaultDims);
  //	attribute = object.createAttribute(atomName.c_str(), strType, defaultDataSpace);
	attribute = object.createAttribute(atomName.c_str(), strType, DataSpace(H5S_SCALAR));
	initialized= true;
}



#define MAKE_TYPED_CREATE(T, predtype) template<> \
	void HDFAtom<T>::TypedCreate(H5Object &object, string &atomName, DataSpace &defaultDataSpace) {				\
  attribute = object.createAttribute(atomName.c_str(), (predtype), defaultDataSpace );	\
}


MAKE_TYPED_CREATE(int, PredType::NATIVE_INT)
MAKE_TYPED_CREATE(unsigned int, PredType::NATIVE_UINT)
MAKE_TYPED_CREATE(unsigned char, PredType::NATIVE_UINT8)
MAKE_TYPED_CREATE(char, PredType::NATIVE_INT8)
MAKE_TYPED_CREATE(float, PredType::NATIVE_FLOAT)
MAKE_TYPED_CREATE(uint64_t, PredType::STD_I64LE)


template<>
void HDFAtom<vector<int> >::Write(const vector<int> vect) {
  hsize_t  length = vect.size();
  DataType baseType = PredType::NATIVE_INT;
  ArrayType arrayDataType(baseType, 1, &length);
  attribute.write(arrayDataType, &((vect)[0]));
}


template<>
void HDFAtom<string>::Write(string value) {
  //	StrType strType(0, value.size());
	StrType strType(0, H5T_VARIABLE);
	attribute.write(strType, H5std_string(value.c_str()));
}

template<>
void HDFAtom<uint64_t>::Write(uint64_t value) {
	attribute.write( PredType::STD_I64LE, &value);
}

template<>
void HDFAtom<int>::Write(int value) {
	attribute.write( PredType::NATIVE_INT, &value);
}

template<>
void HDFAtom<unsigned int>::Write(unsigned int value) {
	attribute.write( PredType::NATIVE_INT, &value);
}

template<>
void HDFAtom<unsigned char>::Write(unsigned char value) {
	attribute.write( PredType::NATIVE_UINT8, &value);
}

template<>
void HDFAtom<char>::Write(char value) {
	attribute.write( PredType::NATIVE_INT8, &value);
}

template<>
void HDFAtom<float>::Write(float value) {
	attribute.write( PredType::NATIVE_FLOAT, &value);
}
/*
template<>
void HDFAtom<vector<string> >::Write(vector<string> &values) {
  vector<char*> stringPtrs;
  int i;
  for (i = 0;i < values.size(); i++ ){
    stringPtrs[i] = values[i].c_str();
  }
	// Declare and initialize vector of pointers to string attribute list.
	//	if (nPoints > 1) {
	// Copy the pointers.

	// Copy the strings into memory the main program has control over.
	unsigned int i;
	for (i = 0; i < ptrsToHDFControlledMemory.size(); i++ ){
		values.push_back(ptrsToHDFControlledMemory[i]);
		free(ptrsToHDFControlledMemory[i]);
	}
}
*/

template<>
void HDFAtom<string>::Read(string &value) {
	/*
	 * Read in a string that has been stored either as an array or a
	 * variable length string.  To decide which, query the
	 * isVariableStr() option.
	 */
	StrType stringType = attribute.getStrType();
	bool stringIsVariableLength = stringType.isVariableStr();
	if (stringIsVariableLength) 
		attribute.read(stringType, value);
	else {
		hsize_t stsize = attribute.getStorageSize();
		value.resize(stsize);
		//		char *valueStr = new char[stsize+1];
		attribute.read(stringType, &value[0]);
		if (stsize > 0 and value[stsize-1] == '\0') {
			value.resize(stsize-1);
		}
		//		valueStr[stsize] = '\0';
		//		value = valueStr;
		// This read an extra '\0', which is handled by the string class
		//		if (stsize > 0) {
			//			value = valueStr;
		//			delete[] valueStr;
		//		}
	}
}

template<>
void HDFAtom<int>::Read(int &value) {
	DataType intType(PredType::NATIVE_INT);
	attribute.read(intType, &value);
}

template<>
void HDFAtom<uint64_t>::Read(uint64_t &value) {
	DataType intType(PredType::STD_I64LE);
	attribute.read(intType, &value);
}

template<>
void HDFAtom<unsigned int>::Read(unsigned int &value) {
	DataType uintType(PredType::NATIVE_UINT);
	attribute.read(uintType, &value);
}

template<>
void HDFAtom<float>::Read(float &value) {
	DataType type(PredType::NATIVE_FLOAT);
	attribute.read(type, &value);
}

template<>
void HDFAtom<vector<string> >::Read(vector<string> &values) {
	string value;

	/*
	 * This attribute is an array of strings. They are read in by
	 * storing pointers to strings in memory controlled by HDF.  To read
	 * the strings, read the pointers into a temporary array, then copy
	 * those strings to the values array. This way when the values array
	 * is destroyed, it will not try and get rid of space that is under
	 * HDF control.
	 */
  DataSpace attributeSpace = attribute.getSpace();
	hsize_t nPoints;
	nPoints = attributeSpace.getSelectNpoints();
	DataType attrType = attribute.getDataType(); // necessary for attr.read()
	hsize_t stsize = attribute.getStorageSize();

	// Declare and initialize vector of pointers to string attribute list.
	//	if (nPoints > 1) {
	vector<char*> ptrsToHDFControlledMemory;
	ptrsToHDFControlledMemory.resize(nPoints);
	// Copy the pointers.
	attribute.read(attrType, &ptrsToHDFControlledMemory[0]);
	// Copy the strings into memory the main program has control over.
	unsigned int i;
	for (i = 0; i < ptrsToHDFControlledMemory.size(); i++ ){
		values.push_back(ptrsToHDFControlledMemory[i]);
		free(ptrsToHDFControlledMemory[i]);
	}
}


#endif
