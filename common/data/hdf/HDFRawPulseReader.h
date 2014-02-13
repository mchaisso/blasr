#ifndef DATA_HDF_HDF_RAW_PULSE_READER_H_
#define DATA_HDF_HDF_RAW_PULSE_READER_H_

//
// Load pulse information from the 
//DEFAULT_PULSE_FIELDS = [ "PulseWidth", "pkmid", "StartTime", "IPD", "Light", "ClassifierQV", "QualityValue" ]
//                                     ,        ,  u32       ,      
//BASEMAP = numpy.array(['-','A','C','-','G','-','-','-','T', '-','-','-','-','-','-','-','N'])


#include "data/hdf/HDFArray.h"
#include "data/hdf/HDFFile.h"

class HDFRawPulseReader {
 public:
	vector<HDFArray*> fields;
	vector<string> fieldNames;
	H5File hdfRawPulseFile;
	Group  rootGroup;
	Initialize(string fileName, const H5::FileAccPropList & fileAccPropList=H5::FileAccPropList::DEFAULT) {
		/*
		 * Open the file for reading.
		 */
		try {
			hdfRawPulseFile.openFile(fileName.c_str(), H5F_ACC_RDONLY, fileAccPropList);
		}
		catch (Exception &e) {
			cout << e.getDetailMsg() << endl;
			return 0;
		}
		rootGroup = hdfRawPulseFile.openGroup("/");
		
		fieldNames.push_back("PulseWidth");
		fieldNames.push_back("pkmid");
		fieldNames.push_back("StartTime");
		fieldNames.push_back("IPD");
		fieldNames.push_back("Light");
		fieldNames.push_back("ClassifierQV");
		fieldNames.push_back("QualityValue");
		
		int filedIndex;
		for (fieldIndex = 0; fieldIndex < fieldNames.size(); fieldIndex++) {
			fields.push_back(new HDFArray);
		}
	}
		


};


#endif


