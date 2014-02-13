#include "data/hdf/HDFAtom.h"
#include "data/hdf/HDFArray.h"
#include <string>

using namespace std;

int main(int argc, char* argv[]) {

	string hdfFileName   = argv[1];
  /*	string datasetName   = argv[2];
	string attributeName = argv[3];
	*/
  //	hdfFile.openFile(hdfFileName.c_str(),  H5F_ACC_RDWR);
  FileCreatPropList filePropList;
  hsize_t ub = filePropList.getUserblock();
  filePropList.setUserblock(512);
  HDFFile hdfFile;
  hdfFile.Open(hdfFileName, H5F_ACC_TRUNC);
  
	HDFAtom<int> intAtom;
	HDFAtom<string> strAtom;
  //  HDFAtom<vector<string> > strVectAtom;
  
  //	strAtom.Initialize(hdfFile, datasetName, attributeName);
  
  //
  // Make a fake dataset to add attributes to.
  //
  HDFArray<int> data1;
  data1.Create(hdfFile.rootGroup, "data1");
  int value = 3;
  data1.Write(&value, 1);

  intAtom.Create(data1.dataset, "intatom");
  intAtom.Write(300);
  string tempval("hello world");  
  //  strAtom.Create(data1.dataset, "stratomname", tempval);
  strAtom.Create(data1.dataset, "stratomname");
  strAtom.Write(tempval);
  string tempval2("bonjour monde");
  strAtom.Write(tempval2);
  string tempval3("hola mundo");
  strAtom.Write(tempval3);
  
  HDFAtom<vector<string> > strVectAtom;
  vector<string> stringVectValues;
  stringVectValues.push_back(tempval);
  stringVectValues.push_back(tempval2);
  stringVectValues.push_back(tempval3);
  strVectAtom.Create(data1.dataset, "strvectatom", stringVectValues);

  HDFAtom<vector<int>* > intVectAtom;
  vector<int> values;
  values.push_back(100);
  values.push_back(200);
  values.push_back(300);
  values.push_back(400);
  intVectAtom.Create(data1.dataset, "intvectatom", values);

	return 0;
}
