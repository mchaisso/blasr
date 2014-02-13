#include "../common/data/hdf/HDFArray.h"
#include "../common/data/hdf/HDFFile.h"
#include "../common/CommandLineParser.h"


int main(int argc, char* argv[]) {

  CommandLineParser clp;
  int pos = 0;
  int value = 10;
  bool openTrunc = false;
  clp.RegisterFlagOption("trunc", &openTrunc, "Open file in truncated mode.", false);
  clp.RegisterIntOption("value", &value, "Store this value.", CommandLineParser::Integer, false);
  clp.RegisterIntOption("pos", &pos, "Store at this position", CommandLineParser::PositiveInteger, false);
  clp.ParseCommandLine(argc, argv);

	HDFFile hdfFile;
  string name("test_write_int_array.hdf");
  if (openTrunc) {
    hdfFile.Open(name, H5F_ACC_TRUNC);
  }
  else {
    hdfFile.Open(name, H5F_ACC_RDWR);
  }

  HDFArray<int> intArray;
  HDFArray<string> strArray;
  intArray.Initialize(hdfFile.rootGroup, "IntData");
  intArray.WriteToPos(&value, 1, pos);
  hdfFile.Close();
  
  strArray.Initialize(hdfFile.rootGroup, "StringData");
  string hello("hello"), joe("joe");
  strArray.Write(&hello, 1);
  strArray.Write(&joe, 1);
  /*
  //
  // Try this on an already created file.
  hdfFile.Create(name);
  intArray;
  intArray.Create(hdfFile.rootGroup, "IntData");
  value=100;
  intArray.WriteToPos(&value, 1, 10);
  hdfFile.Close();
  */
};
