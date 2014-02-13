#ifndef DATASTRUCTURES_READS_READ_TYPE_H_
#define DATASTRUCTURES_READS_READ_TYPE_H_

#include <string>
using namespace std;
class ReadType {
 public:
  typedef enum E_ReadTypeEnum {NoReadType, Standard, CCS, RCCS} ReadTypeEnum;
};

ReadType::ReadTypeEnum ParseReadType(string &readTypeString) {
  if (readTypeString == "Standard") {
    return ReadType::Standard;
  }
  else if (readTypeString == "CCS") {
    return ReadType::CCS;
  }
  else if (readTypeString == "RCCS") {
    return ReadType::RCCS;
  }
  else {
      return ReadType::NoReadType;
  }
};

#endif
