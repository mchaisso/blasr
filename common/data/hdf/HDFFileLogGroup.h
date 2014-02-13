#ifndef HDF_FILE_LOG_GROUP_H_
#define HDF_FILE_LOG_GROUP_H_
#include "HDFGroup.h"
#include "HDFArray.h"

class HDFFileLogGroup {
 public:
  HDFGroup group;
  HDFArray<string> commandLineArray;
  HDFArray<string> versionArray;
  HDFArray<string> timestampArray;
  HDFArray<unsigned int> idArray;
  HDFArray<string> logArray;
  HDFArray<string> programArray;

  int Initialize(HDFGroup &parentGroup) {
    if (group.Initialize(parentGroup.group, "FileLog") == 0) { return 0; }
    commandLineArray.Initialize(group, "CommandLine");
    versionArray.Initialize(group, "Version");
    timestampArray.Initialize(group, "Timestamp");
    idArray.Initialize(group, "ID");
    logArray.Initialize(group, "Log");
    programArray.Initialize(group, "Program");
  }

  
  void AddEntry(string command, string log, string program, string timestamp, string version) {
    commandLineArray.Write(&command, 1);
    versionArray.Write(&version, 1);
    timestampArray.Write(&timestamp, 1);
    programArray.Write(&program, 1);
    logArray.Write(&log, 1);
    
    unsigned int id = idArray.size();
    id = id + 1;
    idArray.Write(&id, 1);
    
  }

  bool Create(HDFGroup &parent) {
    parent.AddGroup("FileLog");
    if (group.Initialize(parent.group, "FileLog") == 0) { return 0; }
    commandLineArray.Create(group, "CommandLine");
    versionArray.Create(group, "Version");
    timestampArray.Create(group, "Timestamp");
    programArray.Create(group, "Program");
    logArray.Create(group, "Log");
    idArray.Create(group, "ID");
  }    

};

#endif
