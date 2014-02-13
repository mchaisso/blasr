#ifndef DATA_HDF_CMP_WRITER_H_
#define DATA_HDF_CMP_WRITER_H_

#include "FASTASequence.h"
#include "NucConversion.h"
#include <iostream>
#include <iomanip>
#include "qvs/QualityValue.h"
#include "datastructures/matrix/Matrix.h"
#include "cmpseq/CompressedSequence.h"
using namespace std;

class HDFCmpWriter : public HDFCmpData {
  HDFFile outFile;
  string  outFileName;
  static const int NCols=22;
  Group  rootGroup;

  const char  *colNameIds[NCols] = {"00", "01", "02", "03", "04", "05", "06", "07", "08", "09",
                                    "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                    "20", "21"};
  vector<HDFAtom<string> > colNameAtoms;
  HDF2DArray<int> alignmentIndexArray;
  vector<HDFCmpRefAlignmentGroup*> refAlignGroups;
  static const int bufferSize = 256;
 public:
  
  void Create(string &outFileNameP) {
    outFileName = outFileNameP;
    outFile.Create(outFileName);
    outFile.OpenGroup("/", rootGroup);
  }

  void AddCommandLine(string commandLineStr) {
    commandLine.Create(rootGroup, "CommandLine", commandLineStr);
  }

  void CreateAlignmentIndexArray() {
    alignmentIndexArray.Initialize(&rootGroup, "AlignmentIndex", NCols, NCols*10);
  }

  void AddColumnNames(CmpFile &cmpFile) {
    int colIndex;
    colNameAtoms.resize(NCols);
    for (colIndex = 0; colIndex < NCols; colIndex++ ) {
      colNameStrm.str("");
      colNameStrm << "ColName" << colNameIds[colIndex];
      string colNameString = colNameStrm.str();
      colNameAtoms[colIndex].Create(alignmentIndexArray.dataset, colNameString);
      colNameAtoms[colIndex].Write(cmpFile.colNames[colIndex]);
    }
  }

  void WriteTables(CmpFile &cmpFile) {
    movieNameIdArray.Initialize(rootGroup, "MovieId");
    movieNameArray.Initialize(&rootGroup, "MovieName");
    movieNameIdArray.Write(&cmpFile.movieNameTable.movieIds[0], 
                           cmpFile.movieNameTable.movieIds.size());
    
    int movieNameIndex;
    for (movieNameIndex = 0 ; movieNameIndex < cmpFile.movieNameTable.movieIds.size(); movieNameIndex++) {
      movieNameArray.Write(cmpFile.movieNameTable.movieNames[i].c_str(), 1);
    }

    readGroupPathIdArray.Initialize(&rootGroup, "ReadGroupPath");
    readGroupPathIdArray.Write(&cmpFile.readGroupTable.readGroupIds[0],
                               cmpFile.readGroupTable.readGroupIds.size());

    readGroupPathArray.Initialize(&rootGroup, "ReadGroupPath");
    int readGroupIndex;
    for (readGroupIndex = 0; readGroupIndex < cmpFile.readGroupTable.readGroupIds.size(); readGroupIndex++) {
      readGroupPathArray.Write(cmpFile.readGroupTable.readGroupNames[readGroupIndex].c_str(), 1);
    }
    
    refSeqNameIdArray.Initialize(&rootGroup, "RefSeqID");
    refSeqNameIdArray.Write(&cmpFile.refSeqNameTable.refSeqNameIds[0]
                            cmpFile.refSeqNameTable.refSeqNameIds.size());
    
    refSeqNameArray.Initialize(&rootGroup, "RefSeqName");
    int refSeqNameIndex;
    for (refSeqNameIndex = 0; refSeqNameIndex < cmpFile.refSeqNameTable.refSeqNames.size(); refSeqNameIndex++) {
      refSeqNameArray.Write(cmpFile.refSeqNameTable.refSeqNames[refSeqNameIndex].c_str(), 1);
    }
  }

  void WriteAlignmentIndex(CmpFile &cmpFile, HDFPlsFile pulseFile) {
    /*
     * Write out all alignment indices to the file.
     */
    int i;
    vector<int> alignmentIndex;
    for (i = 0; i < cmpFile.alignments.size(); i++ ){
      alignmentIndexArray.WriteRow(cmpFile.alignments.GetAlignmentIndex(), cmpFile.alignments[i].GetAlignmentIndexSize());
    }
  }

  void WriteCmpFile(CmpFile &cmpFile) {
    int i;

  }
  
};
  /*
    def getPulseDataFieldsForHoleNumber(self, fields, hn, rpo=True):
        rtn = {}
        channel = self._getPulseFieldForHoleNumber("Channel", hn, rpo)
        starts = 1.0*self._getPulseFieldForHoleNumber("StartFrame", hn, rpo) #need to cast the unsigned integer type to singed type for IPD calculation
        widths = 1.0*self._getPulseFieldForHoleNumber("WidthInFrames", hn, rpo) #need to cast the unsigned integer type to singed type for IPD calculation
        for field in fields:
            if field in [ "QualityValue", "DeletionQV" ]:
                rtn[field] = self._getBaseCallAnnotationFromHoleNumber(field, hn)[1]
            elif field not in COMPUTED_FIELDS:
                rtn[field] = self._getPulseFieldForHoleNumber(field, hn, rpo)
            else:
                if field == 'StartTime':
                    rtn[field] = starts / self._metaData.FrameRate

                elif field == 'PulseWidth':
                    rtn[field] = widths / self._metaData.FrameRate

                elif field == 'pkmid':
                    midSignal = self._getPulseFieldForHoleNumber("MidSignal", hn, rpo)
                    rtn[field] = self._createArrayByChannel(midSignal, channel)

                elif field == 'pkmax':
                    maxSignal = self._getPulseFieldForHoleNumber("MaxSignal", hn, rpo)
                    rtn[field] = self._createArrayByChannel(maxSignal, channel)

                elif field == 'IPD':
                    rtnArray = numpy.zeros( len(starts), dtype=numpy.float )
                    rtnArray[1:] = (starts[1:] - starts[:-1] - widths[:-1]) / self._metaData.FrameRate
                    rtn[field] = self._flatten(rtnArray)
                    
                elif field == 'Light':
                    meanSignal = self._getPulseFieldForHoleNumber("MeanSignal", hn, rpo)
                    pkmean = self._createArrayByChannel(meanSignal, channel)
                    rtn[field] = self._flatten( pkmean * widths ) 
        
        return rtn
  */
  
};


#endif
