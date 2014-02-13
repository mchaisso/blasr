#define __FAST_MATH__

#include "data/hdf/HDFCmpFile.h"
#include "data/hdf/HDFBasReader.h"
#include "data/hdf/HDFPlsReader.h"
#include "data/hdf/HDFCCSReader.h"
#include "data/hdf/PlatformId.h"
#include "datastructures/alignment/CmpFile.h"
#include "datastructures/alignment/CmpAlignment.h"
#include "datastructures/alignment/ByteAlignment.h"
#include "datastructures/alignment/AlignmentMap.h"
#include "datastructures/reads/BaseFile.h"
#include "datastructures/reads/PulseFile.h"
#include "datastructures/reads/ReadType.h"
#include "datastructures/loadpulses/MetricField.h"
#include "datastructures/loadpulses/MovieAlnIndexLookupTable.h"
#include "utils/FileOfFileNames.h"
#include "utils/TimeUtils.h"
#include "files/BaseSequenceIO.h"
#include "CommandLineParser.h"
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <assert.h>
#include <numeric>
#include <stdio.h>

using namespace std;

typedef map<string, int> MovieNameToArrayIndex;
typedef map<string, bool> MetricOptionsMap;
typedef map<string, vector<string> > RequirementMap;


char VERSION[] = "v1.1.0";
char PERFORCE_VERSION_STRING[] = "$Change: 126407 $";

// define default values for metrics
const float NaN = 0.0/0.0;
const UChar missingQualityValue      = 255;
const unsigned char maxQualityValue  = 100;
const HalfWord missingFrameRateValue = USHRT_MAX;
const unsigned int missingPulseIndex = UINT_MAX;


void CapQualityValue(QualityValueVector<QualityValue> &vect, DNALength length, unsigned char maxQualityValue=100) {
    unsigned int i;
    if (vect.data == NULL) {
        return;
    }
    for (i = 0; i < length; i++) {
        vect.data[i] = min(vect.data[i], maxQualityValue);
    }
}


void CapQualityValues(SMRTSequence &seq, unsigned char maxQualityValue = 100) {
    CapQualityValue(seq.qual, seq.length, maxQualityValue);
    CapQualityValue(seq.deletionQV, seq.length, maxQualityValue);
    CapQualityValue(seq.preBaseDeletionQV, seq.length, maxQualityValue);
    CapQualityValue(seq.insertionQV, seq.length, maxQualityValue);
    CapQualityValue(seq.substitutionQV, seq.length, maxQualityValue);
    CapQualityValue(seq.mergeQV, seq.length, maxQualityValue);
}


int CheckCmpFileFormat(CmpFile &cmpFile) {
    if (cmpFile.readType != ReadType::Standard) {
        cout << "ERROR! Reading pulse information into a cmp.h5 file generated from circular " << endl
            << "consensus called sequences is not supported." << endl;
        exit(1);
    }
    return 1;
}

void BuildRequirementMap(RequirementMap &fieldRequirements) {
    fieldRequirements["StartTimeOffset"].push_back("StartFrame");
    fieldRequirements["StartTimeOffset"].push_back("NumEvent");
    fieldRequirements["StartFrame"].push_back("PreBaseFrames");
    fieldRequirements["StartFrame"].push_back("WidthInFrames");
    fieldRequirements["PulseWidth"].push_back("WidthInFrames");
    fieldRequirements["pkmid"].push_back("MidSignal");
    fieldRequirements["pkmid"].push_back("NumEvent");
    fieldRequirements["IPD"].push_back("StartFrame");
    fieldRequirements["IPD"].push_back("NumEvent");
    fieldRequirements["IPD"].push_back("PreBaseFrames");
    fieldRequirements["IPD"].push_back("WidthInFrames");
    fieldRequirements["Light"].push_back("MeanSignal");
    fieldRequirements["Light"].push_back("NumEvent");
    fieldRequirements["Light"].push_back("WidthInFrames");

    // Build requirementMap for sneaky metrics
    fieldRequirements["StartFrameBase"].push_back("PreBaseFrames");
    fieldRequirements["StartFrameBase"].push_back("WidthInFrames");
    fieldRequirements["StartFramePulse"].push_back("PreBaseFrames");
    fieldRequirements["StartFramePulse"].push_back("WidthInFrames");
}

void ExclusivelyAdd(const char *value, vector<string> &vect) {
    if (find(vect.begin(), vect.end(), value) == vect.end()) {
        vect.push_back(value);
    }
}

bool AnyFieldRequiresFrameRate(vector<string> &fields) {
    int i;
    for (i = 0; i < fields.size(); i++ ) {
        if (fields[i] == "PulseWidth" or
                fields[i] == "IPD" or
                fields[i] == "Light" or
                fields[i] == "StartTimeOffset" or
                fields[i] == "StartFrame" or
                fields[i] == "PulseWidth" or 
                fields[i] == "PreBaseFrames" or 
                fields[i] == "WidthInFrames") {
            return true;
        }
    }
    return false;
}

template<typename T>
void Free(T* &buf) {
    if (buf != NULL){ 
        delete[] buf;
    }
    buf = NULL;
}

// Return all eighteen metrics that can be loaded.
// StartTimeOffset  QualityValue    InsertionQV   MergeQV    
// DeletionQV       DeletionTag     PulseIndex    SubstitutionTag
// SubstitutionQV   ClassifierQV    StartFrame    PulseWidth
// PreBaseFrames    WidthInFrames   pkmid         IPD
// Light            WhenStarted
vector<string> GetAllSupportedMetrics(bool isSneakyMetricsIncluded = true) {
    // The order of metrics matters. With -bymetric option, all fields 
    // which are required for computing a metric are cached before WriteMetric()
    // and cleared afterwards. If two neighboring metrics share a subset of 
    // required fields, then the cached fields can be re-used. Arrange metrics
    // in an order that maximizes reuse of cached fields.
    vector<string> supportedMetrics;
    supportedMetrics.push_back("WhenStarted");

    supportedMetrics.push_back("QualityValue");
    supportedMetrics.push_back("InsertionQV");
    supportedMetrics.push_back("MergeQV");
    supportedMetrics.push_back("DeletionQV");
    supportedMetrics.push_back("DeletionTag");
    supportedMetrics.push_back("SubstitutionTag");
    supportedMetrics.push_back("SubstitutionQV");
    supportedMetrics.push_back("PreBaseFrames");
    // Sneaky metrics for internal use Only
    if (isSneakyMetricsIncluded) {
        supportedMetrics.push_back("StartFrameBase");
    }
    supportedMetrics.push_back("IPD");
    supportedMetrics.push_back("StartFrame");
    if (isSneakyMetricsIncluded) {
        supportedMetrics.push_back("StartFramePulse");
    }
    // Disable metric StartTimeOffset for now.
    // StartTimeOffset is placed at the same level as AlnArray, However, the 
    // size of StartTimeOffset is far less than AlnArray, while cmp.h5 spec
    // requires all datasets at that level to have the same size.  
    
    // supportedMetrics.push_back("StartTimeOffset");

    supportedMetrics.push_back("PulseWidth");
    supportedMetrics.push_back("WidthInFrames");
    supportedMetrics.push_back("Light");
    supportedMetrics.push_back("pkmid");
    supportedMetrics.push_back("ClassifierQV");
    supportedMetrics.push_back("PulseIndex");

    return supportedMetrics;
}

// Return metrics to load by default.
vector<string> GetDefaultMetrics() {
    vector<string> defaultMetrics;
    defaultMetrics.push_back("QualityValue");
    defaultMetrics.push_back("ClassifierQV");
    defaultMetrics.push_back("StartFrame");
    defaultMetrics.push_back("PulseWidth");
    defaultMetrics.push_back("WidthInFrames");
    defaultMetrics.push_back("pkmid");
    defaultMetrics.push_back("IPD");
    return defaultMetrics;
}

// Return metrics that can be computed from PulseCalls.
vector<string> GetPulseMetrics() {
    vector<string> pulseMetrics;
    pulseMetrics.push_back("StartFrame");
    pulseMetrics.push_back("StartTimeOffset");
    pulseMetrics.push_back("ClassifierQV");
    pulseMetrics.push_back("PulseWidth");
    pulseMetrics.push_back("WidthInFrames");
    pulseMetrics.push_back("IPD");
    pulseMetrics.push_back("pkmid");
    pulseMetrics.push_back("Light");
    pulseMetrics.push_back("StartFramePulse");
    return pulseMetrics;
}

// Return true if this metric can be computed from PulseCalls.
bool IsPulseMetric(const string & metric) {
    vector<string> pulseMetrics = GetPulseMetrics();
    for (int i = 0; i < pulseMetrics.size(); i++) {
        if (pulseMetrics[i] == metric) 
            return true;
    }
    return false;
}

// Return all metrics that are 
// (1) supported, 
// (2) requested to load, and 
// (3) computable with all required fields available 
//     in either bas.h5 or pls.h5.
vector<string> GetMetricsToLoad(map<string, bool> & metricOptions) {
    vector<string> metricsToLoad; 
    // Get all supported metrics. 
    vector<string> supportedMetrics = GetAllSupportedMetrics();
    map<string, bool>::iterator metricIt;
    for (int i = 0; i < supportedMetrics.size(); i++) {
        string metric = supportedMetrics[i];
        metricIt = metricOptions.find(metric);
        if (metricIt!=metricOptions.end() and metricIt->second) {
            // Get metrics that are required and computable
            metricsToLoad.push_back(metricIt->first);
        }
    }
    return metricsToLoad; 
}

void StoreDatasetFieldsFromPulseFields(MetricOptionsMap &fieldSet,
        RequirementMap &fieldRequirements, 
        vector<string> &datasetFields) {
    int f;
    int d;
    MetricOptionsMap::iterator optionsIt;
    for (optionsIt = fieldSet.begin(); optionsIt != fieldSet.end(); ++optionsIt) {
        if (optionsIt->second == true) {
            if (fieldRequirements.find(optionsIt->first) == fieldRequirements.end()) {
                ExclusivelyAdd(optionsIt->first.c_str(), datasetFields);
            }
            else {
                for (d = 0; d < fieldRequirements[optionsIt->first].size(); d++) {
                    ExclusivelyAdd(fieldRequirements[optionsIt->first][d].c_str(), datasetFields );
                }
            }
        }
    }
}

void ParseMetricsList(string metricListString, MetricOptionsMap &metricOptions) {
    vector<string> metrics;
    Tokenize(metricListString, ",", metrics);
    int m;
    for  (m = 0; m < metrics.size(); m++) {
        if (metricOptions.find(metrics[m]) != metricOptions.end()) {
            metricOptions[metrics[m]] = true;
        }
        else {
            cout << "ERROR! Metric " << metrics[m] << " is not supported." << endl;
            exit(1);
        }
    }
}

// Set default metric options to true
void SetDefaultMetricOptions(map<string, bool> & metricOptions) {
    vector<string> defaultMetrics = GetDefaultMetrics(); 
    for (int i = 0; i < defaultMetrics.size(); i++) {
        metricOptions[defaultMetrics[i]] = true;
    }
}


// Initialize all supported metric options and set all to false
void CreateMetricOptions(map<string, bool> &metricOptions) {
    vector<string> supportedMetrics = GetAllSupportedMetrics();
    for (int i = 0; i < supportedMetrics.size(); i++) {
        metricOptions[supportedMetrics[i]] = false;
    }
}

// Check whether all fields are available or not.
bool AreAllFieldsAvailable(
        vector <Field>   & requiredFields,
        HDFBasReader     & hdfBasReader,
        HDFPlsReader     & hdfPlsReader, 
        const bool       & useBaseFile,
        const bool       & usePulseFile) {
    bool allAvailable = true;

    for (int i = 0; i < requiredFields.size(); i++) {
        Field field = requiredFields[i];
        if (field.type == BasField) {
            if (!useBaseFile or !hdfBasReader.FieldIsIncluded(field.name)
                or !hdfBasReader.includedFields[field.name]) {
                allAvailable = false;
                break;
            }
        } else if (field.type == PlsField) {
            if (!usePulseFile or !hdfPlsReader.FieldIsIncluded(field.name)
                or !hdfPlsReader.includedFields[field.name]) {
                allAvailable = false;
                break;
            }
        }
    }
    return allAvailable;
}

//
// Check whether a metric is computable or not.
// fieldsToBeUsed = all fields that will be used for computing a metric. 
// If a metric can be computed from both bas and pls files (e.g. 
// StartFrame, IPD, PulseWidth, WidthInFrame), only compute it from pls.
//
bool CanThisMetricBeComputed ( 
        const string     & metricName,
        HDFBasReader     & hdfBasReader,
        HDFPlsReader     & hdfPlsReader, 
        const bool       & useBaseFile,
        const bool       & usePulseFile,
        vector<Field>    & fieldsToBeUsed) {
    fieldsToBeUsed.clear();

    FieldsRequirement fieldsRequirement = FieldsRequirement(metricName);

    bool metricMayBeComputedFromPls = true;
    if (fieldsRequirement.fieldsUsePlsFile.size() != 0 && usePulseFile) {
        metricMayBeComputedFromPls = AreAllFieldsAvailable(
                fieldsRequirement.fieldsUsePlsFile, 
                hdfBasReader, hdfPlsReader,
                useBaseFile, usePulseFile);
    } else {
        metricMayBeComputedFromPls = false;
    }

    bool metricMayBeComputedFromBas = true;
    if (fieldsRequirement.fieldsUseBasFile.size() != 0 && useBaseFile) {
        metricMayBeComputedFromBas = AreAllFieldsAvailable(
                fieldsRequirement.fieldsUseBasFile, 
                hdfBasReader, hdfPlsReader,
                useBaseFile, usePulseFile);
    } else {
        metricMayBeComputedFromBas = false;
    }

    bool metricMayBeComputed = true;
    if (!metricMayBeComputedFromBas and !metricMayBeComputedFromPls) {
        metricMayBeComputed = false;
    }
    
    // Compute from pls if possible
    if (metricMayBeComputedFromPls) {
        fieldsToBeUsed = fieldsRequirement.fieldsUsePlsFile;
    } else if (metricMayBeComputedFromBas) {
        fieldsToBeUsed = fieldsRequirement.fieldsUseBasFile;
    }

    if (metricName == "StartTimeOffset") {
        metricMayBeComputed = false;
        // Disable StartTimeOffset for now.
    }
    if (metricName == "WhenStarted") {
        // WhenStarted requires no fields from neither bas nor pls.
        metricMayBeComputed = true;
    }

    return metricMayBeComputed; 
}

//
// Check whether metrics are computable or not. If a metric is not
// computable, disable it with a warning or exit with an error.
//
void CanMetricsBeComputed(
        MetricOptionsMap & metricOptions,
        HDFBasReader     & hdfBasReader,
        HDFPlsReader     & hdfPlsReader, 
        const bool       & useBaseFile,
        const bool       & usePulseFile,
        const bool       & failOnMissingData, 
        const string     & movieName) {

    map<string,bool>::iterator metricIt;
    for (metricIt = metricOptions.begin(); metricIt != metricOptions.end(); ++metricIt) {
        string metricName = metricIt->first;
        if (metricName == "") {
            metricIt->second == false;
        }

        if (metricIt->second == false) {
            continue;
        }
        vector<Field> fieldsToBeUsed;
        bool metricMayBeComputed = CanThisMetricBeComputed(metricName, 
                hdfBasReader, hdfPlsReader, useBaseFile, usePulseFile, 
                fieldsToBeUsed);

        if (metricMayBeComputed == false) {
            if (failOnMissingData) {
                cout << "ERROR";
            }
            else {
                cout << "WARNING";
            }
            cout << ": There is insufficient data to compute metric: " 
                 << metricName << " in the file " << movieName << " ";
            cout << " It will be ignored." << endl;
            if (failOnMissingData) {
                exit(1);
            }
            metricOptions[metricName] = false;
        }
    }
}

// Return size of a single field in KB.
UInt ComputeRequiredMemoryForThisField(
        Field          & thisField,
        HDFBasReader   & hdfBasReader,
        HDFPlsReader   & hdfPlsReader,
        const bool     & useBaseFile,
        const bool     & usePulseFile) {
    UInt memory = 0;
    if (thisField.type == BasField) {
        assert(useBaseFile);
        return hdfBasReader.GetFieldSize(thisField.name);
    }
    if (thisField.type == PlsField) {
        assert(usePulseFile);
        return hdfPlsReader.GetFieldSize(thisField.name);
    }
    assert(false);
}

//
// Return estimated memory peak (in KB) for buffering all data using -bymetric.
//
UInt ComputeRequiredMemory(
        vector<string> & metricsToLoad, 
        HDFBasReader   & hdfBasReader,
        HDFPlsReader   & hdfPlsReader,
        const bool     & useBaseFile,
        const bool     & usePulseFile,
        HDFCmpFile<CmpAlignment> & cmpReader,
        UInt           & totalAlnLength) {
    UInt maxMemory = 0;
    for (int i = 0; i < metricsToLoad.size(); i++) {
        UInt memoryForThisMetric = 0;
        vector<Field> fieldsToBeUsed;
        bool canBeComputed = CanThisMetricBeComputed(
                metricsToLoad[i], hdfBasReader, hdfPlsReader,
                useBaseFile, usePulseFile, fieldsToBeUsed);

        for (int j = 0; j < fieldsToBeUsed.size(); j++) {
            UInt memoryForThisField = ComputeRequiredMemoryForThisField(
                    fieldsToBeUsed[j], hdfBasReader, hdfPlsReader,
                    useBaseFile, usePulseFile);
            memoryForThisMetric += memoryForThisField;
        }
        maxMemory = max(maxMemory, memoryForThisMetric);
    }
    //
    // AlnIndex will be buffered. Some other datastructures also need  
    // to be buffered for quick look up. Approximately double the size.
    //
    UInt totalAlnIndexMem = 2 * cmpReader.alnInfoGroup.GetAlnIndexSize();

    //
    // AlnArray and metrics to load needs to be buffered in KB.  
    //
    UInt totalAlnArrayMem = totalAlnLength / 1024 *
                            (sizeof(unsigned int) + sizeof(unsigned char));

    //
    // It's diffcult to estimate how much memory will be used by hdf5.
    // Assume memory consumed by hdf5 scales with AlnIndex and AlnArray datasets. 
    //
    UInt hdf5Mem = totalAlnIndexMem / 2 + totalAlnLength / 1024 * sizeof(unsigned int); 

    maxMemory += totalAlnIndexMem + totalAlnArrayMem + hdf5Mem;

    //cout << "The estimated peak memory for buffering fields is " 
    //     << maxMemory << " KB." << endl;
    //cout << "The estimated memory for buffering AlnIndex related data is " 
    //     << totalAlnIndexMem << " KB."<< endl;
    //cout << "The estimated memory for buffering AlnArray related data is " 
    //     << totalAlnArrayMem << " KB." << endl;
    //cout << "The estimated memory for hdf5 is "
    //     << hdf5Mem << " KB." << endl; 
    //cout << "The estimated total memory is " 
    //     << maxMemory << " KB." << endl;
    return maxMemory;
}

//
// Get aligned sequence for this alignment from cmpFile
//
string GetAlignedSequenceFromCmpFile(
        const HDFCmpFile<CmpAlignment> & cmpReader,
        MovieAlnIndexLookupTable       & lookupTable) {
    string alignedSequence;
    vector <unsigned char> byteAlignment;
    int alignedSequenceLength = lookupTable.offsetEnd - lookupTable.offsetBegin;
    if (alignedSequenceLength >= 0 ) {
        alignedSequence.resize(alignedSequenceLength);
        byteAlignment.resize(alignedSequenceLength);
    }
    //
    // Read the alignment string.  All alignments 
    //
    cmpReader.refAlignGroups[lookupTable.refGroupIndex]->readGroups[lookupTable.readGroupIndex]->alignmentArray.Read(
        lookupTable.offsetBegin, 
        lookupTable.offsetEnd, 
        &byteAlignment[0]);
    
    //
    // Convert to something we can compare easily.
    //
    ByteAlignmentToQueryString(&byteAlignment[0], byteAlignment.size(), &alignedSequence[0]);
    return alignedSequence;
}

//
// Store info necessary for loading pulses to lookupTable.
//
void BuildLookupTable(
    const int                        & movieAlignmentIndex,
    CmpFile                          & cmpFile,
    BaseFile                         & baseFile,
    const bool                       & usePulseFile,
    PulseFile                        & pulseFile,
    HDFCmpFile<CmpAlignment>         & cmpReader,
    const vector<int>                & movieAlnIndex, 
    const vector< pair<int,int> >    & toFrom,
    const set<uint32_t>              & moviePartHoleNumbers,
    MovieAlnIndexLookupTable         & lookupTable) {
    //
    // Query the cmp file for a way to look up a read based on
    // coordinate information.  For Astro reads, the coords are
    // based on x and y.  For Springfield, it is read index.  The
    // base files should be able to look up reads by x,y or by
    // index. 
    //
    if (cmpFile.platformId == Astro) {
        cout << "ASTRO pulse loading is deprecated." << endl;
        exit(1);
    }

    int alignmentIndex = movieAlnIndex[toFrom[movieAlignmentIndex].second];

    //
    // Alignments are grouped by ref group id then movie id.
    //
    int refGroupId  = cmpFile.alnInfo.alignments[alignmentIndex].GetRefGroupId();
    int movieId     = cmpFile.alnInfo.alignments[alignmentIndex].GetMovieId();
    UInt holeNumber = cmpFile.alnInfo.alignments[alignmentIndex].GetHoleNumber();
    int alnGroupId  = cmpFile.alnInfo.alignments[alignmentIndex].GetAlnGroupId();

    if (cmpReader.refGroupIdToArrayIndex.find(refGroupId) == 
        cmpReader.refGroupIdToArrayIndex.end()) {
        cout << "ERROR! An alignment " << alignmentIndex 
             << " is specified with reference group " << endl
             << refGroupId << " that is not found as an alignment group." << endl;
        exit(1);
    }
    int refGroupIndex = cmpReader.refGroupIdToArrayIndex[refGroupId];

    //
    // Now find the group containing the alignment.
    //
    if (cmpReader.alnGroupIdToReadGroupName.find(alnGroupId) == 
        cmpReader.alnGroupIdToReadGroupName.end()) {
        cout << "ERROR! An alignment " << alignmentIndex 
             << " is specified with alignment group " << endl
             << alnGroupId << " that is not found." << endl;
        exit(1);
    }

    string readGroupName = cmpReader.alnGroupIdToReadGroupName[alnGroupId];
    if (cmpReader.refAlignGroups[refGroupIndex]->experimentNameToIndex.find(readGroupName) == 
        cmpReader.refAlignGroups[refGroupIndex]->experimentNameToIndex.end()) {
        cout << "ERROR! An alignment " << alignmentIndex 
             << " is specified with read group name " << endl
             << readGroupName << " that is not found." << endl;
        exit(1);
    }

    int readGroupIndex = cmpReader.refAlignGroups[refGroupIndex]->experimentNameToIndex[readGroupName];
   
    UInt offsetBegin = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetBegin();
    UInt offsetEnd   = cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetEnd();

    //
    // First pull out the bases corresponding to this read.
    //
    int queryStart = cmpFile.alnInfo.alignments[alignmentIndex].GetQueryStart();
    int queryEnd   = cmpFile.alnInfo.alignments[alignmentIndex].GetQueryEnd();

    bool skip = false;
    int readIndex, readStart, readLength, plsReadIndex; 
    readIndex = readStart = readLength = plsReadIndex = -1;
    //
    // Since the movie may be split into multiple parts, look to see
    // if this hole number is one of the ones covered by this
    // set. If it is not, just continue. It will be loaded on
    // another pass through a different movie part.
    //
    if (moviePartHoleNumbers.find(holeNumber) == moviePartHoleNumbers.end()) {
        skip = true;
    } else {
        if (!baseFile.LookupReadIndexByHoleNumber(holeNumber, readIndex)) {
            cout << "ERROR! Alignment has hole number " << holeNumber 
                 << " that is not in the movie. " << endl;
            exit(1); 
        }
        readStart  = baseFile.readStartPositions[readIndex];
        readLength = baseFile.readStartPositions[readIndex+1] - 
                     baseFile.readStartPositions[readIndex];
        if (usePulseFile) {
            if (!pulseFile.LookupReadIndexByHoleNumber(holeNumber, plsReadIndex)) {
                cout << "ERROR! Alignment has  hole number " << holeNumber 
                     << " that is not in the movie. " << endl;
                exit(1);
            }
            assert(pulseFile.holeNumbers[plsReadIndex] ==
                   baseFile.holeNumbers[readIndex]);
        }
    }
    // Save info to lookupTable
    lookupTable.SetValue(skip, // Skip processing this or not
        movieAlignmentIndex, 
        alignmentIndex,  
        refGroupIndex,       
        readGroupIndex,
        holeNumber,    // cmp.h5 /AlnInfo/AlnIndex column 7 
        offsetBegin,   // cmp.h5 /AlnInfo/AlnIndex column 18
        offsetEnd,     // cmp.h5 /AlnInfo/AlnIndex column 19
        queryStart,    // cmp.h5 /AlnInfo/AlnIndex column 11 
        queryEnd,      // cmp.h5 /AlnInfo/AlnIndex column 12 
        readIndex,     // hole Index in BaseCalls/ZMW/HoleNumber
        readStart,     // readStart in BaseCalls/* (e.g. *=Basecall)
        readLength,    // readLength in BaseCalls/*
        plsReadIndex); // readIndex in PulseCalls/ZMW/HoleNumber
}

//
// Map bases of a read to pulse indices.
//
void MapBaseToPulseIndex(
        BaseFile                 & baseFile, 
        PulseFile                & pulseFile,
        MovieAlnIndexLookupTable & table,
        vector<int>              & baseToPulseIndexMap) {
    baseToPulseIndexMap.resize(table.readLength);

    int pulseStart = pulseFile.pulseStartPositions[table.plsReadIndex];
    //
    // Copy the subset of pulses that correspond to the ones called as bases.
    //
    int i;
    for (i = 0; i < table.readLength; i++) {
        baseToPulseIndexMap[i] = pulseStart + baseFile.pulseIndex[table.readStart + i];
    }
}

//
// Get source read from the bas/pls file.
//
void GetSourceRead(CmpFile      & cmpFile,      
                   BaseFile     & baseFile, 
                   PulseFile    & pulseFile,  
                   HDFBasReader & hdfBasReader, 
                   HDFPlsReader & hdfPlsReader,
                   HDFCCSReader<SMRTSequence> & hdfCcsReader,
                   const bool   & useBaseFile,  
                   const bool   & usePulseFile,
                   const bool   & useCcsOnly,
                   //const bool   & byRead, 
                   MovieAlnIndexLookupTable & table, 
                   const string & alignedSequence,
                   SMRTSequence & sourceRead,   
                   unsigned int & numPasses) {

    assert(!table.skip);
    //
    // These are not allocated in the regular allocate function
    // since they are only used in loadPulses. (maybe I should
    // subclass SMRTSequence here).
    //

    // Read in the data from the bas file if it exsts.
    if (useBaseFile) {
        hdfBasReader.GetReadAt(table.holeNumber, sourceRead);
        if (cmpFile.readType == ReadType::CCS or useCcsOnly) {
            numPasses = hdfCcsReader.GetNumPasses(table.holeNumber);
        }
    }
    // Read in the data from the pls file if it exists.
    if (usePulseFile) {
        hdfPlsReader.GetReadAt(table.plsReadIndex, sourceRead.pulseIndex, sourceRead);
    }
    CapQualityValues(sourceRead);
}

//
// Build lookup tables for all alignments whose indices in 
// AlnArray are saved in movieAlnIndex.
// Also check whether the bas file and the cmp file match.
//
void BuildLookupTablesAndMakeSane(
        CmpFile                          & cmpFile,
        BaseFile                         & baseFile,
        PulseFile                        & pulseFile,
        HDFCmpFile<CmpAlignment>         & cmpReader,
        HDFBasReader                     & hdfBasReader,
        HDFPlsReader                     & hdfPlsReader,
        HDFCCSReader<SMRTSequence>       & hdfCcsReader,
        const bool                       & useBaseFile,
        const bool                       & usePulseFile,
        const bool                       & useCcsOnly,
        const vector<int>                & movieAlnIndex, 
        const vector< pair<int,int> >    & toFrom,
        const set<uint32_t>              & moviePartHoleNumbers,
        vector<MovieAlnIndexLookupTable> & lookupTables) {

    lookupTables.resize(movieAlnIndex.size());
    int movieAlignmentIndex = 0;
    for (movieAlignmentIndex = 0; movieAlignmentIndex < movieAlnIndex.size(); movieAlignmentIndex++) {
        BuildLookupTable(movieAlignmentIndex,
            cmpFile, 
            baseFile,
            usePulseFile,
            pulseFile,
            cmpReader, 
            movieAlnIndex,
            toFrom,
            moviePartHoleNumbers, 
            lookupTables[movieAlignmentIndex]);
    }

    //
    // Load entire Basecall from pls/bas to memory, and
    // check whether aligned sequences in cmp.h5 matches 
    // sequences in pls/bas or not
    //
    hdfBasReader.ReadField(baseFile, "Basecall");

    // 
    // For each alignment, do sanity check and 
    // cache aligned sequence in MovieAlnIndexLookupTable
    //
    for (movieAlignmentIndex = 0; movieAlignmentIndex < movieAlnIndex.size(); movieAlignmentIndex++) {
        MovieAlnIndexLookupTable & table = lookupTables[movieAlignmentIndex];
        if (table.skip) continue;
        //
        // Get aligned sequence for this alignment from cmpFile 
        //
        string alignedSequence = GetAlignedSequenceFromCmpFile(cmpReader, table);
        
        // Save the aligned sequence in the table
        table.alignedSequence = alignedSequence;

        RemoveGaps(alignedSequence, alignedSequence);

        //
        // Get sequence for this alignment from baseFile
        //
        Nucleotide * seq = new Nucleotide[table.readLength]; 
        baseFile.CopyArray(baseFile.baseCalls, table.readStart, table.readLength, seq);
        
        string readSequence;
        readSequence.resize(table.queryEnd - table.queryStart);
        copy((char*) (seq + table.queryStart), 
             (char*) (seq + table.queryEnd), 
             readSequence.begin());
        delete seq;
        
        //
        // Do a sanity check to make sure the pulses and the alignment
        // make sense.  The main check is to see if the query sequence
        // in the alignment is the same as the query sequence in the
        // read.
        //
        if (alignedSequence.size() != readSequence.size() or alignedSequence != readSequence) {
            cout << "ERROR, the query sequence does not match the aligned query sequence." << endl
                 << "HoleNumber: "       << cmpFile.alnInfo.alignments[table.alignmentIndex].GetHoleNumber()
                 << ", MovieName: "      << baseFile.GetMovieName()
                 << ", ReadIndex: "      << table.readIndex
                 << ", qStart: "         << table.queryStart 
                 << ", qEnd: "           << table.queryEnd << endl
                 << "Aligned sequence: " << endl
                 << alignedSequence      << endl
                 << "Original sequence: "<< endl
                 << readSequence         << endl;
            exit(1);
        }
    }

    hdfBasReader.ClearField(baseFile, "Basecall");
}

// Given a vector of lookupTables in which items with the same 
// refGroupIndex and readGroupIndex are grouped, find index boundaries
// of each group and save these boundaries to groupedLookupTablesIndexPairs
// The index boundary of each group consists of:
//   1, index (0 based, inclusive) of the very first item of a group
//   2, index (0 based, exclusive) of the very last item of a group 
//
// Assume that lookupTables satisfy the following criteria.
//   1, items are already grouped by refGroupIndex and readGroupIndex 
//   2, items which have the same alnGroupIndex, should have 
//      the same refGroupIndex and readGroupIndex 
// Note that: 
//   1, alnGroupIndex represents index of AlnGroupID, (i.e. dataset
//      /AlnInfo/AlnIndex column 1);
//      refGroupIndex represents index of RefGroupID, (i.e. dataset
//      /AlnInfo/AlnIndex column 3);
//      readGroupIndex represents index of an experiment group within
//      a refGroup (e.g. if a refGroup /ref0001 contains two experiment
//      groups /ref0001/movie1 and /ref0001/movie2, then readGroupIndex
//      for these two groups are 0 and 1.).
//   2, within each grouped item, offsetBegin may not begin from 0, 
//      and offsets may not be continugous. 
//
void GroupLookupTables(
        vector<MovieAlnIndexLookupTable>  & lookupTables,
        vector<pair<UInt, UInt> > & groupedLookupTablesIndexPairs) {

    vector<pair<UInt, UInt> > refGroupIndexReadGroupIndexPairs;
    UInt movieAlignmentIndex = 0;
    UInt preRefGroupIndex    = 0;
    UInt preReadGroupIndex   = 0;
    UInt pairFirst           = 0;
    bool isVeryFirstGroup    = true;

    for (movieAlignmentIndex = 0; movieAlignmentIndex < lookupTables.size(); movieAlignmentIndex++) {
        MovieAlnIndexLookupTable & lookupTable = lookupTables[movieAlignmentIndex];
        
        if (isVeryFirstGroup or
            (lookupTable.refGroupIndex  != preRefGroupIndex or
             lookupTable.readGroupIndex != preReadGroupIndex)) {
            // Find a new group 
            if (isVeryFirstGroup) {
                // This is the very first group
                isVeryFirstGroup = false;
            } else 
            if (lookupTable.refGroupIndex  == preRefGroupIndex && 
                lookupTable.readGroupIndex != preReadGroupIndex) {
                // Assumption (1) has been violated 
                cout << "ERROR! lookupTables should have been sorted by reference"
                     << "group index and read group index." << endl;
                exit(1);
            } else {
                // Find the first lookupTable of a new group, save indices of [first and last)
                // lookupTables of the last group.
                groupedLookupTablesIndexPairs.push_back(pair<UInt,UInt> (pairFirst, movieAlignmentIndex));
                // Save refGroupIndex and readGroupIndex of the last group
                pair<UInt,UInt> refGroupIndexReadGroupIndexPair(preRefGroupIndex, preReadGroupIndex);
                refGroupIndexReadGroupIndexPairs.push_back(refGroupIndexReadGroupIndexPair);
            }

            // Store index of the first lookupTable of the new group in lookupTables
            pairFirst = movieAlignmentIndex; 
            // Store refGroupIndex and readGroupIndex of the new group
            preRefGroupIndex  = lookupTable.refGroupIndex;
            preReadGroupIndex = lookupTable.readGroupIndex;
        } 
    }
    if (not isVeryFirstGroup) {
        // Save indices of [first and last) lookupTables of the very last group
        groupedLookupTablesIndexPairs.push_back(pair<UInt,UInt> (pairFirst, movieAlignmentIndex));
        // Save refGroupIndex and readGroupIndex of the very last group
        pair<UInt,UInt> refGroupIndexReadGroupIndexPair(preRefGroupIndex, preReadGroupIndex);
        refGroupIndexReadGroupIndexPairs.push_back(refGroupIndexReadGroupIndexPair);
    } // Do nothing, if no lookupTable exists


    // Double check all assumptions are met
    for (int i = 0; i < refGroupIndexReadGroupIndexPairs.size(); i++) {
        for (int j = i+1; j < refGroupIndexReadGroupIndexPairs.size(); j++) {
            // Assure that assumption (1) is met. If this assertion fails, 
            // then alignments in the input cmp.h5 are not grouped by
            // reference. Check /AlnInfo/AlnIndex dataset column 3.
            assert(refGroupIndexReadGroupIndexPairs[i] != refGroupIndexReadGroupIndexPairs[j]);
        }
    }
    assert(groupedLookupTablesIndexPairs.size() == refGroupIndexReadGroupIndexPairs.size());
    int i ;
    for (i = 0; i < groupedLookupTablesIndexPairs.size(); i++) {
        UInt firstIndex     = groupedLookupTablesIndexPairs[i].first;
        UInt lastIndex      = groupedLookupTablesIndexPairs[i].second;
        UInt refGroupIndex  = refGroupIndexReadGroupIndexPairs[i].first;
        UInt readGroupIndex = refGroupIndexReadGroupIndexPairs[i].second;
        for(UInt index = firstIndex; index < lastIndex; index++) {
            assert(lookupTables[index].refGroupIndex  == refGroupIndex);
            assert(lookupTables[index].readGroupIndex == readGroupIndex);
        }
    }
}


//
// Read all required fields for computing the specified metric into memory, 
// unless the fields have been cached.
//
void CacheRequiredFieldsForMetric(
        BaseFile                    & baseFile,
        PulseFile                   & pulseFile,
        HDFBasReader                & hdfBasReader,
        HDFPlsReader                & hdfPlsReader,
        HDFCCSReader<SMRTSequence>  & hdfCcsReader,
        const bool                  & useBaseFile,
        const bool                  & usePulseFile,
        const bool                  & useCcsOnly,
        vector<Field>               & cachedFields,
        const string                & curMetric) {

    vector<Field> fieldsToBeUsed;
    bool canBeComputed = CanThisMetricBeComputed( 
        curMetric, hdfBasReader, hdfPlsReader,
        useBaseFile, usePulseFile, fieldsToBeUsed); 
    assert(canBeComputed);

    // Cache all required fields 
    for (int i = 0; i < fieldsToBeUsed.size(); i++) {
        bool isFieldCached = false;
        for (int j = 0; j < cachedFields.size(); j++) {
            if (fieldsToBeUsed[i] == cachedFields[j]) {
                isFieldCached = true;
                break;
            }
        }
        if (isFieldCached) {
            continue;
        }
        string    & curField = fieldsToBeUsed[i].name;
        FieldType & fieldType= fieldsToBeUsed[i].type;

        if (fieldType == BasField and useBaseFile
            and hdfBasReader.FieldIsIncluded(curField)
            and hdfBasReader.includedFields[curField]) {
            hdfBasReader.ReadField(baseFile, curField);
            cachedFields.push_back(fieldsToBeUsed[i]);
        } else 
        if (fieldType == PlsField and usePulseFile 
            and hdfPlsReader.FieldIsIncluded(curField)
            and hdfPlsReader.includedFields[curField]) {
            hdfPlsReader.ReadField(pulseFile, curField);
            cachedFields.push_back(fieldsToBeUsed[i]);
        }
    }
}

// 
// Clear cached fields unless they are also required for computing 
// the next metric.
//
void ClearCachedFields(
        BaseFile                   & baseFile,
        PulseFile                  & pulseFile,
        HDFBasReader               & hdfBasReader,
        HDFPlsReader               & hdfPlsReader,
        HDFCCSReader<SMRTSequence> & hdfCcsReader,
        const bool                 & useBaseFile,
        const bool                 & usePulseFile,
        const bool                 & useCcsOnly,
        vector<Field>              & cachedFields,
        const string               & curMetric,
        const string               & nextMetric) {

 
    vector<Field> nextRequiredFields;
    if (nextMetric != "") {
        bool canBeComputed = CanThisMetricBeComputed( 
            nextMetric, hdfBasReader, hdfPlsReader,
            useBaseFile, usePulseFile, nextRequiredFields); 
        assert(canBeComputed);
    }
    for (int i = 0; i < cachedFields.size(); i++) {
        bool isRequiredForNextMetric = false;
        for (int j = 0; j < nextRequiredFields.size(); j++) {
            if (cachedFields[i] == nextRequiredFields[j]) {
                isRequiredForNextMetric = true;
                break;
            }
        }
        if (isRequiredForNextMetric) {
            continue;
        }
        string    & curField = cachedFields[i].name;
        FieldType & fieldType= cachedFields[i].type;

        if (fieldType == BasField and useBaseFile  and 
            hdfBasReader.FieldIsIncluded(curField) and
            hdfBasReader.includedFields[curField]) {
            hdfBasReader.ClearField(baseFile, curField);
            // Remove it from cachedFields
            cachedFields.erase(cachedFields.begin()+i);
            i--;
        } else 
        if (fieldType == PlsField and usePulseFile and 
            hdfPlsReader.FieldIsIncluded(curField) and
            hdfPlsReader.includedFields[curField]) {
            if (curField == "NumEvent") {
                // Always keep NumEvent 
                continue;
            }
            hdfPlsReader.ClearField(pulseFile, curField);
            // Remove it from cachedFields
            cachedFields.erase(cachedFields.begin()+i);
            i--;
        }
    }        
}

// Compute StartFrame from BaseCalls only. 
// Return true if succeed, false otherwise.
bool ComputeStartFrameFromBase(
        BaseFile           & baseFile,
        HDFBasReader       & hdfBasReader,
        const bool         & useBaseFile,
        MovieAlnIndexLookupTable & lookupTable,
        vector<UInt>       & newStartFrame) {
    newStartFrame.resize(lookupTable.readLength);
    if (useBaseFile and hdfBasReader.FieldIsIncluded("PreBaseFrames") 
        and hdfBasReader.includedFields["PreBaseFrames"]
        and baseFile.preBaseFrames.size() > 0) {
        // baseFile.preBaseFrame data type = uint16
        // startFrame data type = uint32
        for (int i = 0; i < lookupTable.readLength; i++) {
            newStartFrame[i] = baseFile.preBaseFrames[lookupTable.readStart+i];
        }
        for (int i = 0; i < lookupTable.readLength-1; i++) {
            newStartFrame[i+1] += baseFile.basWidthInFrames[lookupTable.readStart+i];
        }
        partial_sum(&newStartFrame[0], &newStartFrame[lookupTable.readLength], &newStartFrame[0]);
        return true;
    } 
    return false;
}

// Compute StartFrame from PulseCalls only.
// Return true if succeed, false otherwise.
bool ComputeStartFrameFromPulse(
        PulseFile          & pulseFile,
        HDFPlsReader       & hdfPlsReader,
        const bool         & usePulseFile,
        MovieAlnIndexLookupTable & lookupTable,
        vector<int>        & baseToPulseIndexMap,
        vector<UInt>       & newStartFrame) {
    newStartFrame.resize(lookupTable.readLength);
    if (usePulseFile) {
        assert(pulseFile.startFrame.size() > 0);
        hdfPlsReader.CopyFieldAt(pulseFile, "StartFrame", lookupTable.plsReadIndex,
                &baseToPulseIndexMap[0], &newStartFrame[0],
                lookupTable.readLength);
        return true;
    }
    return false;
}


// Compute StartFrame from either (1) BaseCalls or (2) PulseCalls.
//    (1) Uses baseFile.preBaseFrames and baseFile.basWidthInFrames
//    (2) Uses pulseFile.startFrame
// In theory, the generated results using both methods should 
// be exactly the same. However, they can be different in practice
// because PreBaseFrames is of data type uint_16, while its
// value can exceed maximum uint_16 (65535). 
// When possible, always use PulseCalls. 
void ComputeStartFrame(
        BaseFile           & baseFile,
        PulseFile          & pulseFile,
        HDFBasReader       & hdfBasReader,
        HDFPlsReader       & hdfPlsReader,
        bool                 useBaseFile,
        bool                 usePulseFile,
        MovieAlnIndexLookupTable & lookupTable,
        vector<int>        & baseToPulseIndexMap,
        vector<UInt>       & newStartFrame) {

    if (!ComputeStartFrameFromPulse(pulseFile, hdfPlsReader, usePulseFile,
                lookupTable, baseToPulseIndexMap, newStartFrame)) {
        if (!ComputeStartFrameFromBase(baseFile, hdfBasReader, useBaseFile,
                lookupTable, newStartFrame)) {
            cout << "ERROR! There is insufficient data to compute metric: StartFrame." << endl;
            exit(1);
        }
    }
}

//
// Compute and write an entire metric to cmp.h5.
// Assume that all required fields have been loaded.
//
void WriteMetric(
        CmpFile                          & cmpFile,
        BaseFile                         & baseFile,
        PulseFile                        & pulseFile,
        HDFCmpFile<CmpAlignment>         & cmpReader,
        HDFBasReader                     & hdfBasReader,
        HDFPlsReader                     & hdfPlsReader,
        HDFCCSReader<SMRTSequence>       & hdfCcsReader,
        const bool                       & useBaseFile,
        const bool                       & usePulseFile,
        const bool                       & useCcsOnly,
        vector<MovieAlnIndexLookupTable> & lookupTables,
        vector<pair<UInt, UInt> >        & groupedLookupTablesIndexPairs,
        const string                     & curMetric ) {

    int movieAlignmentIndex = 0;
    for (int index = 0; index < groupedLookupTablesIndexPairs.size(); index++) {
        // Group[index] contains all items in lookupTables[firstIndex...lastIndex)
        UInt firstIndex = groupedLookupTablesIndexPairs[index].first;
        UInt lastIndex  = groupedLookupTablesIndexPairs[index].second;

        assert(lookupTables.size() > firstIndex);
        UInt refGroupIndex  = lookupTables[firstIndex].refGroupIndex;
        UInt readGroupIndex = lookupTables[firstIndex].readGroupIndex;
        // Obtain alignment array length from *.cmp.h5/refGroup/readGroup/AlnArray.
        HDFCmpExperimentGroup* expGroup = cmpReader.refAlignGroups[refGroupIndex]->readGroups[readGroupIndex];
        UInt alnArrayLength = expGroup->alignmentArray.size();

        //
        // Compute any necessary data fields.  These usually involve
        // using differences of pulse indices, pulse widths, etc..
        // Missing fields are stored as 0's. 
        //
        vector<UInt>     startTimeOffsetMetric;
        // pulseIndex's data type is uint16 in ICD, 
        // but I have seen it defined as uint32 in a bas file.
        vector<UInt>     pulseMetric; 
        vector<UChar>    qvMetric;
        vector<HalfWord> frameRateMetric;
        vector<UInt>     timeMetric;
        vector<char>     tagMetric; 
        vector<float>    floatMetric;
       
        /*
        if (curMetric == "StartTimeOffset") {
            startTimeOffsetMetric.resize(alnNum);
            HDFArray<UInt> * data = (HDFArray<UInt>*) expGroup->fields[curMetric];
            if (data->IsInitialized()) {
                assert(data->size() == alnNum);
                data->UpdateH5Dataspace();
                data->Read(0, alnNum-1, &StartTimeOffsetMetric[0]);
            } else {
                data->Initialize(expGroup->experimentGroup, curMetric);
                fill(startTimeOffsetMetric.begin(), startTimeOffsetMetric.end(), );
            }

        } else */
        if (curMetric == "QualityValue" || curMetric == "InsertionQV" ||
            curMetric == "DeletionQV"   || curMetric == "MergeQV"     ||
            curMetric == "SubstitutionQV") {
            qvMetric.resize(alnArrayLength);
            HDFArray<UChar> * data = (HDFArray<UChar>*) expGroup->fields[curMetric];
            if (data->IsInitialized()) {
                assert(data->size() == alnArrayLength);
                data->UpdateH5Dataspace();
                data->Read(0, alnArrayLength-1, &qvMetric[0]);
            } else {
                data->Initialize(expGroup->experimentGroup, curMetric);
                //fill(qvMetric.begin(), qvMetric.end(), missingQualityValue);
            }

        } else if (curMetric == "ClassifierQV" || curMetric == "pkmid" ) {
            // Note that data type of pkmid=midSignal, which is uint_8 in bas/pls files,
            // has been changed to float in cmp.h5. Why?
            floatMetric.resize(alnArrayLength);
            HDFArray<float> * data = (HDFArray<float>*) expGroup->fields[curMetric];
            if (data->IsInitialized()) {
                assert(data->size() == alnArrayLength);
                data->UpdateH5Dataspace();
                data->Read(0, alnArrayLength-1, &floatMetric[0]);
            } else {
                data->Initialize(expGroup->experimentGroup, curMetric);
                //fill(floatMetric.begin(), floatMetric.end(), NaN);
            }

        } else if (curMetric == "PulseIndex"     ) {
            pulseMetric.resize(alnArrayLength);
            HDFArray<UInt> * data = (HDFArray<UInt>*) expGroup->fields[curMetric];
            if (data->IsInitialized()) {
                assert(data->size() == alnArrayLength);
                data->UpdateH5Dataspace();
                data->Read(0, alnArrayLength-1, &pulseMetric[0]);
            } else {
                data->Initialize(expGroup->experimentGroup, curMetric);
                //fill(pulseMetric.begin(), pulseMetric.end(), 0);
            }

        } else if (curMetric == "DeletionTag"  || curMetric == "SubstitutionTag") {
            tagMetric.resize(alnArrayLength);
            HDFArray<char> * data = (HDFArray<char>*) expGroup->fields[curMetric];
            if (data->IsInitialized()) {
                assert(data->size() == alnArrayLength);
                data->UpdateH5Dataspace();
                data->Read(0, alnArrayLength-1, &tagMetric[0]);
            } else {
                data->Initialize(expGroup->experimentGroup, curMetric);
                //fill(tagMetric.begin(), tagMetric.end(), '-'); 
            }

        } else if (curMetric == "StartFrame"   || curMetric == "StartFrameBase" ||
                   curMetric == "StartFramePulse") {
            timeMetric.resize(alnArrayLength);
            HDFArray<UInt> * data = (HDFArray<UInt>*) expGroup->fields[curMetric];
            if (data->IsInitialized()) {
                assert(data->size() == alnArrayLength);
                data->UpdateH5Dataspace();
                data->Read(0, alnArrayLength-1, &timeMetric[0]);
            } else {
                data->Initialize(expGroup->experimentGroup, curMetric);
                //fill(timeMetric.begin(), timeMetric.end(), missingPulseIndex);
            }

        } else if (curMetric == "PulseWidth"   || curMetric == "PreBaseFrames" ||
                   curMetric == "WidthInFrames"|| curMetric == "IPD"           ||
                   curMetric == "Light") {
            frameRateMetric.resize(alnArrayLength);
            HDFArray<HalfWord> * data = (HDFArray<HalfWord>*) expGroup->fields[curMetric];
            if (data->IsInitialized()) {
                assert(data->size() == alnArrayLength);
                data->UpdateH5Dataspace();
                data->Read(0, alnArrayLength-1, &frameRateMetric[0]);
            } else {
                data->Initialize(expGroup->experimentGroup, curMetric);
                //fill(frameRateMetric.begin(), frameRateMetric.end(), missingFrameRateValue);
            }
        } else {
            cout << "ERROR, metric " << curMetric << " is not supported." << endl;
            exit(1);
        }

        for (movieAlignmentIndex = firstIndex; movieAlignmentIndex < lastIndex; movieAlignmentIndex++) {
            MovieAlnIndexLookupTable & lookupTable   = lookupTables[movieAlignmentIndex];
            if (lookupTable.skip) continue;

            const UInt alignedSequenceLength         = lookupTable.offsetEnd - lookupTable.offsetBegin; 
            const UInt ungappedAlignedSequenceLength = lookupTable.queryEnd  - lookupTable.queryStart;
            const UInt   & readIndex                 = lookupTable.readIndex;
            const UInt   & plsReadIndex              = lookupTable.plsReadIndex;
            const UInt   & readStart                 = lookupTable.readStart;
            const UInt   & readLength                = lookupTable.readLength;
            const UInt   & queryStart                = lookupTable.queryStart;
            const UInt   & offsetBegin               = lookupTable.offsetBegin;
            const UInt   & offsetEnd                 = lookupTable.offsetEnd;
            assert (offsetEnd <= alnArrayLength); 
            assert (offsetBegin+alignedSequenceLength <= alnArrayLength);

            // Condense gaps and get ungapped aligned sequence.
            string ungappedAlignedSequence = lookupTable.alignedSequence;
            RemoveGaps(ungappedAlignedSequence, ungappedAlignedSequence);

            vector<int> baseToAlignmentMap;
            // Map bases in the aligned sequence to their positions in the alignment.
            CreateSequenceToAlignmentMap(lookupTable.alignedSequence, baseToAlignmentMap);

            vector<int> baseToPulseIndexMap;
            if (usePulseFile && IsPulseMetric(curMetric)) {
                // Map bases in the read to pulse indices. 
                MapBaseToPulseIndex(baseFile, pulseFile, lookupTable, baseToPulseIndexMap);
            }

            UInt i;
            if (curMetric == "QualityValue") {
                assert(baseFile.qualityValues.size() > 0 && 
                       baseFile.qualityValues.size() >= readStart + readLength);
                fill(&qvMetric[offsetBegin], &qvMetric[offsetEnd], missingPulseIndex);
                for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                    // cap quality value
                    qvMetric[offsetBegin+baseToAlignmentMap[i]] = min(maxQualityValue, baseFile.qualityValues[readStart+queryStart+i]);
                }
                qvMetric[offsetBegin+alignedSequenceLength] = 0;

            } else if (curMetric == "InsertionQV") {
                assert(baseFile.insertionQV.size() > 0 && 
                       baseFile.insertionQV.size() >= readStart + readLength);
                fill(&qvMetric[offsetBegin], &qvMetric[offsetEnd], missingPulseIndex);
                for (i = 0; i < ungappedAlignedSequenceLength; i++) {
                    // cap quality value
                    qvMetric[offsetBegin+baseToAlignmentMap[i]] = min(maxQualityValue, baseFile.insertionQV[readStart+queryStart+i]);
                }
                qvMetric[offsetBegin+alignedSequenceLength] = 0;
            
            } else if (curMetric == "MergeQV") {
                assert(baseFile.mergeQV.size() > 0 && 
                       baseFile.mergeQV.size() >= readStart + readLength);
                fill(&qvMetric[offsetBegin], &qvMetric[offsetEnd], missingPulseIndex);
                for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                    // cap quality value
                    qvMetric[offsetBegin+baseToAlignmentMap[i]] = min(maxQualityValue, baseFile.mergeQV[readStart+queryStart+i]);
                }
                qvMetric[offsetBegin+alignedSequenceLength] = 0;

            } else if (curMetric == "DeletionQV") {
                assert(baseFile.deletionQV.size() > 0 && 
                       baseFile.deletionQV.size() >= readStart + readLength);
                fill(&qvMetric[offsetBegin], &qvMetric[offsetEnd], missingPulseIndex);
                for (i = 0; i < ungappedAlignedSequenceLength; i++) {
                    // cap quality value
                    qvMetric[offsetBegin+baseToAlignmentMap[i]] = min(maxQualityValue, baseFile.deletionQV[readStart+queryStart+i]);
                }
                qvMetric[offsetBegin+alignedSequenceLength] = 0;

            } else if (curMetric == "DeletionTag") {
                assert(baseFile.deletionTag.size() > 0 && 
                       baseFile.deletionTag.size() >= readStart + readLength);
                fill(&tagMetric[offsetBegin], &tagMetric[offsetEnd], '-'); 
                for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                    assert(offsetBegin+baseToAlignmentMap[i] < tagMetric.size());
                    tagMetric[offsetBegin+baseToAlignmentMap[i]] = baseFile.deletionTag[readStart+queryStart+i];
                }
                tagMetric[offsetBegin+alignedSequenceLength] = 0;

            } else if (curMetric == "PulseIndex") {
                assert(baseFile.pulseIndex.size() > 0 &&
                       baseFile.pulseIndex.size() >= readStart + readLength);
                fill(&pulseMetric[offsetBegin], &pulseMetric[offsetEnd], 0);
                for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                    pulseMetric[offsetBegin+baseToAlignmentMap[i]] = baseFile.pulseIndex[readStart+queryStart+i];
                }
                pulseMetric[offsetBegin+alignedSequenceLength] = 0;

            } else if (curMetric == "SubstitutionTag") {
                assert(baseFile.substitutionTag.size() > 0 &&
                       baseFile.substitutionTag.size() >= readStart + readLength);
                fill(&tagMetric[offsetBegin], &tagMetric[offsetEnd], '-'); 
                for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                    tagMetric[offsetBegin+baseToAlignmentMap[i]] = baseFile.substitutionTag[readStart+queryStart+i];
                }
                tagMetric[offsetBegin+alignedSequenceLength] = 0;

            } else if (curMetric == "SubstitutionQV") {
                assert(baseFile.substitutionQV.size() > 0 &&
                       baseFile.substitutionQV.size() >= readStart + readLength);
                fill(&qvMetric[offsetBegin], &qvMetric[offsetEnd], missingPulseIndex);
                for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                    qvMetric[offsetBegin+baseToAlignmentMap[i]] = min(maxQualityValue, baseFile.substitutionQV[readStart+queryStart+i]);
                }
                qvMetric[offsetBegin+alignedSequenceLength] = 0;

            } else if (curMetric == "ClassifierQV") {
                assert(pulseFile.classifierQV.size() > 0 && 
                       pulseFile.classifierQV.size() >= readStart + readLength);
                vector<float> newClassifierQV; 
                newClassifierQV.resize(ungappedAlignedSequenceLength);
                // For the data used for this table, it is possible to simply
                // reference the data for the bas file,  but for the pls file,
                // it is necessary to copy since there is a packing of data.
                hdfPlsReader.CopyFieldAt(pulseFile, "ClassifierQV", plsReadIndex, 
                        &baseToPulseIndexMap[queryStart], &newClassifierQV[0], 
                        ungappedAlignedSequenceLength);
                
                fill(&floatMetric[offsetBegin], &floatMetric[offsetEnd], NaN);
                for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                    floatMetric[offsetBegin+baseToAlignmentMap[i]] = newClassifierQV[i];
                }
                floatMetric[offsetBegin+alignedSequenceLength] = 0;
            
/*            } else if (curMetric == "StartTimeOffset") {
                // StartTimeOffset is a subset of StartFrame.
                vector<UInt> newStartFrame;
                ComputeStartFrame(baseFile, pulseFile, hdfBasReader, hdfPlsReader,
                                  useBaseFile, usePulseFile, lookupTable, 
                                  baseToPulseIndexMap, newStartFrame);
        
                startTimeOffsetMetric[offsetBegin] = newStartFrame[queryStart];
*/
           } else if (curMetric == "StartFrame") {
                vector<UInt> newStartFrame;
                ComputeStartFrame(baseFile, pulseFile, hdfBasReader, hdfPlsReader,
                                  useBaseFile, usePulseFile, lookupTable, 
                                  baseToPulseIndexMap, newStartFrame);
                fill(&timeMetric[offsetBegin], &timeMetric[offsetEnd], missingPulseIndex);
                for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                    timeMetric[offsetBegin+baseToAlignmentMap[i]] = newStartFrame[queryStart+i];
                }
                timeMetric[offsetBegin+alignedSequenceLength] = 0;

            } else if (curMetric == "StartFrameBase") {
                // Sneaky metric, compute StartFrame from BaseCalls only. 
                vector<UInt> newStartFrame;
                ComputeStartFrameFromBase(baseFile, hdfBasReader, useBaseFile, 
                                  lookupTable, newStartFrame);
                fill(&timeMetric[offsetBegin], &timeMetric[offsetEnd], missingPulseIndex);
                for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                    timeMetric[offsetBegin+baseToAlignmentMap[i]] = newStartFrame[queryStart+i];
                }
                timeMetric[offsetBegin+alignedSequenceLength] = 0;

            } else if (curMetric == "StartFramePulse") {
                // Sneaky metric, compute StartFrame from PulseCalls only. 
                vector<UInt> newStartFrame;
                ComputeStartFrameFromPulse(pulseFile, hdfPlsReader, usePulseFile, 
                                  lookupTable, baseToPulseIndexMap, newStartFrame);
                fill(&timeMetric[offsetBegin], &timeMetric[offsetEnd], missingPulseIndex);
                for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                    timeMetric[offsetBegin+baseToAlignmentMap[i]] = newStartFrame[queryStart+i];
                }
                timeMetric[offsetBegin+alignedSequenceLength] = 0;

            } else if (curMetric == "PreBaseFrames") {
                // Directly load baseFile.PreBaseFrames.
                // DON'T compute it from PulseCalls even if you can.
                assert(baseFile.preBaseFrames.size() > 0 &&
                       baseFile.preBaseFrames.size() >= readStart + readLength); 
                fill(&frameRateMetric[offsetBegin], &frameRateMetric[offsetEnd], missingFrameRateValue);
                for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                    frameRateMetric[offsetBegin+baseToAlignmentMap[i]] = baseFile.preBaseFrames[readStart+queryStart+i];
                }
                frameRateMetric[offsetBegin+alignedSequenceLength] = 0;

            } else if (curMetric == "WidthInFrames" || curMetric == "PulseWidth") {
                // For legacy reasons, it's possible the width in frames is
                // stored in the bas file. If this is the case, use the width
                // in frames there.  Otherwise, use the width in frames stored
                // in the pls file.
                vector<uint16_t> newWidthInFrames;
                newWidthInFrames.resize(ungappedAlignedSequenceLength);
                if (usePulseFile) {
                    hdfPlsReader.CopyFieldAt(pulseFile, "WidthInFrames", plsReadIndex,
                            &baseToPulseIndexMap[queryStart], &newWidthInFrames[0], 
                            ungappedAlignedSequenceLength);
                } else
                if (useBaseFile) {
                    // basWidthInFrames data type uint16
                    copy(&baseFile.basWidthInFrames[readStart+queryStart], 
                         &baseFile.basWidthInFrames[readStart+queryStart+ungappedAlignedSequenceLength],
                         &newWidthInFrames[0]);
                } 
                
                fill(&frameRateMetric[offsetBegin], &frameRateMetric[offsetEnd], missingFrameRateValue);
                for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                    frameRateMetric[offsetBegin+baseToAlignmentMap[i]] = newWidthInFrames[i];
                }
                frameRateMetric[offsetBegin+alignedSequenceLength] = 0;

            } else if (curMetric == "pkmid") {
                // pkmid in cmp.h5 is MidSignal in pls.h5, but
                // data type of MidSignal is uint16 in pls files, 
                // data type of pkmid is float in cmp files.
                assert(usePulseFile);
                vector<HalfWord> newMidSignal; 
                newMidSignal.resize(ungappedAlignedSequenceLength);
                hdfPlsReader.CopyFieldAt(pulseFile, "MidSignal", plsReadIndex,
                        &baseToPulseIndexMap[queryStart], &newMidSignal[0],
                        ungappedAlignedSequenceLength, ungappedAlignedSequence);

                fill(&floatMetric[offsetBegin], &floatMetric[offsetEnd], NaN);
                for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                    floatMetric[offsetBegin+baseToAlignmentMap[i]] = newMidSignal[i];
                }
                floatMetric[offsetBegin+alignedSequenceLength] = 0;

            } else if (curMetric == "IPD") {
                fill(&frameRateMetric[offsetBegin], &frameRateMetric[offsetEnd], missingFrameRateValue);

                // IPD can be either (1) copied from baseFile.preBaseFrames 
                // or (2) computed from pulseFile.StartFrame and pulseFile.WidthInFrames
                // Always use method (2) when possible as it is more accurate.
                if (usePulseFile) {
                    // Need to read StartFrame & WidthInFrames for the entire read,
                    // not only for a subset of bases in the alignment
                    assert(pulseFile.startFrame.size() > 0);
                    assert(pulseFile.plsWidthInFrames.size() > 0);
                    vector<UInt> newStartFrame;
                    newStartFrame.resize(readLength);
                    hdfPlsReader.CopyFieldAt(pulseFile, "StartFrame", plsReadIndex,
                         &baseToPulseIndexMap[0], &newStartFrame[0], readLength);

                    vector<uint16_t> newWidthInFrames;
                    newWidthInFrames.resize(readLength);
                    hdfPlsReader.CopyFieldAt(pulseFile, "WidthInFrames", plsReadIndex,
                         &baseToPulseIndexMap[0], &newWidthInFrames[0], readLength);

                    for (i = 0; i < ungappedAlignedSequenceLength; i++) {
                        // The IPD is undefined for the first base in a read.
                        if (queryStart == 0 and i == 0) {
                            frameRateMetric[offsetBegin+baseToAlignmentMap[i]] = 0;
                        } else {
                            frameRateMetric[offsetBegin+baseToAlignmentMap[i]] = newStartFrame[queryStart+i]  
                                - newStartFrame[i+queryStart-1] - newWidthInFrames[i+queryStart-1];
                        }
                    }
                } else 
                if (useBaseFile) {
                    assert(baseFile.preBaseFrames.size() > 0);
                    assert(baseFile.preBaseFrames.size() >= readStart + readLength); 

                    for (i = 0; i < ungappedAlignedSequenceLength; i++) {
                        frameRateMetric[offsetBegin+baseToAlignmentMap[i]] = 
                            baseFile.preBaseFrames[readStart+queryStart+i];
                    }
                }
                frameRateMetric[offsetBegin+alignedSequenceLength] = 0;

            } else if (curMetric == "Light") {
                // Light can be computed from pulseFile.meanSignal and 
                // pulseFile.plsWidthInFrames. Might have been deprecated.
                assert(usePulseFile);
                fill(&frameRateMetric[offsetBegin], &frameRateMetric[offsetEnd], missingFrameRateValue);

                vector<uint16_t> newMeanSignal; 
                newMeanSignal.resize(ungappedAlignedSequenceLength);
                hdfPlsReader.CopyFieldAt(pulseFile, "MeanSignal", plsReadIndex,
                        &baseToPulseIndexMap[queryStart], &newMeanSignal[0], 
                        ungappedAlignedSequenceLength, ungappedAlignedSequence);

                vector<uint16_t> newWidthInFrames;
                newWidthInFrames.resize(ungappedAlignedSequenceLength);
                hdfPlsReader.CopyFieldAt(pulseFile, "WidthInFrames", plsReadIndex,
                        &baseToPulseIndexMap[queryStart], &newWidthInFrames[0], 
                        ungappedAlignedSequenceLength);

                for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                    frameRateMetric[offsetBegin+baseToAlignmentMap[i]] =  newMeanSignal[i] * newWidthInFrames[i];
                }
                frameRateMetric[offsetBegin+alignedSequenceLength] = 0;

            } else { 
                cout << "ERROR, unknown metric " << curMetric << endl;
                exit(1);
            } 
        }

        // Write the computed metric to cmp.h5.
        /*if (curMetric == "StartTimeOffset") {
            expGroup->startTimeOffset.WriteToPos(&startTimeOffsetMetric[0], startTimeOffsetMetric.size(), 0);
        } else */
        if (curMetric == "QualityValue" || curMetric == "InsertionQV" ||
                   curMetric == "DeletionQV"   || curMetric == "MergeQV"     ||
                   curMetric == "SubstitutionQV") {
            HDFArray<UChar> * data = (HDFArray<UChar> *) expGroup->fields[curMetric];
            data->WriteToPos(&qvMetric[0], qvMetric.size(), 0);
 
        } else if (curMetric == "ClassifierQV" || curMetric == "pkmid" ) {
            HDFArray<float> * data = (HDFArray<float> *) expGroup->fields[curMetric];
            data->WriteToPos(&floatMetric[0], floatMetric.size(), 0);

        } else if (curMetric == "PulseIndex") {
            HDFArray<UInt> * data = (HDFArray<UInt> *) expGroup->fields[curMetric];
            data->WriteToPos(&pulseMetric[0], pulseMetric.size(), 0);
        
        } else if (curMetric == "DeletionTag"  || curMetric == "SubstitutionTag") {
            HDFArray<char> * data = (HDFArray<char> *) expGroup->fields[curMetric];
            data->WriteToPos(&tagMetric[0], tagMetric.size(), 0);
      
        } else if (curMetric == "StartFrame"   || curMetric == "StartFrameBase"||
                   curMetric == "StartFramePulse") {
            HDFArray<UInt> * data = (HDFArray<UInt>*) expGroup->fields[curMetric];
            data->WriteToPos(&timeMetric[0], timeMetric.size(), 0);

        } else if (curMetric == "PulseWidth"   || curMetric == "PreBaseFrames" ||
                   curMetric == "WidthInFrames"|| curMetric == "IPD"           ||
                   curMetric == "Light") {
            HDFArray<HalfWord> * data = (HDFArray<HalfWord>*) expGroup->fields[curMetric];
            data->WriteToPos(&frameRateMetric[0], frameRateMetric.size(), 0);

        } else {
            cout << "ERROR, unknown metric " << curMetric << endl;
            exit(1);
        }
    }
} 

//
// Write "WhenStarted" from pls.h5 and write to cmp.h5
//
void WriteMetricWhenStarted(
        HDFCmpFile<CmpAlignment>         & cmpReader,
        HDFPlsReader                     & hdfPlsReader,
        const string                     & movieName) {
    string metric = "WhenStarted";
    string whenStarted;
    if (hdfPlsReader.scanDataReader.useWhenStarted == false) {
        cout << "ERROR! Attempting to read WhenStarted from " 
             << movieName 
             << " but the attriubte does not exist." << endl;
        exit(1);
    }
    hdfPlsReader.scanDataReader.ReadWhenStarted(whenStarted);

    if (!cmpReader.movieInfoGroup.whenStartedArray.IsInitialized()) {
        cmpReader.movieInfoGroup.whenStartedArray.Initialize(cmpReader.movieInfoGroup.movieInfoGroup, metric);
    }
    cmpReader.movieInfoGroup.whenStartedArray.Write(&whenStarted, 1);
}

//
// Print metrics.
//
string MetricsToString(const vector<string> & metrics) {
    string ret = ""; 
    int j = 0; 
    for (int i = 0; i < metrics.size(); i++) {
        ret += metrics[i]; 
        if (i != metrics.size()-1) ret += ","; 
        if (i % 4 == 3) ret += "\n";
    }
    return ret; 
}

// 
// Print usage.
//
void PrintUsage() {
    cout << "  loadPulses - Load pulse information and quality values into a Compare file" << endl;
    cout << "usage: loadPulses movieFile cmpFile [-metrics m1,m2,...] [-byread]" << endl;
    cout << "  movieFile may be a movie file or a fofn of movie file names." << endl;
    cout << "  metrics m1,m2,... is a comma-separated list (without spaces) of metrics " << endl
         << "  to print to the pulse file." << endl;
    cout << "  Valid metrics are: " << endl;
    cout << MetricsToString(GetAllSupportedMetrics(false)) << endl;
    //     << "    QualityValue, ClassifierQV, MergeQV," << endl
    //     << "    StartFrame, PulseWidth, pkmid, IPD, Light" << endl
    //     << "    WhenStarted, StartTimeOffset, PreBaseFrames," << endl
    //     << "    InsertionQV, DeletionQV, DeletionTag, SubstitutionQV" << endl
    //     << "    SubstitutionTag, PulseIndex, WidthInFrames" << endl;
    cout << "  By default, " << MetricsToString(GetDefaultMetrics()) << " are added" << endl;
    // Deprecate -useccs, an option for old data. 
    // cout << "  -useccs  This option is for older cmp.h5 files that do not have the read type " << endl
    //      << "    stored.  Newer cmp.h5 files have a read type that indicates the cmp.h5 file " << endl
    //      << "    has alignments generated from de novo ccs sequences.  Using this flag assuems"<<endl
    //      << "    ALL alignments in the cmp.h5 file are from ccs sequences, and loads the "<< endl
    //      << "    quality values from ccs instead of the raw sequence. "<<endl 
    //      << "  The only metrics that are allowed for de novo ccs sequences are QualityValue, " << endl
    //      << "  InsertionQV, DeletionQV, and SubstitutionQV" << endl;
    cout << "  -byread  Reads pulse/base fields by read, rather than reading an entire " << endl
         << "    movie first.  This uses considerably less memory than the defualt mode" << endl
         << "    but is slow." << endl;
    cout << "  -byMetric  Loads every pls/base field for each movie entirely before loading " << endl
         << "    another field. This uses more memory than -byread, but can be faster." << endl
         << "    This opiton is experimental. " << endl;
    cout << "  Using hdf version " << H5_VERS_MAJOR << "." << H5_VERS_MINOR << "." << H5_VERS_RELEASE << endl;
}


//
// The main function.
//
int main(int argc, char* argv[]) {
    string program = "loadPulses";
    string versionStr(VERSION);
    AppendPerforceChangelist(PERFORCE_VERSION_STRING, versionStr);

    string cmpFileName, movieFileName;
    int argi = 3;
    int numMetrics = 8;
    map<string,bool> metricOptions;
    int maxElements = 0;
    //Maximum Memory allowed for bymetric is 6 GB
    int maxMemory = 4; 
    //
    // Default is all options are false
    //
    CreateMetricOptions(metricOptions);
    string metricList = "";
    bool useCcsOnly = false;
    bool byRead = false;
    bool byMetric = false;
    bool failOnMissingData = false;

    CommandLineParser clp;
    clp.SetProgramName(program);
    clp.SetVersion(versionStr);

    clp.RegisterStringOption("basFileName", &movieFileName, 
            "The input {bas,pls}.h5 or input.fofn.", true);
    clp.RegisterStringOption("cmpFileName", &cmpFileName, 
            "The cmp.h5 file to load pulse information into.", true);
    clp.RegisterPreviousFlagsAsHidden();

    string metricsDescription = "A comma separated list of metrics (with no spaces).\nValid options are:\n";
    metricsDescription += MetricsToString(GetAllSupportedMetrics(false));
    metricsDescription += "\nDefault options are:\n";
    metricsDescription += MetricsToString(GetDefaultMetrics());

    clp.RegisterStringOption("metrics", &metricList, metricsDescription);
    clp.RegisterFlagOption("failOnMissingData", &failOnMissingData, 
            "Exit if any data fields are missing from the bas.h5 or pls.h5 "
            "input that are required to load a metric. Defualt is a warning.");
    clp.RegisterFlagOption("byread", &byRead, 
            "Load pulse information by read rather than buffering metrics.");
    clp.RegisterFlagOption("bymetric", & byMetric, 
            "Load pulse information by metric rather than by read. "
            "This uses more memory than -byread, but can be faster.");
    clp.RegisterIntOption("maxElements", &maxElements, 
            "Set a limit on the size of pls/bas file to buffer in with -bymetric "
            "(default value: maximum int). Use -byread if the limit is exceeded.", 
            CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("maxMemory", & maxMemory, 
            "Set a limit (in GB) on the memory to buffer data with -bymetric "
            "(default value: 4 GB). Use -byread if the limit is exceeded.",
             CommandLineParser::PositiveInteger);
    string progSummary = ("Loads pulse information such as inter pulse "
            "distance, or quality information into the cmp.h5 file. This allows " 
            "one to analyze kinetic and quality information by alignment column.");
    clp.SetProgramSummary(progSummary);
    clp.ParseCommandLine(argc, argv);

    cerr << "[INFO] " << GetTimestamp() << " [" << program << "] started." << endl;
    //use byMetric by default unless byRead is specified.
    byMetric = true;
    if (byRead) {
        byMetric = false;
    } 

    if (metricList == "") {
        SetDefaultMetricOptions(metricOptions);
    }
    else {
        ParseMetricsList(metricList, metricOptions);
    }

    // 
    // Always read in basecalls since they are used to check the sanity
    // of the alignment indices.
    //
    metricOptions["Basecall"] = true;

    //
    // Translate from the metrics to be loaded to the ones that are
    // required to compute them.
    // Need to be refactored.
    //
    vector<string> datasetFields;
    RequirementMap fieldRequirements;
    BuildRequirementMap(fieldRequirements);
    StoreDatasetFieldsFromPulseFields(metricOptions, fieldRequirements, datasetFields);
    

    //e.g. /PATH_TO_FILE/m120321_032600_42142_c100310572550000001523013208061210_s1_p0.bas.h5
    //     /PATH_TO_FILE/m120321_032600_42142_c100310572550000001523013208061210_s2_p0.bas.h5
    vector<string> movieFileNames;

    //e.g. m120321_032600_42142_c100310572550000001523013208061210_s1_p0
    //     m120321_032600_42142_c100310572550000001523013208061210_s2_p0
    vector<string> fofnMovieNames;

    FileOfFileNames::StoreFileOrFileList(movieFileName, movieFileNames);

    HDFBasReader hdfBasReader;
    HDFPlsReader hdfPlsReader;
    HDFCCSReader<SMRTSequence> hdfCcsReader;

    vector<string> baseFileFields, pulseFileFields;
    int fieldIndex;
    bool useBaseFile = false, usePulseFile = false;
    for (fieldIndex = 0; fieldIndex < datasetFields.size(); fieldIndex++) {
        if (hdfBasReader.ContainsField(datasetFields[fieldIndex])) {
            useBaseFile = true;
            baseFileFields.push_back(datasetFields[fieldIndex]);
        }
    }

    if (maxElements != 0) {
        hdfBasReader.maxAllocNElements = maxElements;
        hdfPlsReader.maxAllocNElements = maxElements;
    }

    //
    // For now, all runs will attempt to use information from a .bas
    // file, since it's assumed that if one has alignments, one has a
    // .bas file.
    //
    useBaseFile = true;
    //
    // Add some default fields.
    //
    hdfBasReader.IncludeField("Basecall");
    hdfBasReader.IncludeField("PulseIndex");
    hdfBasReader.InitializeFields(baseFileFields);

    for (fieldIndex = 0; fieldIndex < datasetFields.size(); fieldIndex++) {
        if (hdfPlsReader.ContainsField(datasetFields[fieldIndex])) {
            usePulseFile = true;
            pulseFileFields.push_back(datasetFields[fieldIndex]);
        }
    }
    if (usePulseFile) {
        // set hdfPlsReader.includedFields[fieldX] to true if fieldX is
        // in pulseFileFields 
        hdfPlsReader.InitializeFields(pulseFileFields);
    }
    hdfPlsReader.IncludeField("NumEvent");

    int nMovies = movieFileNames.size();
    int movieIndex;
    MovieNameToArrayIndex movieNameMap;
    //
    // Initialize movies. This accomplishes two tasks.  First, all movie
    // files are opened and initialized, so that if there are data
    // fields missing the program will exit now rather than in the
    // middle of loading pulses.  
    // Next, a list of movie names is created in fofnMovieNames.  The
    // cmp file does not necessarily index movies in the order of the
    // fofn, and so when loading pulses from a movie indexed by a cmp
    // file, one needs to look up the file name of the movie.  This is
    // done by scanning the fofnMovieNames list in order until the movie
    // is found. 

    //
    // h5 file access property list can be customized here.
    // 
    H5::FileAccPropList fileAccPropList = H5::FileAccPropList::DEFAULT;
    int    mdc_nelmts  = 4096; // h5: number of items in meta data cache
    size_t rdcc_nelmts = 4096; // h5: number of items in raw data chunk cache
    size_t rdcc_nbytes = 9192; // h5: raw data chunk cache size (in bytes) per dataset
    double rdcc_w0     = 0.75;    // h5: preemption policy
    // fileAccPropList.getCache(mdc_nelmts, rdcc_nelmts, rdcc_nbytes, rdcc_w0);
    fileAccPropList.setCache(mdc_nelmts, rdcc_nelmts, rdcc_nbytes, rdcc_w0);
    // fileAccPropList.setCache(4096, 4096, 8388608, rdcc_w0);


    // If one of the h5 in the fofn is a ccs.h5 file, then only load pulse
    // information from group /PulseData/ConsensusBaseCalls.
    for (movieIndex = 0; movieIndex < nMovies; movieIndex++) {
        FileType fileType;
        BaseSequenceIO::DetermineFileTypeByExtension(movieFileNames[movieIndex], fileType, true);
        if (fileType == HDFCCSONLY) {
            useCcsOnly = true;
        }
    }

    for (movieIndex = 0; movieIndex < nMovies; movieIndex++) {
        if (useCcsOnly) {
            hdfCcsReader.SetReadBasesFromCCS();
            hdfBasReader.SetReadBasesFromCCS();
        }
        if (!hdfBasReader.Initialize(movieFileNames[movieIndex], fileAccPropList)) {
            cout << "ERROR, could not initialize HDF file "
                << movieFileNames[movieIndex] << " for reading bases." << endl;
            exit(1);
        }
        else {
            fofnMovieNames.push_back(hdfBasReader.GetMovieName());
            movieNameMap[hdfBasReader.GetMovieName()] = movieIndex;
            hdfBasReader.Close();
        }

        // 
        // The pulse file is optional.  
        //
        if (usePulseFile) {
            if (hdfPlsReader.Initialize(movieFileNames[movieIndex], fileAccPropList) == 0) {
                usePulseFile = false;
            }
        }
   }

    CmpFile cmpFile;

    //
    // These readers pull information from the same pls file.
    //
    HDFCmpFile<CmpAlignment> cmpReader;

    if (cmpReader.Initialize(cmpFileName, H5F_ACC_RDWR) == 0) {
        cout << "ERROR, could not open the cmp file." << endl;
        exit(1);
    }

    if (cmpReader.HasNoAlignments()) {
        cout << "WARNING, there is no alignment in the cmp file." << endl;
        if (useBaseFile) {
            hdfBasReader.Close();
        }
        if (usePulseFile) {
            hdfPlsReader.Close();
        }
        cmpReader.Close();
        cerr << "[INFO] " << GetTimestamp() << " [" << program << "] ended." << endl;
        exit(0);
    }

    cmpReader.Read(cmpFile, false);

    // Sanity check: if there is a ccs.h5 file in the fofn and 
    // cmp.h5 file's readType is not CCS, something is wrong. 
    if (cmpFile.readType != ReadType::CCS and useCcsOnly) {
        cout << "ERROR, there is a ccs.h5 file in the fofn, while read type of" 
             << " the cmp.h5 file is not CCS." << endl;
        exit(1);
    }

    string commandLine;
    clp.CommandLineToString(argc, argv, commandLine);
    cmpReader.fileLogGroup.AddEntry(commandLine, "Loading pulse metrics", program, GetTimestamp(), versionStr);

    //
    // Group alignment indices by movie so that they may be processed one movie at a time
    // later on.  The movie indices set keeps track of all indices
    // listed in alignment files.  This keeps a reference to all
    // alignments in memory at once.   At the time of writing this, most
    // projects will have at most a few million alignments, and so the
    // size of this structure is modest.
    // Each movieIndexSets[$movieId] contains indices of all the alignments, which
    // are associated with a movie whose id in dataset /MovieInfo/ID equals $movieId
    //
    UInt alignmentIndex;
    map<int, vector<int> > movieIndexSets;

    for (alignmentIndex = 0; alignmentIndex < cmpFile.alnInfo.alignments.size(); alignmentIndex++) {
        movieIndexSets[cmpFile.alnInfo.alignments[alignmentIndex].GetMovieId()].push_back(alignmentIndex);
    }

    //
    // Load pulses from movies in order they appear in the input fofn.
    //
    int m;
    int fofnMovieIndex;
    for (fofnMovieIndex = 0; fofnMovieIndex < fofnMovieNames.size(); fofnMovieIndex++) {
        bool byMetricForThisMovie = byMetric;

        if (cmpFile.readType == ReadType::CCS or useCcsOnly) {
            hdfBasReader.SetReadBasesFromCCS();
            hdfCcsReader.Initialize(movieFileNames[fofnMovieIndex], fileAccPropList);
        }
        hdfBasReader.Initialize(movieFileNames[fofnMovieIndex], fileAccPropList);
        
        BaseFile  baseFile;
        PulseFile pulseFile;

        //
        // Deprecate reading the entire bas.h5 file. Reads are scanned
        // one by one or by metric, instead of caching all. 
        // It is still necessary to read in some of the datasets entirely,
        // in particular the start positions and hole numbers.
        //
        hdfBasReader.ReadBaseFileInit(baseFile);
				hdfBasReader.PrepareForRandomAccess();
        set<uint32_t> moviePartHoleNumbers;
        copy(baseFile.holeNumbers.begin(), baseFile.holeNumbers.end(), inserter(moviePartHoleNumbers, moviePartHoleNumbers.begin()));

        if (usePulseFile) {
            hdfPlsReader.Initialize(movieFileNames[fofnMovieIndex], fileAccPropList);
            hdfPlsReader.IncludeField("NumEvent");
            hdfPlsReader.IncludeField("StartFrame");
            //
            // Deprecate reading the entire pls.h5 file.
            // Reads are scanned by read or by metric instead of caching all.  
            // It is still necessary to read in some of the datasets entirely,
            // in particular the start positions and hole numbers.
            //
            hdfPlsReader.ReadPulseFileInit(pulseFile); 

        }

        string cmpFileMovieName;

        for (m = 0; m < cmpFile.movieInfo.name.size(); m++) {
            //
            // First find the file name for the movie 'm'
            //
            cmpFileMovieName = cmpFile.movieInfo.name[m];

            if (baseFile.GetMovieName() == cmpFileMovieName) {
                break;
            }
        }

        //
        // If the movie specified in the input.fofn is not found in the
        // cmp file, that indicates something bad is happeing.  Either the
        // input.fofn was not used to generate the cmp.h5 file, or no
        // alignments were found between the input bas.h5 and the
        // reference.  That shouldn't happen.
        // 
        if (m == cmpFile.movieInfo.name.size()) {
            cout << "WARNING: Could not find any alignments for file " << movieFileNames[fofnMovieIndex] << endl;
            continue;
        }

        //
        // Open the movie and load its pulses into memory.
        //
        movieIndex = cmpFile.movieInfo.id[m];
        UInt movieAlignmentIndex;
        
        //
        // Since usePulseFile is set when the input file is a pulseFile,
        // and ReadType::CCS becomes the read type when the alignments are
        // ccs, when pulse files are specified for de novo ccs alignments,
        // they will be opened as pulse files.  Since the de novo ccs
        // sequences do not have pulse file information, the auto-reading
        // of pulse files needs to be disabled.  Do that here.
        //
        if (cmpFile.readType == ReadType::CCS or useCcsOnly) {
            usePulseFile = false;
        }

        // Check whether all metrics are computable or not. 
        CanMetricsBeComputed(metricOptions, hdfBasReader, hdfPlsReader, useBaseFile,
            usePulseFile, failOnMissingData, movieFileNames[fofnMovieIndex]);

        // Get all metrics that are (1) supported, (2) required and (3) can be loaded.
        vector<string> metricsToLoad = GetMetricsToLoad(metricOptions);

        //
        // An index set is a set of indices into the alignment array that
        // are of reads generated by this movie.  Load pulses for all
        // alignments generated for this movie.
        //
        // Movie index sets should be sorted by alignment index. Build a lookup table for this.
        //
        
        std::vector<std::pair<int,int> > toFrom;
        UInt totalAlnLength = 0;
        for (movieAlignmentIndex = 0; movieAlignmentIndex < movieIndexSets[movieIndex].size(); movieAlignmentIndex++) {
            alignmentIndex = movieIndexSets[movieIndex][movieAlignmentIndex];
            totalAlnLength += cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetEnd() - \
                              cmpFile.alnInfo.alignments[alignmentIndex].GetOffsetBegin();
            toFrom.push_back(std::pair<int,int>(cmpFile.alnInfo.alignments[alignmentIndex].GetAlignmentId(), movieAlignmentIndex));
        }

        // orders by first by default.
        std::sort(toFrom.begin(), toFrom.end());

        // 
        // Check metric dataset size in this movie and the required memory 
        // consumption, if either limit is exceeded, switch to byread.
        //
        if (byMetricForThisMovie)
        {
            UInt requiredMem = ComputeRequiredMemory(metricsToLoad, hdfBasReader, 
                    hdfPlsReader, useBaseFile, usePulseFile, cmpReader, totalAlnLength);
            if (hdfBasReader.baseArray.arrayLength > hdfBasReader.maxAllocNElements or 
                (usePulseFile and 
                 hdfPlsReader.GetStartFrameSize() > hdfPlsReader.maxAllocNElements) or
                ((float)requiredMem / 1024 / 1024) > maxMemory) {
                cout << "Either the number of elements exceeds maxElement (" 
                     << hdfPlsReader.maxAllocNElements << "). Or the estimated memory " << endl
                     << "consumption exceeds maxMemory ("  << maxMemory << " GB)." << endl
                     << "Loading pulses from " << movieFileNames[fofnMovieIndex] 
                     << " by read." << endl;
                byMetricForThisMovie = false;
            }
        }

        if (((metricOptions.find("StartFrameBase") != metricOptions.end() and 
              metricOptions["StartFrameBase"]) or 
             (metricOptions.find("StartFramePulse")!= metricOptions.end() and 
              metricOptions["StartFramePulse"])) and !byMetricForThisMovie) {
            // Sneaky metrics StartFrameBase and StartFramePulse can used 
            // with -bymetric only
            cout << "ERROR: Internal metrics StartFrameBase and StartFramePulse "
                 << "can only be loaded with -bymetric." << endl;
            exit(1);
        }

        // Load "WhenStarted" before processing the others.
        if (metricOptions["WhenStarted"]) {
            WriteMetricWhenStarted(cmpReader, hdfPlsReader, movieFileNames[fofnMovieIndex]);
        }

        // Now load frame rate.
        if (AnyFieldRequiresFrameRate(datasetFields)) {
            if (useBaseFile) {
                cmpReader.movieInfoGroup.StoreFrameRate(m, baseFile.GetFrameRate());
            } else
            if (usePulseFile) {
                cmpReader.movieInfoGroup.StoreFrameRate(m, pulseFile.GetFrameRate());
            }
        }

        //
        // Load metrics for alignments from movie 'movieIndex'.
        //
        cout << "loading " <<  movieIndexSets[movieIndex].size() << " alignments for movie " << movieIndex << endl;

        UInt i;
        if (byMetricForThisMovie) {
            //
            // Build lookup tables for all alignments which 
            // are generated by the movie and check whether 
            // pls/bas.h5 and cmp.h5 match.
            //
            vector<MovieAlnIndexLookupTable> lookupTables;
           
            BuildLookupTablesAndMakeSane(cmpFile,     
                baseFile,     
                pulseFile,
                cmpReader,   
                hdfBasReader, 
                hdfPlsReader, 
                hdfCcsReader,
                useBaseFile, 
                usePulseFile,   
                useCcsOnly,
                movieIndexSets[movieIndex],
                toFrom,
                moviePartHoleNumbers,
                lookupTables);
            
            //
            // Group lookup tables by refGroupIndex and readGroupIndex.
            //
            vector<pair<UInt, UInt> >  groupedLookupTablesIndexPairs;
            GroupLookupTables(lookupTables, groupedLookupTablesIndexPairs);
            
            if (cmpFile.readType == ReadType::CCS or useCcsOnly) {
                vector<unsigned int> numPassesMetric;
                numPassesMetric.resize(lookupTables.size());
                UInt index = 0;  
                for (index = 0; index < lookupTables.size(); index++) {
                    if (lookupTables[index].skip) {
                        continue;
                    }
                    numPassesMetric[index] = hdfCcsReader.GetNumPasses(lookupTables[index].readIndex);
                }
                if (!cmpReader.alnInfoGroup.numPasses.IsInitialized()) {
                    cmpReader.alnInfoGroup.InitializeNumPasses();
                    // Clear /AlnInfo/NumPasses dataset.
                    cmpReader.alnInfoGroup.numPasses.Resize(0);
                }
                // Append numPasses of this movie to the end of /AlnInfo/NumPasses.
                UInt numPassesSize = cmpReader.alnInfoGroup.numPasses.size();
                cmpReader.alnInfoGroup.numPasses.WriteToPos(
                        &numPassesMetric[0], 
                        numPassesMetric.size(),
                        numPassesSize);
            }

            // Keep a list of currently cached fields. 
            vector<Field> cachedFields;
            if (usePulseFile) {
                // PulseCalls/ZMW/NumEvent is always cached in plsFile.
                cachedFields.push_back(Field("NumEvent", PlsField));
            } 

            for (int metricsToLoadIndex = 0; metricsToLoadIndex < metricsToLoad.size(); metricsToLoadIndex++) {
                string curMetric = metricsToLoad[metricsToLoadIndex];
                // Metric "WhenStarted" should have been loaded before getting here.
                if (curMetric == "WhenStarted") {
                    continue;
                }
                // Get the next metric to load.
                string nextMetric = "";
                if (metricsToLoadIndex+1 < metricsToLoad.size()) {
                    nextMetric = metricsToLoad[metricsToLoadIndex+1];
                }

                // Cache all required data for computing this metric.
                CacheRequiredFieldsForMetric(baseFile,
                    pulseFile,
                    hdfBasReader,
                    hdfPlsReader,
                    hdfCcsReader,
                    useBaseFile,
                    usePulseFile,
                    useCcsOnly,
                    cachedFields,
                    curMetric);

                // Compute the metric and write it to cmp.h5.
                WriteMetric(cmpFile, 
                    baseFile, 
                    pulseFile,
                    cmpReader, 
                    hdfBasReader,
                    hdfPlsReader,
                    hdfCcsReader,
                    useBaseFile,
                    usePulseFile, 
                    useCcsOnly,
                    lookupTables, 
                    groupedLookupTablesIndexPairs,
                    curMetric);

                // Clear cached fields unless they are required by the next metric.
                ClearCachedFields(baseFile,
                    pulseFile,
                    hdfBasReader,
                    hdfPlsReader,
                    hdfCcsReader,
                    useBaseFile,
                    usePulseFile,
                    useCcsOnly,
                    cachedFields,
                    curMetric,
                    nextMetric);
            }

            // Clear the default field "NumEvent"
            if (usePulseFile) {
                hdfPlsReader.ClearField(pulseFile, "NumEvent");
            }

        } else { // byRead for this movie
            for (movieAlignmentIndex = 0; movieAlignmentIndex < movieIndexSets[movieIndex].size(); movieAlignmentIndex++) {
                MovieAlnIndexLookupTable lookupTable; 
                BuildLookupTable(movieAlignmentIndex,
                    cmpFile,
                    baseFile,
                    usePulseFile,
                    pulseFile,
                    cmpReader,
                    movieIndexSets[movieIndex], 
                    toFrom,
                    moviePartHoleNumbers,
                    lookupTable);

                // Skip this alignment if it is not generated by this movie
                if (lookupTable.skip) {
                    continue;
                }

                UInt & alignmentIndex = lookupTable.alignmentIndex;
                int  & refGroupIndex  = lookupTable.refGroupIndex;
                int  & readGroupIndex = lookupTable.readGroupIndex;
                UInt & holeNumber     = lookupTable.holeNumber;
                int  & readIndex      = lookupTable.readIndex;
                int  & queryStart     = lookupTable.queryStart;
                int  & queryEnd       = lookupTable.queryEnd;
                int  & readStart      = lookupTable.readStart;
                int  & readLength     = lookupTable.readLength;
                UInt & offsetBegin    = lookupTable.offsetBegin;
                UInt & offsetEnd      = lookupTable.offsetEnd;

                string alignedSequence = GetAlignedSequenceFromCmpFile(cmpReader, lookupTable);

                // Create a map of where.
                vector<int> baseToAlignmentMap; 
                CreateSequenceToAlignmentMap(alignedSequence, baseToAlignmentMap);

			    // Condense gaps in the alignment for easy comparison.
                RemoveGaps(alignedSequence, alignedSequence);
                
                // Get source read.
                unsigned int numPasses;
                SMRTSequence sourceRead;
                GetSourceRead(cmpFile, 
                    baseFile    ,
                    pulseFile   ,
                    hdfBasReader, 
                    hdfPlsReader,
                    hdfCcsReader,
                    useBaseFile , 
                    usePulseFile,
                    useCcsOnly  , 
                    lookupTable , 
                    alignedSequence,
                    sourceRead  , 
                    numPasses);
                
                string readSequence;
                readSequence.resize(queryEnd - queryStart);
                copy((char*) (sourceRead.seq + queryStart),
                     (char*) (sourceRead.seq + queryEnd),
                     readSequence.begin());

                if (alignedSequence.size() != readSequence.size() or alignedSequence != readSequence) {
                    cout << "ERROR, the query sequence does not match the aligned query sequence." << endl;
                    cout << "HoleNumber: "<< holeNumber << ", MovieName: " << cmpFileMovieName;
                    cout << ", ReadIndex: " << (int) readIndex << 
                    cout << ", qStart: "<< queryStart << ", qEnd: " << queryEnd << endl;
                    cout << "Aligned sequence: "<< endl;
                    cout << alignedSequence << endl;
                    cout << "Original sequence: " << endl;
                    cout << readSequence << endl;
                    assert(0);
                }

                //
                // Compute any necessary data fields.  These usually involve
                // using differences of pulse indices, pulse widths, etc..
                // Missing fields are stored as 0's. 
                //

                vector<float> readPulseMetric;
                vector<float> floatMetric;
                vector<UChar> qvMetric;
                vector<HalfWord> frameRateMetric;
                vector<uint32_t> timeMetric;
                int ungappedAlignedSequenceLength = alignedSequence.size();
                assert(ungappedAlignedSequenceLength == queryEnd - queryStart);

                int alignedSequenceLength = offsetEnd - offsetBegin;
                readPulseMetric.resize(alignedSequenceLength+1);
                qvMetric.resize(alignedSequenceLength+1);
                frameRateMetric.resize(alignedSequenceLength+1);
                timeMetric.resize(alignedSequenceLength+1);

                UInt i;
                UInt pi;

                HDFCmpExperimentGroup* expGroup = cmpReader.refAlignGroups[refGroupIndex]->readGroups[readGroupIndex];
                UInt alnArrayLength = expGroup->alignmentArray.size();

                if (cmpFile.readType == ReadType::CCS or useCcsOnly) {
                    if (!cmpReader.alnInfoGroup.numPasses.IsInitialized()) {
                        cmpReader.alnInfoGroup.InitializeNumPasses();
                    }
                    cmpReader.alnInfoGroup.numPasses.WriteToPos(&numPasses, 1, alignmentIndex);
                }

                if (metricOptions["StartTimeOffset"] == true) {
                    if (!expGroup->startTimeOffset.IsInitialized()) {
                        expGroup->startTimeOffset.Initialize(expGroup->experimentGroup, "StartTimeOffset");
                    }
                    unsigned int readStartTimeOffset = sourceRead.startFrame[queryStart];
                    expGroup->startTimeOffset.WriteToPos(&readStartTimeOffset, 1, alignmentIndex);
                }

                if (metricOptions["QualityValue"] == true) {
                    if (!expGroup->qualityValue.IsInitialized()) {
                        expGroup->qualityValue.Initialize(expGroup->experimentGroup, "QualityValue", true, alnArrayLength);
                    }
                    // Store QualityValue. 
                    fill(qvMetric.begin(), qvMetric.end(), missingQualityValue);
                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        qvMetric[baseToAlignmentMap[i]] = sourceRead.qual[queryStart + i];
                    }
                    qvMetric[qvMetric.size()-1] = 0;
                    expGroup->qualityValue.WriteToPos(&qvMetric[0], qvMetric.size(), offsetBegin);
                }

                if (metricOptions["InsertionQV"] == true) {
                    if (!expGroup->insertionQV.IsInitialized()) {
                        expGroup->insertionQV.Initialize(expGroup->experimentGroup, "InsertionQV", true, alnArrayLength);
                    }
                    // Store InsertionQV.
                    fill(qvMetric.begin(), qvMetric.end(), missingQualityValue);
                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        qvMetric[baseToAlignmentMap[i]] = sourceRead.insertionQV[queryStart+ i];
                    }
                    qvMetric[qvMetric.size()-1] = 0;
                    expGroup->insertionQV.WriteToPos(&qvMetric[0], qvMetric.size(), offsetBegin);
                }

                if (metricOptions["MergeQV"] == true) {
                    if (!expGroup->mergeQV.IsInitialized()) {
                        expGroup->mergeQV.Initialize(expGroup->experimentGroup, "MergeQV", true, alnArrayLength);
                    }
                    // Store MergeQV. 
                    fill(qvMetric.begin(), qvMetric.end(), missingQualityValue);
                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        qvMetric[baseToAlignmentMap[i]] = sourceRead.mergeQV[queryStart+ i];
                    }
                    qvMetric[qvMetric.size()-1] = 0;
                    expGroup->mergeQV.WriteToPos(&qvMetric[0], qvMetric.size(), offsetBegin);
                }

                if (metricOptions["DeletionQV"] == true) {
                    if (!expGroup->deletionQV.IsInitialized()) {
                        expGroup->deletionQV.Initialize(expGroup->experimentGroup, "DeletionQV", true, alnArrayLength);
                    }
                    // Store DeletionQV. 
                    fill(qvMetric.begin(), qvMetric.end(), missingQualityValue);
                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        qvMetric[baseToAlignmentMap[i]] = sourceRead.deletionQV[queryStart+i];
                    }
                    qvMetric[qvMetric.size()-1] = 0;
                    expGroup->deletionQV.WriteToPos(&qvMetric[0], qvMetric.size(), offsetBegin);
                }

                if (metricOptions["DeletionTag"] == true) {
                    if (!expGroup->deletionTag.IsInitialized()) {
                        expGroup->deletionTag.Initialize(expGroup->experimentGroup, "DeletionTag", true, alnArrayLength);
                    }
                    vector<char> readDeletionTagMetric;
                    readDeletionTagMetric.resize(readPulseMetric.size());
                    // Store DeletionTag.
                    for (i = 0; i < readDeletionTagMetric.size()-1; i++ ) {
                        readDeletionTagMetric[i] = '-';
                    }
                    readDeletionTagMetric[i] = '\0';
                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        assert(baseToAlignmentMap[i] < readDeletionTagMetric.size());
                        readDeletionTagMetric[baseToAlignmentMap[i]] = sourceRead.deletionTag[queryStart+i];
                    }
                    readDeletionTagMetric[readDeletionTagMetric.size()-1] = 0;
                    expGroup->deletionTag.WriteToPos(&readDeletionTagMetric[0], readDeletionTagMetric.size(), offsetBegin);
                }

                if (metricOptions["PulseIndex"] == true) {
                    if (!expGroup->pulseIndex.IsInitialized()) {
                        expGroup->pulseIndex.Initialize(expGroup->experimentGroup, "PulseIndex", true, alnArrayLength);
                    }
                    vector<uint32_t> readPulseIndexMetric;
                    fill(readPulseIndexMetric.begin(), readPulseIndexMetric.end(), missingPulseIndex);
                    readPulseIndexMetric.resize(readPulseMetric.size());
                    // Store Pulse Index.
                    assert(readPulseIndexMetric.size() > 0);
                    for (i = 0; i < readPulseIndexMetric.size(); i++ ) {
                        readPulseIndexMetric[i] = 0;
                    }
                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        readPulseIndexMetric[baseToAlignmentMap[i]] = sourceRead.pulseIndex[queryStart+i];
                    }
                    readPulseIndexMetric[readPulseIndexMetric.size()-1] = 0;
                    expGroup->pulseIndex.WriteToPos(&readPulseIndexMetric[0], readPulseIndexMetric.size(), offsetBegin);
                } 

                if (metricOptions["SubstitutionTag"] == true) {
                    if (!expGroup->substitutionTag.IsInitialized()) {
                        expGroup->substitutionTag.Initialize(expGroup->experimentGroup, "SubstitutionTag", true, alnArrayLength);
                    }
                    vector<char> readSubstitutionTagMetric;
                    readSubstitutionTagMetric.resize(readPulseMetric.size());
                    // Store substitutionTag
                    for (i = 0; i < readSubstitutionTagMetric.size()-1; i++ ) {
                        readSubstitutionTagMetric[i] = '-';
                    }
                    readSubstitutionTagMetric[i] = '\0';
                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        readSubstitutionTagMetric[baseToAlignmentMap[i]] = sourceRead.substitutionTag[queryStart+i];
                    }
                    readSubstitutionTagMetric[readSubstitutionTagMetric.size()-1] = 0;
                    expGroup->substitutionTag.WriteToPos(&readSubstitutionTagMetric[0], readSubstitutionTagMetric.size(), offsetBegin);
                }

                if (metricOptions["SubstitutionQV"] == true) {
                    if (!expGroup->substitutionQV.IsInitialized()) {
                        expGroup->substitutionQV.Initialize(expGroup->experimentGroup, "SubstitutionQV", true, alnArrayLength);
                    }

                    // Store start time normalized to frame rate.
                    fill(qvMetric.begin(), qvMetric.end(), missingQualityValue);

                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        qvMetric[baseToAlignmentMap[i]] = sourceRead.substitutionQV[queryStart+i];
                    }
                    qvMetric[qvMetric.size()-1] = 0;
                    expGroup->substitutionQV.WriteToPos(&qvMetric[0], qvMetric.size(), offsetBegin);
                }

                if (metricOptions["ClassifierQV"] == true) {
                    if (!expGroup->classifierQV.IsInitialized()) {
                        expGroup->classifierQV.Initialize(expGroup->experimentGroup, "ClassifierQV", true, alnArrayLength);
                    }
                    fill(floatMetric.begin(), floatMetric.end(), NaN);

                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        floatMetric[baseToAlignmentMap[i]] = sourceRead.classifierQV[i+queryStart];
                    }
                    qvMetric[qvMetric.size()-1] = 0;
                    expGroup->classifierQV.WriteToPos(&floatMetric[0], floatMetric.size(), offsetBegin);
                }

                if (metricOptions["StartFrame"] == true) {
                    if (!expGroup->startTime.IsInitialized()) {
                        expGroup->startTime.Initialize(expGroup->experimentGroup, "StartFrame", true, alnArrayLength);	
                    }

                    // StartFrame used to be computed from baseFile.preBaseFrame and 
                    // baseFile.basWidthInFrames, whenever possible. But a more accurate
                    // way is to obtain StartFrame directly from pulseFile.StartFrame 
                    // when a pulseFile is provided.
                    if (usePulseFile) {
                        assert(sourceRead.startFrame);
                    } else if (useBaseFile) {
                        if (sourceRead.startFrame) {
                            Free(sourceRead.startFrame);
                        }
                        sourceRead.startFrame = new unsigned int[sourceRead.length];
                        copy(sourceRead.preBaseFrames, &sourceRead.preBaseFrames[sourceRead.length], sourceRead.startFrame);
                        for (i = 0; i < sourceRead.length-1; i++) {
                            sourceRead.startFrame[i+1] += sourceRead.widthInFrames[i];
                        }
                        partial_sum(sourceRead.startFrame, &sourceRead.startFrame[sourceRead.length],  sourceRead.startFrame);
                    }
                    
                    fill(timeMetric.begin(), timeMetric.end(), missingPulseIndex);
                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        timeMetric[baseToAlignmentMap[i]] = sourceRead.startFrame[i+queryStart];
                    }
                    timeMetric[timeMetric.size()-1] = 0;
                    expGroup->startTime.WriteToPos(&timeMetric[0], timeMetric.size(), offsetBegin);
                }

                if (metricOptions["PulseWidth"] == true) {
                    if (!expGroup->pulseWidth.IsInitialized()) {
                        expGroup->pulseWidth.Initialize(expGroup->experimentGroup, "PulseWidth", true, alnArrayLength);
                    }
                    fill(frameRateMetric.begin(), frameRateMetric.end(), missingFrameRateValue);
                    
                    // For legacy reasons, it's possible the width in frames is
                    // stored in the bas file. If this is the case, use the width
                    // in frames there.  Otherwise, use the width in frames stored
                    // in the pls file.
                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        frameRateMetric[baseToAlignmentMap[i]] = sourceRead.widthInFrames[queryStart + i];
                    }
                    frameRateMetric[frameRateMetric.size()-1] = 0;
                    expGroup->pulseWidth.WriteToPos(&frameRateMetric[0], frameRateMetric.size(), offsetBegin);
                }

                if (metricOptions["PreBaseFrames"] == true) {
                    if (!expGroup->preBaseFrames.IsInitialized()) {
                        expGroup->preBaseFrames.Initialize(expGroup->experimentGroup, "PreBaseFrames", true, alnArrayLength);
                    }
                    fill(frameRateMetric.begin(), frameRateMetric.end(), missingFrameRateValue);
                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        frameRateMetric[baseToAlignmentMap[i]] = sourceRead.preBaseFrames[i+queryStart];
                    }
                    frameRateMetric[frameRateMetric.size()-1] = 0;
                    expGroup->preBaseFrames.WriteToPos(&frameRateMetric[0], frameRateMetric.size(), offsetBegin);
                }

                if (metricOptions["WidthInFrames"] == true) {
                    if (!expGroup->widthInFrames.IsInitialized()) {
                        expGroup->widthInFrames.Initialize(expGroup->experimentGroup, "WidthInFrames", true, alnArrayLength);
                    }
                    // Compute width in frames. 
                    fill(frameRateMetric.begin(), frameRateMetric.end(), missingFrameRateValue);

                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        frameRateMetric[baseToAlignmentMap[i]] = sourceRead.widthInFrames[i+queryStart];
                    }
                    frameRateMetric[frameRateMetric.size()-1] = 0;
                    expGroup->widthInFrames.WriteToPos(&frameRateMetric[0], frameRateMetric.size(), offsetBegin);
                }

                if (metricOptions["pkmid"] == true) {
                    if (!expGroup->pkmid.IsInitialized()) {
                        expGroup->pkmid.Initialize(expGroup->experimentGroup, "pkmid", true, alnArrayLength);
                    }

                    for (i = 0; i < readPulseMetric.size(); i++ ) {
                        readPulseMetric[i] = NaN;
                    }

                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        readPulseMetric[baseToAlignmentMap[i]] = sourceRead.midSignal[i+queryStart];
                    }
                    readPulseMetric[readPulseMetric.size()-1] = 0;
                    expGroup->pkmid.WriteToPos(&readPulseMetric[0], readPulseMetric.size(), offsetBegin);
                }

                if (metricOptions["IPD"] == true) {
                    if (!expGroup->ipd.IsInitialized()) {
                        expGroup->ipd.Initialize(expGroup->experimentGroup, "IPD", true, alnArrayLength);
                    }
                    fill(frameRateMetric.begin(), frameRateMetric.end(), missingFrameRateValue);				

                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        //
                        // The IPD is undefined for the first base in a read.
                        //
                        if (usePulseFile ) {
                            if (queryStart == 0 and i == 0) {
                                frameRateMetric[baseToAlignmentMap[i]] = 0;
                            }
                            else {
                                frameRateMetric[baseToAlignmentMap[i]] = (sourceRead.startFrame[i+queryStart]  
                                        - sourceRead.startFrame[i+queryStart-1]
                                        - sourceRead.widthInFrames[i+queryStart-1]);
                            }
                        }
                        else if (useBaseFile) {
                            frameRateMetric[baseToAlignmentMap[i]] = sourceRead.preBaseFrames[i + queryStart];
                        }
                    }
                    frameRateMetric[frameRateMetric.size()-1] = 0;
                    expGroup->ipd.WriteToPos(&frameRateMetric[0], frameRateMetric.size(), offsetBegin);			
                }

                if (metricOptions["Light"] == true) {
                    if (!expGroup->light.IsInitialized()) {
                        expGroup->light.Initialize(expGroup->experimentGroup, "Light", true, alnArrayLength);
                    }
                    fill(frameRateMetric.begin(), frameRateMetric.end(), missingFrameRateValue);
                    for (i = 0; i < ungappedAlignedSequenceLength; i++ ) {
                        frameRateMetric[baseToAlignmentMap[i]] = sourceRead.meanSignal[i+queryStart];
                        frameRateMetric[baseToAlignmentMap[i]] = (frameRateMetric[baseToAlignmentMap[i]] * 
                                sourceRead.widthInFrames[i+queryStart]);
                    }
                    frameRateMetric[frameRateMetric.size()-1] = 0;
                    expGroup->light.WriteToPos(&frameRateMetric[0], frameRateMetric.size(), offsetBegin);			
                }

                sourceRead.Free();
                Free(sourceRead.meanSignal);
                Free(sourceRead.maxSignal);
                Free(sourceRead.midSignal);
                Free(sourceRead.startFrame);
                Free(sourceRead.classifierQV);
                Free(sourceRead.widthInFrames);
            }
        }

        if (useBaseFile) {
            hdfBasReader.Close();
        }
        if (cmpFile.readType == ReadType::CCS or useCcsOnly) {
            hdfCcsReader.Close();
        }
        if (usePulseFile) {
            hdfPlsReader.Close();
        }
    } // Done loading movies.
    
    cmpReader.Close();
    cerr << "[INFO] " << GetTimestamp() << " [" << program << "] ended." << endl;
}
