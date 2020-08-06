#ifndef _CONFIGURATIONDATA_H_

#define _CONFIGURATIONDATA_H_

#include <string>
#include <iostream>
#include <fstream>
#include <ostream>
#include "json/json.h"
using namespace std;


class ConfigurationData {


public:
    void load(const string &);

    string toString();


public:
    string InputDataDirectory;
    string McsFilePath;
    string TargetSearchSequencesFilePath;
    string OutputResultDirectory;
    string LogDirectory;
    string Alphabet;
    int NumOfThreads;
    int SizeOfSeed;
    int NumOfAllowedMismatches;
    int PcnLevels;
    int FormBufferSize;
    int SequencesPerChunk;
    int MaxMatchesInPcn;
    double ComplexityThreshold;
};

#endif
