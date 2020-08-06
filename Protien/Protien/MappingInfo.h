#pragma once

#include <string>
#include <vector>
#include "FileSys.h"
#include "Logger.h"


using namespace std;
using namespace filesystem;

class MappingInfo {
public:
    MappingInfo(int filesCount);

    void AddTime(int fileIndex, int chunkIndex, int time);

    void AddToLog(const vector<path> fileNames);

private:
    vector<vector<vector<double>>> chunkInfo;
};

