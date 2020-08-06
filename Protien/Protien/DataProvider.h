#pragma once

#include "Structures.h"
#include <chrono>
#include <ostream>
#include <map>
#include "Logger.h"
#include "Logger.h"
#include "Structures.h"
#include "FileSys.h"
using namespace filesystem;
using namespace std;


class DataProvider {
	int maxMcsFormLength = 0;

public:
	DataProvider(int sizeOfSeed);

	list<McsData> readMcsFile(const string& mcsFileName);

	void getFiles(const string& directoryPath, vector<path>* files);

	void readSequanceNames(map<string, string>* seqDataList, const string& dataFileName, size_t* fileOffset, int sequencesPerChunk);

	void readDataFile(vector<FastaEntry>* seqDataList, const string& dataFileName, size_t* fileOffset,
		int sequencesPerChunk);

	void loadQueries(vector<QueryData>* pQueries, const string& queriesFileName);

	list<McsData> mcsList;
};