#pragma once

#include<string>
#include <fstream>
#include "Helpers.h"
#include "Structures.h"
#include "Logger.h"

using namespace std;

class Serializer {
private:
    // Disallow creating an instance of this object
    Serializer() {}

    //static unsigned short crc16CalcChecksum(const char* data, unsigned long length, unsigned short crc);
    static unsigned short write(ofstream &stream, const char *str, streamsize count, unsigned short crc);

    static unsigned short read(ifstream &stream, char *str, streamsize count, unsigned short crc);

public:

    static void
    serializePcnMap(unordered_map<QuerySeed, unordered_map<int, vector<PcnEntry>>> &map, string const &outDirPath,
                    int level);

    static void
    deserializePcnMap(string path, unordered_map<QuerySeed, unordered_map<int, vector<PcnEntry>>> *map, int *level);
};
