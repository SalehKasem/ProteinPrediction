#include "Serializer.h"


using namespace logging;

static const int _magicNumber = 1234423;


unsigned short crc16CalcChecksum(const char *data, unsigned long length, unsigned short crc) {
    while (length > 0) {
        crc = (crc >> 8) | (crc << 8);
        crc ^= *data;
        crc ^= (crc & 0xff) >> 4;
        crc ^= crc << 12;
        crc ^= (crc & 0xff) << 5;

        ++data;
        --length;
    }
    //LOG_DEBUG("crc =" + to_string(crc));

    return crc;
}

unsigned short calculateQuerySeedCrc(QuerySeed &qs, unsigned short crc) {
    crc = crc16CalcChecksum((char *) &qs.queryIndex, sizeof(int), crc);
    int s = qs.seedContent.size();
    crc = crc16CalcChecksum((char *) &s, sizeof(int), crc);
    crc = crc16CalcChecksum((char *) qs.seedContent.c_str(), s, crc);
    //LOG_DEBUG("QuerySeed CRC: " + to_string(crc));
    return crc;
}

unsigned short calculatePcnEntryCrc(PcnEntry &e, unsigned short crc) {
    int s = e.seedContent.size();
    crc = crc16CalcChecksum((char *) &s, sizeof(int), crc);
    crc = crc16CalcChecksum((char *) e.seedContent.c_str(), s, crc);
    s = e.sequenceName.size();
    crc = crc16CalcChecksum((char *) &s, sizeof(int), crc);
    crc = crc16CalcChecksum((char *) e.sequenceName.c_str(), s, crc);
    crc = crc16CalcChecksum((char *) &e.offsetInSeq, sizeof(int), crc);
    crc = crc16CalcChecksum((char *) &e.ancestorIndex, sizeof(int), crc);
    crc = crc16CalcChecksum((char *) &e.hasDuplicatedContent, sizeof(bool), crc);
    //LOG_DEBUG("PcnEntry CRC: " + to_string(crc));

    return crc;
}


unsigned short Serializer::write(ofstream &stream, const char *str, streamsize count, unsigned short crc) {
    stream.write(str, count);
    return crc16CalcChecksum(str, count, crc);
}

unsigned short Serializer::read(ifstream &stream, char *str, streamsize count, unsigned short crc) {
    stream.read(str, count);
    return crc16CalcChecksum(str, count, crc);
}


void Serializer::serializePcnMap(unordered_map<QuerySeed, unordered_map<int, vector<PcnEntry>>> &map,
                                 string const &outDirPath, int level) {

    string lvlDir = outDirPath + "/Level_" + to_string(level);
    Helpers::createDirectory(lvlDir.c_str());

    string path = lvlDir + "/" + "pcnDump_lvl_" + to_string(level) + ".pcnDump";
    LOG_INFO_BROADCAST(("Serialization to " + path + ". Level " + to_string(level)).c_str());
    ofstream dump(path, ios::binary);

    unsigned short crc;

    //serialize magic number
    crc = write(dump, (char *) &_magicNumber, sizeof(int), 0);

    //serialize current pcn level
    crc = write(dump, (char *) &level, sizeof(int), crc);

    int mapSize = (int) map.size();
    //serialize map size
    crc = write(dump, (char *) &mapSize, sizeof(int), crc);

    for (auto it = map.begin(); it != map.end(); ++it) {
        unordered_map<int, vector<PcnEntry>> &innerMap = it->second;
        //QuerySeed contains string inside. Use its serializer to serialize it
        QuerySeed qs = it->first;
        qs.serialize(dump);
        crc = calculateQuerySeedCrc(qs, crc);

        int innerMapSize = (int) innerMap.size();
        crc = write(dump, (char *) &innerMapSize, sizeof(int), crc);


        for (auto innerIt = innerMap.begin(); innerIt != innerMap.end(); ++innerIt) {
            crc = write(dump, (char *) &(innerIt->first), sizeof(int), crc);

            vector<PcnEntry> &vec = innerIt->second;
            int vecSize = (int) vec.size();
            crc = write(dump, (char *) &vecSize, sizeof(int), crc);

            for (int i = 0; i < vecSize; i++) {
                vec[i].serialize(dump);
                crc = calculatePcnEntryCrc(vec[i], crc);

            }
        }
    }
    //LOG_DEBUG("final crc = " + to_string(crc));
    dump.write((char *) &crc, sizeof(unsigned short));
    dump.flush();
    dump.close();
}

void Serializer::deserializePcnMap(string path, unordered_map<QuerySeed, unordered_map<int, vector<PcnEntry>>> *map,
                                   int *level) {
    LOG_INFO_BROADCAST(("Deserialization from " + path).c_str());

    if (!Helpers::isExists(path.c_str())) {
        throw runtime_error("cannot locate dump file");
    }
    unsigned short crc;
    int magicNumber;
    ifstream dump(path, ios::binary);

    crc = read(dump, (char *) &magicNumber, sizeof(int), 0);

    if (magicNumber != _magicNumber) {
        throw runtime_error("unexpected dump file content (magic number mismatch)");
    }

    int targetLevel;
    crc = read(dump, (char *) &targetLevel, sizeof(int), crc);
    unordered_map<QuerySeed, unordered_map<int, vector<PcnEntry>>> result;
    int mapSize;
    crc = read(dump, (char *) &mapSize, sizeof(int), crc);

    for (int i = 0; i < mapSize; i++) {
        QuerySeed seed;
        seed.deserialize(dump);
        crc = calculateQuerySeedCrc(seed, crc);
        int innerMapSize;
        crc = read(dump, (char *) &innerMapSize, sizeof(int), crc);

        unordered_map<int, vector<PcnEntry>> innerMap;
        for (int j = 0; j < innerMapSize; j++) {
            int level;
            int vecSize;

            crc = read(dump, (char *) &level, sizeof(int), crc);
            crc = read(dump, (char *) &vecSize, sizeof(int), crc);

            vector<PcnEntry> vec(vecSize);
            for (int k = 0; k < vecSize; k++) {
                PcnEntry entry;
                entry.deserialize(dump);
                crc = calculatePcnEntryCrc(entry, crc);
                vec[k] = entry;
            }
            innerMap.insert(make_pair(level, vec));
        }
        result.insert(make_pair(seed, innerMap));
    }
    //LOG_DEBUG("final crc = " + to_string(crc));

    unsigned short mapCrc;
    dump.read((char *) &mapCrc, sizeof(unsigned short));
    //LOG_DEBUG("deserialized crc = " + to_string(mapCrc));

    if (mapCrc != crc) {
        throw runtime_error("dump file is corrupted (CRC check failed)");
    }
    *map = result;
    *level = targetLevel;
}