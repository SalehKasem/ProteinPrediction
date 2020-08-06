#include "Helpers.h"


using namespace std::experimental::filesystem;
using namespace std;
using namespace logging;


#if defined(__GNUC__) || defined(__GNUG__)

#include <sys/stat.h>
#include <cstring>

int Helpers::create_dir(const char* filePath) {
	return mkdir(filePath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}

#elif defined(_MSC_VER)
#include <direct.h>
int Helpers::create_dir(const char* filepath) {
	return _mkdir(filepath);
}
#endif

int Helpers::getHammingDistance(const char* source, const char* target) {
	int srcLength = strlen(source);
	int trgLength = strlen(target);
	int distance = 0;
	int length = min(srcLength, trgLength);
	int difference = abs((int)(srcLength - trgLength));

	for (int i = 0; i < length; i++) {
		if (source[i] != target[i]) {
			distance++;
		}
	}

	return distance + difference;
}

int Helpers::getHammingDistance(const char* src, int srcOffset, int srcLength, const char* target, int targetOffset,
	int trgLength) {
	int distance = 0;
	int length = min(srcLength, trgLength);
	int difference = abs((int)(srcLength - trgLength));

	for (int i = 0; i < length; i++) {
		if (src[srcOffset + i] != target[targetOffset + i]) {
			distance++;
		}
	}

	return distance + difference;
}

string Helpers::substringOfCString(const char* cstr, size_t start, size_t length) {
	assert(start + length <= strlen(cstr));
	return std::string(cstr + start, length);
}

bool Helpers::strCmp(const char* str1, const char* str2) {
	int length = min(strlen(str1), strlen(str2));
	return memcmp(str1, str2, sizeof(char) * length) == 0;
}

string Helpers::GetSeedName(SimilarityCandidate* seed) {
	return seed->entry.seqPtr->seqData.substr(seed->entry.offset, seed->length);
}

list<string> Helpers::divideToSeeds(const string& pSequence, int sizeOfSeed) {
	if (pSequence.length() < sizeOfSeed) {
		return list<string>{};
	}
	if (pSequence.length() == sizeOfSeed) {
		return list<string>{pSequence};
	}

	list<string> seeds;

	for (int i = 0; i <= pSequence.length() - sizeOfSeed; i++) {
		seeds.push_back(pSequence.substr(i, sizeOfSeed));
	}

	return seeds;
}

void Helpers::createOutputResultDirectory(string const& path) {
	if (!experimental::filesystem::exists(path.c_str())) {
		Helpers::create_dir(path.c_str());
	}
}

void Helpers::createPcnFile(map<int, vector<PcnEntry>>* pcnPtr, int level, string const& outDirPath,
	int pcnLevelsCount) {
	vector<PcnEntry>& v = (*pcnPtr)[0];
	PcnEntry& p = v[0];

	string lvlDir = outDirPath + "/Level_" + to_string(level);

	if (!experimental::filesystem::exists(lvlDir.c_str())) {
		Helpers::create_dir(lvlDir.c_str());
	}

	string path = lvlDir + "/" + p.sequenceName + "_Seed-" + p.seedContent + ".pcn";
	ofstream pcnFile(path);

	pcnFile << "Index, SeedContent, LevelIndex, AncestorIndex, OffsetInSequence, SequenceName" << endl;

	int* levelsSize = new int[pcnPtr->size()];
	for (const auto& pcnLevel : *pcnPtr) {
		if (pcnLevel.first == 0) {
			levelsSize[pcnLevel.first] = pcnLevel.second.size();
			continue;
		}

		if (pcnLevel.first > pcnLevelsCount) {
#if defined(__GNUC__) || defined(__GNUG__)
			throw runtime_error("unexpected level in PCN");
#else
			throw exception("unexpected level in PCN");
#endif
		}
		levelsSize[pcnLevel.first] = levelsSize[pcnLevel.first - 1] + pcnLevel.second.size();
	}

	for (const auto& pcnLevel : *pcnPtr) {
		if (pcnLevel.first == 0) {
			continue;
		}

		for (int i = 0; i < pcnLevel.second.size(); i++) {
			int ancestorInd = (pcnLevel.first == 1) ? 0 :
				pcnLevel.second[i].ancestorIndex + levelsSize[pcnLevel.first - 2];

			pcnFile << i + levelsSize[pcnLevel.first - 1] << ","    //index in file
				<< pcnLevel.second[i].seedContent << ","        //Seed content
				<< pcnLevel.first << ","                        //PCN level
				<< ancestorInd << ","                            //ancestor index
				<< pcnLevel.second[i].offsetInSeq << ","        //offset in sequence
				<< pcnLevel.second[i].sequenceName << endl;        //sequence name
		}
	}
	pcnFile.close();
}

void Helpers::createPcnFileFromSimpleFormat(map<int, vector<PcnEntry>>* pcnPtr, int level, string const& outDirPath,
	int pcnLevelsCount) {
	vector<PcnEntry>& v = (*pcnPtr)[0];
	PcnEntry& p = v[0];

	string lvlDir = outDirPath + "/Level_" + to_string(level);

	if (!experimental::filesystem::exists(lvlDir.c_str())) {
		Helpers::create_dir(lvlDir.c_str());
	}

	string path = lvlDir + "/" + p.sequenceName + "_Seed-" + p.seedContent +"_Offset-"+ ".pcn";
	ofstream pcnFile(path);

	pcnFile << "Index, SeedContent, LevelIndex, AncestorIndex, OffsetInSequence, SequenceName" << endl;

	int* levelsSize = new int[pcnPtr->size()];
	for (const auto& pcnLevel : *pcnPtr) {
		if (pcnLevel.first == 0) {
			levelsSize[pcnLevel.first] = pcnLevel.second.size();
			continue;
		}

		if (pcnLevel.first > pcnLevelsCount) {
#if defined(__GNUC__) || defined(__GNUG__)
			throw runtime_error("unexpected level in PCN");
#else
			throw exception("unexpected level in PCN");
#endif
		}
		levelsSize[pcnLevel.first] = levelsSize[pcnLevel.first - 1] + pcnLevel.second.size();
	}

	for (const auto& pcnLevel : *pcnPtr) {
		if (pcnLevel.first == 0) {
			continue;
		}

		for (int i = 0; i < pcnLevel.second.size(); i++) {
			int ancestorInd = (pcnLevel.first == 1) ? 0 :
				pcnLevel.second[i].ancestorIndex + levelsSize[pcnLevel.first - 2];

			pcnFile << i + levelsSize[pcnLevel.first - 1] << ","    //index in file
				<< pcnLevel.second[i].seedContent << ","        //Seed content
				<< pcnLevel.first << ","                        //PCN level
				<< ancestorInd << ","                            //ancestor index
				<< pcnLevel.second[i].offsetInSeq << ","        //offset in sequence
				<< pcnLevel.second[i].sequenceName << endl;        //sequence name
		}
	}
	pcnFile.close();
}

void Helpers::buildBlackList(map<int, vector<PcnEntry>>* pcnPtr, list<string>* blackList, int level, int pcnLevelsCount) {
	vector<PcnEntry>& v = (*pcnPtr)[0];
	PcnEntry& p = v[0];


	int* levelsSize = new int[pcnPtr->size()];
	for (const auto& pcnLevel : *pcnPtr) {
		if (pcnLevel.first == 0) {
			levelsSize[pcnLevel.first] = pcnLevel.second.size();
			continue;
		}

		if (pcnLevel.first > pcnLevelsCount) {
#if defined(__GNUC__) || defined(__GNUG__)
			throw runtime_error("unexpected level in PCN");
#else
			throw exception("unexpected level in PCN");
#endif
		}
		levelsSize[pcnLevel.first] = levelsSize[pcnLevel.first - 1] + pcnLevel.second.size();
	}

	for (const auto& pcnLevel : *pcnPtr) {
		if (pcnLevel.first == 0) {
			continue;
		}
		for (int i = 0; i < pcnLevel.second.size(); i++) {
			string	name = pcnLevel.second[i].sequenceName;
			bool inBlackList = (std::find(blackList->begin(), blackList->end(), name) != blackList->end());
			if (pcnLevel.first == 1 && !inBlackList)
				blackList->push_back(name);
		}
	}
}

void Helpers::createPcnSimpleFilter(map<int, vector<PcnEntry>>* pcnPtr, int level, int pcnLevelsCount, SimpleFilteredNodesList* nodesLst, list<string>* blackList,
	bool noCollision = true, bool removeDoublicate = true) {
	vector<PcnEntry>& v = (*pcnPtr)[0];
	PcnEntry& p = v[0];

	string name;
	pair<int, int> tmp_indexs;


	int* levelsSize = new int[pcnPtr->size()];
	for (const auto& pcnLevel : *pcnPtr) {
		if (pcnLevel.first == 0) {
			levelsSize[pcnLevel.first] = pcnLevel.second.size();
			continue;
		}

		if (pcnLevel.first > pcnLevelsCount) {
#if defined(__GNUC__) || defined(__GNUG__)
			throw runtime_error("unexpected level in PCN");
#else
			throw exception("unexpected level in PCN");
#endif
		}
		levelsSize[pcnLevel.first] = levelsSize[pcnLevel.first - 1] + pcnLevel.second.size();
	}

	for (const auto& pcnLevel : *pcnPtr) {
		if (pcnLevel.first == 0) {
			continue;
		}
		for (int i = 0; i < pcnLevel.second.size(); i++) {
			list<int> pcn_levels;
			list<pair <int, int>> pcn_indexs;
			name = pcnLevel.second[i].sequenceName;

			bool inBlackList = (std::find(blackList->begin(), blackList->end(), name) != blackList->end());

			if (!inBlackList) {
				tmp_indexs.second = p.offsetInSeq;
				tmp_indexs.first = pcnLevel.second[i].offsetInSeq;
				SimpleFilteredNodes _node;
				_node.level = pcnLevel.first;
				_node.offsetInDest = p.offsetInSeq;
				_node.OffsetInSrc = pcnLevel.second[i].offsetInSeq;
				_node.seqName = pcnLevel.second[i].sequenceName;
				if (pcnLevel.first == 1)
					_node.lvl1Connection = true;
				nodesLst->add(_node,noCollision,removeDoublicate);
			}
		}
	}
}


bool Helpers::checkIfExist(vector<SimilarityCandidate>* candidates, int offset, const string& protName) {
	for (std::vector<SimilarityCandidate>::iterator candidoz = candidates->begin();
		candidoz != candidates->end(); ++candidoz) {
		if (candidoz->entry.offset == offset && candidoz->entry.seqPtr->seqName == protName) {
			return true;
		}
	}
	return false;
}

void Helpers::createDirectory(const char* path) {
	if (!experimental::filesystem::exists(path)) {
		Helpers::create_dir(path);
	}
}

bool Helpers::isExists(const char* path) {
	return experimental::filesystem::exists(path);
}
