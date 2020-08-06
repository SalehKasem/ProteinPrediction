#include "PcnBuilder.h"


#if defined(__GNUC__) || defined(__GNUG__)

#include <cstring>

#endif
using namespace std;
using namespace logging;

typedef steady_clock timer;

PcnBuilder::PcnBuilder(list<McsData>* mcsListPtr, unordered_map<string, vector<MappedEntry>>* allTextMapPtr,
	ConfigurationData const& confData) {
	mcsList = mcsListPtr;
	allTextMap = allTextMapPtr;
	sizeOfSeed = confData.SizeOfSeed;
	numOfAllowedMismatches = confData.NumOfAllowedMismatches;
	alphabet = confData.Alphabet;
	complexityThreshold = confData.ComplexityThreshold;

}

void ::PcnBuilder::appendPcnEntries(const QuerySeed* querySeedPtr, unordered_map<int, vector<PcnEntry>>* pcnPtr,
	vector<SimilarityCandidate>* subPcn, int curLevel) {
	bool isExists;
	bool duplicatedSeedContent;
	_cur_Level = curLevel;

	unordered_map<int, vector<PcnEntry>>::const_iterator it = pcnPtr->find(curLevel);
	if (it == pcnPtr->end()) {
		//we here, i.e. the curLevel is still not created for pcnPtr;
		(*pcnPtr)[curLevel] = vector<PcnEntry>();
	}

	vector<PcnEntry>& curLvlPcn = (*pcnPtr)[curLevel];

	for (size_t i = 0; i < subPcn->size(); i++) {
		isExists = false;
		duplicatedSeedContent = false;

		//perform search on the last level of Pcn
		//vector<PcnEntry> lastLevel = pcnPtr[lastLevel];
		for (size_t j = 0; j < curLvlPcn.size(); j++) {
			if (Helpers::strCmp(&((*subPcn)[i].entry.seqPtr->seqData)[(*subPcn)[i].entry.offset],
				curLvlPcn[j].seedContent.c_str())) {
				//we here i.e. - similarityCandidata and pcnEntry from curLevel have same seedContent
				if (curLvlPcn[j].offsetInSeq == (*subPcn)[i].entry.offset &&
					Helpers::strCmp(curLvlPcn[j].sequenceName.c_str(), (*subPcn)[i].entry.seqPtr->seqName.c_str())) {
					//we here, i.e. pcnEntryv with same address already exists in curret PCN. No need to add it twice.
					isExists = true;
					break;
				}

				duplicatedSeedContent = true;
			}
		}

		if (isExists) {
			continue;
		}

		//go over levels of pcn
		for (const auto& pcnLevel : *pcnPtr) {
			if (pcnLevel.first == curLevel) {
				//we already performed search in curlevel in separate "for" above
				continue;
			}

			for (size_t j = 0; j < pcnLevel.second.size(); j++) {
				if (pcnLevel.second[j].offsetInSeq == (*subPcn)[i].entry.offset &&
					Helpers::strCmp(pcnLevel.second[j].sequenceName.c_str(),
					(*subPcn)[i].entry.seqPtr->seqName.c_str())) {
					//we here, i.e. pcnEntryvwith same address already exists in curret PCN. No need to add it twice.
					isExists = true;
					break;
				}
			}

			if (isExists) {
				break;
			}
		}

		if (!isExists) {
			PcnEntry p;
			p.offsetInSeq = (*subPcn)[i].entry.offset;
			p.seedContent = Helpers::substringOfCString((*subPcn)[i].entry.seqPtr->seqData.c_str(), p.offsetInSeq,
				sizeOfSeed);
			p.sequenceName = (*subPcn)[i].entry.seqPtr->seqName;
			p.ancestorIndex = (*subPcn)[i].ancestorIndex;
			p.hasDuplicatedContent = duplicatedSeedContent;
			//p.ancestorPtr = cand->ancestor;
			curLvlPcn.push_back(p);
		}
	}
}

bool PcnBuilder::checkComplexity(const MappedEntry& candidate, int length) {

	for (int i = 0; i < alphabet.length(); i++) {
		int count = 0;
		for (int j = 0; j < length; j++) {
			if (candidate.seqPtr->seqData[(int)candidate.indexInSeq + j] == alphabet[i]) {
				count++;
			}
		}

		if (count / length > complexityThreshold) {
			LOG_INFO_BROADCAST(("Candidate: " + candidate.seqPtr->seqData.substr(candidate.indexInSeq, length)
				+ " rejected with complexity threshold " + std::to_string(count / length) +
				" by letter " + alphabet[i] +
				+". Sequence: " + candidate.seqPtr->seqName + ". Index: " +
				std::to_string(candidate.indexInSeq) + ".").c_str());
			return false;
		}
	}
	return true;
}

string PcnBuilder::createPcnName(string pProteinName, string pSeedName) {
	ostringstream result;
	result << pProteinName << "_" << pSeedName;
	return result.str();
}

string PcnBuilder::createProteinDirName(string pProteinName) {
	return pProteinName;
}

bool PcnBuilder::canCreatePcn(string pSeed, string pProteinName, PcnData* pcnData) {
	if (pSeed.empty() || pProteinName.empty()) {
		return false;
	}

	string pcnName = createPcnName(pProteinName, pSeed);
	string dirName = createProteinDirName(pProteinName);
	vector<PcnEntry> vec;
	pcnData->seedFileName = pcnName;
	pcnData->proteinDirName = dirName;
	return true;
}


void
PcnBuilder::findSimilarToSeed(const char* seedSeqName, const char* seedPtr,
	int seedIndex, int bucket, vector<SimilarityCandidate>* candidates, int current_level)
{
	int frameSize = 0;

	for (list<McsData>::iterator currentMcs = (*mcsList).begin(); currentMcs != (*mcsList).end(); ++currentMcs) {
		bool seedIsAlive = true;
		int currentFrameIndex = 0;
		int length = 0;
		frameSize = (*currentMcs).formLength;

		//while not reached end of current seed
		while (seedIsAlive) {
			std::ofstream outfile;
			outfile.open(PcnBuilder::OutputResultDirectory + "/Results.txt", std::ios_base::app);
			string seedCurrentFrame = Helpers::substringOfCString(seedPtr, currentFrameIndex, frameSize);

			//dvigaemsya mcs'om po seed'u
			if ((*currentMcs).missIndexes.size() > 0) {
				//MCS form with dots(misses)
				for (std::list<int>::iterator it = (*currentMcs).missIndexes.begin();
					it != (*currentMcs).missIndexes.end(); ++it) {
					seedCurrentFrame[*it] = '.';
				}
			}

			seedIsAlive = ((int)currentFrameIndex + 1 + frameSize) <= strlen(seedPtr);

			if (allTextMap[bucket].find(seedCurrentFrame) == allTextMap[bucket].end()) {
				currentFrameIndex++;
				continue;
			}

			for (std::vector<MappedEntry>::iterator candidate = allTextMap[bucket][seedCurrentFrame].begin();
				candidate != allTextMap[bucket][seedCurrentFrame].end(); ++candidate) {
				//length - is number of characters since matched frame offset to end of candidate sequence
				//Sometimes, the length of candidate sequence may be shorter than size of seed.
				length = min((int)(candidate->seqPtr->seqData.length() - candidate->indexInSeq), sizeOfSeed);

				//The candidate should be compared against seed in such way that position of current MCS element frame is the same in both seed AND candidate
				//i.e. before calculation of candidate we always need check that its offset in sequence is greater or equal to currentFrameIndex of seed.
				//Else we need calculate hamming distance between substring of seed and candidate that starts from zero offset.
				//In addition we need add to calculated distance the difference between the seed and its substring taken for hamming distance calculation

				int candidateStartOffset = candidate->indexInSeq - currentFrameIndex;

				if (candidateStartOffset < 0 || length < sizeOfSeed) {
					continue;
				}

				if (Helpers::checkIfExist(candidates, candidateStartOffset, candidate->seqPtr->seqName)) {
					continue;
				}

				int hammDistance = Helpers::getHammingDistance(seedPtr, 0, strlen(seedPtr),
					candidate->seqPtr->seqData.c_str(), candidateStartOffset,
					length);

				if (hammDistance <= numOfAllowedMismatches && checkComplexity(*candidate, length)) {
					SimilarityCandidate ocCand;
					ocCand.entry.offset = candidateStartOffset;
					ocCand.entry.seqPtr = candidate->seqPtr;
					ocCand.length = length;
					ocCand.ancestorIndex = seedIndex;
					candidates->push_back(ocCand);


					outfile << seedSeqName << "," << seedPtr <<
						"," << candidate->seqPtr->seqName <<
						"," << candidate->seqPtr->seqData.substr(candidateStartOffset, 20) <<
						"," << candidateStartOffset << "," << current_level << endl;
				}
			}
			currentFrameIndex++;
		}
	}
}