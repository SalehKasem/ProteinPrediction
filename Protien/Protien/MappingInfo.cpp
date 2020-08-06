#include "MappingInfo.h"

using namespace logging;
using namespace filesystem;

MappingInfo::MappingInfo(int filesCount)
{
	for(int i = 0; i< filesCount; i++)
	{
		chunkInfo.emplace_back(vector<vector<double>>());
	}
}

void MappingInfo::AddTime(int fileIndex, int chunkIndex, int time)
{

	if(chunkInfo[fileIndex].size() <= chunkIndex)
	{
		chunkInfo[fileIndex].push_back(vector<double>());
	}
	chunkInfo[fileIndex][chunkIndex].push_back(time);
}

void MappingInfo::AddToLog(vector<path> filenames)
{
	LOG_INFO_BROADCAST("Mapping summary:");
	for (size_t fileIndex = 0; fileIndex < filenames.size(); fileIndex++)
	{
		LOG_INFO_BROADCAST(("File: " + filenames[fileIndex].string()).c_str());

		for (size_t i = 0; i < chunkInfo[fileIndex].size(); i++)
		{
			ostringstream ss;
			ss << "\tChunk#" << i << ":";
			for (size_t j = 0; j < chunkInfo[fileIndex][i].size(); j++)
			{
				ss << chunkInfo[fileIndex][i][j];
				if (j < chunkInfo[fileIndex][i].size() - 1)
				{
					ss << ", ";
				}
				else
				{
					ss << endl;
				}
			}
			LOG_INFO_BROADCAST(ss);
		}
	}
}


