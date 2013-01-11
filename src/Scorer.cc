#include "Scorer.h"

#include "MatchResult.h"

Score Scorer::scoreMatchedChunk(const MatchedChunk& chunk)
{
    if (chunk.isContigGap())
    {
        return scoreGap(chunk.getContigMatchLengthBp());
    }
    else
    {
        return scoreAlignment(chunk.getContigFragB(), chunk.getContigFragE(), chunk.getOpticalFragB(),
                              chunk.getOpticalFragE(), chunk.isBoundaryChunk());
    }
}

// Score the matched chunks in a MatchResult
void Scorer::scoreMatchResult(MatchResult * pResult)
{
    std::vector<MatchedChunk>& chunkList = pResult->matchedChunkList_;
    std::vector<Score>& scoreList = pResult->scoreList_;
    scoreList.clear();

    std::vector<MatchedChunk>::iterator iter = chunkList.begin();
    std::vector<MatchedChunk>::iterator E = chunkList.end();
    for(; iter != E; iter++)
    {
        scoreList.push_back(scoreMatchedChunk(*iter));
    }
}
