#include "Scorer.h"

#include "MatchResult.h"
/*
Score Scorer::scoreMatchedChunk(MatchedChunk& chunk)
{
    Score score;

    if (chunk.isContigGap())
    {
        score = scoreGap(chunk.getContigMatchLengthBp());
    }
    else
    {
        score = scoreAlignment(chunk.getContigFragB(), chunk.getContigFragE(), chunk.getOpticalFragB(),
                               chunk.getOpticalFragE(), chunk.isBoundaryChunk());
    }

    chunk.setScore(score);
    return score;
}

// Score the matched chunks in a MatchResult
void Scorer::scoreMatchResult(MatchResult * pResult)
{
    std::vector<MatchedChunk>& chunkList = pResult->matchedChunkList_;
    std::vector<MatchedChunk>::iterator iter = chunkList.begin();
    std::vector<MatchedChunk>::iterator E = chunkList.end();
    for(; iter != E; iter++)
        scoreMatchedChunk(*iter);
}
*/
