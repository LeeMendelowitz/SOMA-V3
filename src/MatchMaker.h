// Class to build and select MatchResults
// from a ScoreMatrix.

// This is a pure abstract class
#ifndef MATCHMAKER_H
#define MATCHMAKER_H

#include <vector>

#include "ScoreMatrix.h"
#include "MatchResult.h"

typedef std::vector<MatchResult *> MatchResultPtrVec;

class MatchMaker
{
    public:

    virtual ~MatchMaker() {};

    // Given a pScoreMatrix, populate a vector of MatchResults.
    virtual bool makeMatches(const ScoreMatrix_t * pScorematrix, MatchResultPtrVec& matches,
                     const MapData * referenceMap, const MapData * queryMap, bool queryIsForward) = 0;
};

#endif
