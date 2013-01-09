#ifndef STANDARDMATCHMAKER_H
#define STANDARDMATCHMAKER_H

#include "MatchMaker.h"
#include "Scorer.h"

// The StandardMatchMaker builds global alignments
// The best non-overlapping matches for the ScoreMatrix are selected.
// A filter is applied to determine if the match is of acceptable quality.
class StandardMatchMaker : public MatchMaker
{
    public: 
    StandardMatchMaker(int maxMatches, Scorer * pScorer) :
        maxMatches_(maxMatches),
        pScorer_(pScorer)
        {};

    // Build MatchResults from pScoreMatrix. Place MatchResults in matches.
    bool makeMatches(const ScoreMatrix_t * pScoreMatrix, MatchResultPtrVec& matches,
                     const MapData * pOpticalMap, const MapData * pContigMap, bool contigIsForward);


    private:
    MatchResult * buildMatch(const Index_t& end_index, const ScoreMatrix_t * pScoreMatrix, const MapData * pOpticalMap,
                             const MapData * pContigMap, bool contigIsForward);
    bool filterFunction(const MatchResult * pMatch);
    void scoreMatch(MatchResult * pMatch); // Fill in the score attributes of the MatchResult
    int maxMatches_;
    Scorer * pScorer_;
};

#endif
