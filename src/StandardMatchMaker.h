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
    StandardMatchMaker(Scorer * pScorer, int maxMatches, size_t minContigHits, double minLengthRatio, double maxMissRateContig, double avgChi2Threshold) :
        pScorer_(pScorer),
        maxMatches_(maxMatches),
        minContigHits_(minContigHits),
        minLengthRatio_(minLengthRatio),
        maxMissRateContig_(maxMissRateContig),
        avgChi2Threshold_(avgChi2Threshold)
    { };

    // Build MatchResults from pScoreMatrix. Place MatchResults in matches.
    bool makeMatches(const ScoreMatrix_t * pScoreMatrix, MatchResultPtrVec& matches,
                     const MapData * pOpticalMap, const MapData * pContigMap, bool contigIsForward);


    private:
    MatchResult * buildMatch(const Index_t& end_index, const ScoreMatrix_t * pScoreMatrix, const MapData * pOpticalMap,
                             const MapData * pContigMap, bool contigIsForward);
    bool filterFunction(const MatchResult * pMatch);
    void scoreMatch(MatchResult * pMatch); // Fill in the score attributes of the MatchResult

    Scorer * pScorer_;
    int maxMatches_;
    size_t minContigHits_; // minimum number of contig hits to accept for an alignment
    double minLengthRatio_; // minimum length ratio between aligned contig and optical fragments
    double maxMissRateContig_;
    double avgChi2Threshold_;
};

#endif
