#ifndef LOCALMATCHMAKER_H
#define LOCALMATCHMAKER_H

#include "MatchMaker.h"
#include <set>

class Scorer;

// The LocalMatchMaker builds local alignments.
class LocalMatchMaker : public MatchMaker
{

    public:
    LocalMatchMaker(Scorer * pScorer, int maxMatches, size_t minContigHits, double minLengthRatio, double maxMissRateContig, double avgChi2Threshold) :
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

    bool filterFunction(const MatchResult * pMatch);
    MatchResult * buildMatch(const Index_t& end_index, const ScoreMatrix_t * pScoreMatrix, const MapData * pOpticalMap,
                             const MapData * pContigMap, bool contigIsForward, set<Index_t>& usedCells);

    Scorer * pScorer_;
    int maxMatches_;
    size_t minContigHits_; // minimum number of contig hits to accept for an alignment
    double minLengthRatio_; // minimum length ratio between aligned contig and optical fragments
    double maxMissRateContig_;
    double avgChi2Threshold_;
};




#endif
