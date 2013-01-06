#ifndef MATCHMAKERS_H
#define MATCHMAKERS_H

#include "MatchMaker.h"

#include <set>

// The StandardMatchMaker builds global alignments
// The best non-overlapping matches for the ScoreMatrix are selected.
// A filter is applied to determine if the match is of acceptable quality.
class StandardMatchMaker : public MatchMaker
{
    public: 
    StandardMatchMaker(int maxMatches = -1) : maxMatches_(maxMatches) {};

    // Build MatchResults from pScoreMatrix. Place MatchResults in matches.
    bool makeMatches(const ScoreMatrix_t * pScoreMatrix, MatchResultPtrVec& matches,
                     const MapData * pOpticalMap, const MapData * pContigMap, bool contigIsForward);


    private:
    MatchResult * buildMatch(const Index_t& end_index, const ScoreMatrix_t * pScoreMatrix, const MapData * pOpticalMap,
                             const MapData * pContigMap, bool contigIsForward);
    bool filterFunction(const MatchResult * pMatch);
    int maxMatches_;
};


// The LocalMatchMaker builds local alignments.
class LocalMatchMaker : public MatchMaker
{

    public:
    LocalMatchMaker(int maxMatches, size_t minContigHits, double minLengthRatio, double maxMissRateContig, double avgChi2Threshold) :
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

    int maxMatches_;
    size_t minContigHits_; // minimum number of contig hits to accept for an alignment
    double minLengthRatio_; // minimum length ratio between aligned contig and optical fragments
    double maxMissRateContig_;
    double avgChi2Threshold_;
};




#endif
