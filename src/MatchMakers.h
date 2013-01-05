#ifndef MATCHMAKERS_H
#define MATCHMAKERS_H

#include "MatchMaker.h"

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

    void filterResults(MatchResultPtrVec& matches);
    bool filterFunction(const MatchResult * pMatch);
    int maxMatches_;
};

class LocalMatchMaker : public MatchMaker
{

    public:
    LocalMatchMaker(size_t minChunks) :
        minChunks_(minChunks) 
    { };


    private:
    size_t minChunks_; // minimum number of chunks to accept in an alignment
};




#endif
