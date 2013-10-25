#ifndef STANDARDMATCHMAKER_H
#define STANDARDMATCHMAKER_H

#include <algorithm>
#include <iostream>
#include <vector>
using std::vector;

//#include "MatchMaker.h"
#include "MatchResult.h"
#include "MatchedChunk.h"
#include "debugUtils.h"
#include "Scorer.h"

#define BUILDMATCH_DEBUG 0
#define FILTER_DEBUG 0
#define MAKEMATCH_DEBUG 0

// The StandardMatchMaker builds global alignments
// The best non-overlapping matches for the ScoreMatrix are selected.
// A filter is applied to determine if the match is of acceptable quality.
class StandardMatchMaker
{
    public: 
    StandardMatchMaker(int maxMatches, size_t minContigHits, double minLengthRatio, double maxMissRateContig, double avgChi2Threshold) :
        maxMatches_(maxMatches),
        minContigHits_(minContigHits),
        minLengthRatio_(minLengthRatio),
        maxMissRateContig_(maxMissRateContig),
        avgChi2Threshold_(avgChi2Threshold)
    { };

    // Build MatchResults from pScoreMatrix. Place MatchResults in matches.
    template <typename Scorer>
    bool makeMatches(const ScoreMatrix_t * pScoreMatrix, Scorer * pScorer, MatchResultPtrVec& matches,
                     const MapData * pOpticalMap, const MapData * pContigMap, bool contigIsForward);

    // Build MatchResults from a seededScoreMatrix.
    //bool makeMatches(const seeded::ScoreMatrix& scoreMatrix, MatchResultVec& matches,
     //                bool contigIsForward);


    private:
    MatchResult * buildMatch(const Index_t& end_index, const ScoreMatrix_t * pScoreMatrix, const MapData * pOpticalMap,
                             const ContigMapData * pContigMap, bool contigIsForward);
    bool filterFunction(const MatchResult * pMatch);
    void scoreMatch(MatchResult * pMatch); // Fill in the score attributes of the MatchResult

    int maxMatches_;
    size_t minContigHits_; // minimum number of contig hits to accept for an alignment
    double minLengthRatio_; // minimum length ratio between aligned contig and optical fragments
    double maxMissRateContig_;
    double avgChi2Threshold_;

};


// Return true of the match result has an overlap with another match result in the iterator range from B to E.
bool hasOverlap(const MatchResult* pMatch, MatchResultPtrVec::const_iterator B, MatchResultPtrVec::const_iterator E);

template <typename Scorer>
bool StandardMatchMaker::makeMatches(const ScoreMatrix_t * pScoreMatrix, Scorer * pScorer, MatchResultPtrVec& matches,
                                     const MapData * pOpticalMap, const MapData * pContigMap, bool contigIsForward)
{

    matches.clear();

    const ContigMapData * pContigMapData = dynamic_cast<const ContigMapData*>(pContigMap);
    assert(pContigMapData->isForward() == contigIsForward);

    const int m = pScoreMatrix->m_; // num rows
    const int n = pScoreMatrix->n_; // num cols
    const int lr = m-1;
    const int lro = lr*n;
    const ScoreElement_t * pE;
    bool foundMatch = false;

    typedef std::pair<double, Index_t> ScoreIndexPair;
    typedef std::vector<ScoreIndexPair> ScoreIndexPairVec;
    ScoreIndexPairVec trailSeeds;

    // Check the last row in the dynamic programming table
    // for potential matches
    for(int j=1; j < n; j++)
    {
        pE = &pScoreMatrix->d_[lro + j];
        if (pE->score_ > -Constants::INF)
        {
            trailSeeds.push_back(ScoreIndexPair( pE->score_, Index_t(lr,j) ) );
        }
    }

    if (trailSeeds.empty()) return false;

    // Sort the scores in descending order. (The higher the score, the better).
    sort(trailSeeds.begin(), trailSeeds.end());
    reverse(trailSeeds.begin(), trailSeeds.end());

    const ScoreIndexPairVec::const_iterator E = trailSeeds.end();
    for( ScoreIndexPairVec::const_iterator iter = trailSeeds.begin();
         iter != E;
         iter++)
    {
        if ((maxMatches_ >= 0) && matches.size() == (size_t) maxMatches_) break;
        Index_t end_index = iter->second;

        #if MAKEMATCH_DEBUG > 0
        std::cout << "Making match for index: " << end_index <<
                     " score: " << iter->first << std::endl;
        #endif

        MatchResult * pMatch = buildMatch(end_index, pScoreMatrix, pOpticalMap, pContigMapData, contigIsForward);
        if (pMatch == NULL) continue;

        // Check that the match is acceptable. A match is OK if:
        // 1. It passes the filter funtion.
        // 2. It does not overlap a higher scoring match result.
        bool matchOK = filterFunction(pMatch);
        matchOK = matchOK && !hasOverlap(pMatch, matches.begin(), matches.end());
        if (!matchOK)
            delete pMatch;
        else
        {

            #if MAKEMATCH_DEBUG > 0
            std::cout << "Accepted Match!" << std::endl; 
            #endif

            pScorer->scoreMatchResult(pMatch);
            matches.push_back(pMatch);
            foundMatch = true;
        }
    }

    return foundMatch;
}


#endif
