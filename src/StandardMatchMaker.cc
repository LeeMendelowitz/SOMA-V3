#include "MatchMakers.h"
#include "MatchResult.h"

#include <algorithm>
#include <iostream>

#define MATCHBUILDER_DEBUG 0

using namespace std;


// Return true of the match result has an overlap with another match result in the iterator range from B to E.
bool hasOverlap(const MatchResult* pMatch, MatchResultPtrVec::const_iterator B, MatchResultPtrVec::const_iterator E)
{
    bool hasOverlap = false;
    for(MatchResultPtrVec::const_iterator iter = B;
        iter != E;
        iter++)
    {
        if (pMatch->overlaps(*iter))
        {
            hasOverlap = true;
            break;
        }
    }
    return hasOverlap;
}

bool StandardMatchMaker::makeMatches(const ScoreMatrix_t * pScoreMatrix, MatchResultPtrVec& matches,
                                     const MapData * pOpticalMap, const MapData * pContigMap, bool contigIsForward)
{

    matches.clear();
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
        MatchResult * pMatch = buildMatch(end_index, pScoreMatrix, pOpticalMap, pContigMap, contigIsForward);
        if (pMatch == NULL) continue;

        // Check that the match is acceptable. A match is OK if:
        // 1. It passes the filter funtion.
        // 2. It does not overlap a high scoring match result.
        bool matchOK = filterFunction(pMatch);
        matchOK = matchOK && !hasOverlap(pMatch, matches.begin(), matches.end());
        if (!matchOK)
            delete pMatch;
        else
        {
            matches.push_back(pMatch);
            foundMatch = true;
        }
    }

    return foundMatch;
}


MatchResult * StandardMatchMaker::buildMatch(const Index_t& end_index, const ScoreMatrix_t * pScoreMatrix, const MapData * pOpticalMap,
                                             const MapData * pContigMap, bool contigIsForward)
{

    const vector<FragData>& contigFrags = pContigMap->getFrags(contigIsForward);
    const vector<FragData>& opticalFrags = pOpticalMap->getFrags();
    const int n = pScoreMatrix->n_;
    double score = pScoreMatrix->d_[n*end_index.first + end_index.second].score_;

    // Build the trail through scoreMatrix
    vector<Index_t> trail;
    trail.reserve(pScoreMatrix->m_);
    Index_t ind = end_index;
    trail.push_back(ind);
    while (true)
    {
        ScoreElement_t * pE = &(pScoreMatrix->d_[n*ind.first + ind.second]);
        ind = make_pair(pE->pi_, pE->pj_);
        if ((ind.first < 0) || (ind.second < 0) )
            break;
        trail.push_back(ind);
    }
    reverse(trail.begin(), trail.end());

    const vector<Index_t>::iterator tb = trail.begin();
    const vector<Index_t>::iterator te = trail.end();

    MatchResult * pMatch = new MatchResult(pContigMap->getId(), pOpticalMap->getId(),
                                           pContigMap->getLength(), contigIsForward, score);

    #if MATCHBUILDER_DEBUG > 0
    std::cout << "Building match for: "
              << pContigMap->getId() << " "
              << pOpticalMap->getId() << " "
              << "forward: " << contigIsForward
              << "end_index: " << end_index.first << " , " << end_index.second
              << std::endl;
    # endif

    // Build the vector of matched fragments. Add them directory to the match.
    vector<MatchedChunk>& matchedChunkList = pMatch->matchedChunkList_;
    matchedChunkList.reserve(trail.size());
    int pc, po; // Contig and optical end index (exclusive) of predecessor matched fragment
    int os, oe, cs, ce; // optical and contig start index inclusive and end index (exclusive) for matched fragment
    pc = tb->first;
    po = tb->second;
    for(vector<Index_t>::const_iterator ti = tb+1; ti != te; ti++)
    {
        // Add this matched fragment to the matchedChunkList
        // Reminder: matched fragment are given by the range of indexes
        // from opStart (inclusive) to opEnd (exclusive)

        // Here, indices refer to the Fragment index (not the site index)
        // ------|----------|-----------|-------- ...
        //   0         1           2         3    ...

        cs = pc; // contig start index (inclusive, 0 based)
        ce = ti->first; // contig end index (exclusive)
        os = po; // optical start index (inclusive, 0 based)
        oe = ti->second; // optical end index (exclusive, 0 based)

        int cStartBp = pContigMap->getStartBp(cs, contigIsForward);
        int cEndBp = pContigMap->getEndBp(ce-1, contigIsForward); // substract 1 since end is exclusive
        bool gapAlignment = (oe == os);
        int opStartBp, opEndBp;

        #if MATCHBUILDER_DEBUG > 0
        std::cout << "Contig: "
                  << "(" << cs << "," << ce << ")"
                  << " = " << "(" << cStartBp << "," << cEndBp << ")"
                  << std::endl;
        #endif
        if (gapAlignment)
        {
            // If this is a gap alignment, use the end of the last aligned
            // chunk as start & end optical map location
            opStartBp = pOpticalMap->getEndBp(os-1); // 
            opEndBp = opStartBp;
        }
        else
        {
            opStartBp = pOpticalMap->getStartBp(os);
            opEndBp = pOpticalMap->getEndBp(oe-1);
        }

        #if MATCHBUILDER_DEBUG > 0
        std::cout << "Optical: "
                  << "(" << os << "," << oe << ")"
                  << " = " << "(" << opStartBp << "," << opEndBp << ")"
                  << std::endl;
        #endif

        // For boundary cases where the matched chunk includes the first or last contig fragment,
        // the chunk is not bounded by two matched restriction sites. Adjust the optical start/end positions
        // accordingly.
        if (cs == 0)
        {
            // Adjust the optical start position
            opStartBp = opEndBp - (cEndBp - cStartBp);
            #if MATCHBUILDER_DEBUG > 0
            std::cout << "cs==0 adjustment: "
                      << "(" << os << "," << oe << ")"
                      << " = " << "(" << opStartBp << "," << opEndBp << ")"
                      << std::endl;
            #endif
        }
        else if ((size_t) ce == contigFrags.size())
        {
            opEndBp = opStartBp + (cEndBp - cStartBp);
            #if MATCHBUILDER_DEBUG > 0
            std::cout << "ce==contigFrags.size() adjustment: "
                      << "(" << os << "," << oe << ")"
                      << " = " << "(" << opStartBp << "," << opEndBp << ")"
                      << std::endl;
            #endif
        }

        MatchedChunk chunk = MatchedChunk(os, oe, opStartBp, opEndBp, opticalFrags,
                                          cs, ce, cStartBp, cEndBp, contigFrags);
        matchedChunkList.push_back(chunk);
        pc = ce; po = oe;
    }

    //////////////////////////////////////////////////////
    //#ifdef DEBUG
    // Check that there are at most two boundary chunks.
    int boundaryCount = 0;
    const vector<MatchedChunk>::const_iterator E = matchedChunkList.end();
    for (vector<MatchedChunk>::const_iterator iter = matchedChunkList.begin();
         iter != E;
         iter++)
    {
        if (iter->boundaryChunk_)
        {
            boundaryCount++;
            assert( iter == matchedChunkList.begin() || iter == (matchedChunkList.end()-1)) ;
        }
    }
    assert(boundaryCount <= 2);
    //#endif
    ////////////////////////////////////////////////////////

    pMatch->buildAlignmentAttributes();

    return pMatch;
}

// Return true if MatchResult is acceptable
bool StandardMatchMaker::filterFunction(const MatchResult * pResult)
{

    //Chi2 Filter
    bool failChi2Filter = false;
    if (pResult->numAlignedInnerBlocks_ > 0)
        failChi2Filter = (pResult->chi2_/((double) pResult->numAlignedInnerBlocks_)) > opt::avgChi2Threshold;
       
    bool failLengthRatio = pResult->alignedLengthRatio_ < opt::minLengthRatio;
    bool failContigMissRate = pResult->contigMissRate_ > opt::maxMissRateContig;
    bool failContigHitsCheck = pResult->contigHits_ <= 0;

    bool fail = (failChi2Filter || 
                 failLengthRatio ||
                 failContigMissRate ||
                 failContigHitsCheck);

    #if DEBUG_FILTER > 0
    std::cout << "MatchResult: " << *pResult << "\nt"
              << "FailChi2: " << failChi2Filter
              << " FailLengthRatio: " << failLengthRatio
              << " FailContigMissRate: " << failContigMissRate
              << " FailContigHitsCheck: " << failContigHitsCheck
              << std::endl;
    #endif

    return (!fail);
}
