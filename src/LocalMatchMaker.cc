#include "LocalMatchMaker.h"
#include "MatchResult.h"

#include <algorithm>
#include <iostream>
#include <set>

#include "debugUtils.h"

#define MATCHBUILDER_DEBUG 0
#define FILTER_DEBUG 0

using namespace std;

// The ScoreMatrix for local alignment is filled with entries with scores of zero or greater
bool LocalMatchMaker::makeMatches(const ScoreMatrix_t * pScoreMatrix, MatchResultPtrVec& matches,
                                     const MapData * pOpticalMap, const MapData * pContigMap, bool contigIsForward)
{

    matches.clear();
    const int m = pScoreMatrix->m_; // num rows
    const int n = pScoreMatrix->n_; // num cols
    const ScoreElement_t * pE;
    bool foundMatch = false;

    typedef std::pair<double, Index_t> ScoreIndexPair;
    typedef std::vector<ScoreIndexPair> ScoreIndexPairVec;
    ScoreIndexPairVec trailSeeds;

    // Check the entire ScoreMatrix for potential matches
    // These are cells with a positive score
    for (int i=1; i < m; i++)
    {
        const int ro = n*i;
        for(int j=1; j < n; j++)
        {
            pE = &pScoreMatrix->d_[ro + j];
            if (pE->score_ > 0.0)
                trailSeeds.push_back(ScoreIndexPair( pE->score_, Index_t(i,j) ) );
        }
    }

    if (trailSeeds.empty()) return false;

    // Sort the scores in descending order. (The higher the score, the better).
    sort(trailSeeds.begin(), trailSeeds.end());
    reverse(trailSeeds.begin(), trailSeeds.end());

    // Keep track of which cells of the ScoreMatrix have already been used to build matches.
    set<Index_t> usedCells;

    const ScoreIndexPairVec::const_iterator E = trailSeeds.end();
    for( ScoreIndexPairVec::const_iterator iter = trailSeeds.begin();
         iter != E;
         iter++)
    {
        if ((maxMatches_ >= 0) && matches.size() == (size_t) maxMatches_) break;
        Index_t end_index = iter->second;
        MatchResult * pMatch = buildMatch(end_index, pScoreMatrix, pOpticalMap, pContigMap, contigIsForward, usedCells);
        if (pMatch == NULL) continue;
        pScorer_->scoreMatchResult(pMatch);
        matches.push_back(pMatch);
        foundMatch = true;
    }

    return foundMatch;
}


// Build a Local Alignment.
// Check that this alignment does not use any cells in the ScoreMatrix that
// have already been used in a higher scoring alignment.
// If the Match is not acceptable, return NULL
MatchResult * LocalMatchMaker::buildMatch(const Index_t& end_index, const ScoreMatrix_t * pScoreMatrix, const MapData * pOpticalMap,
                                          const MapData * pContigMap, bool contigIsForward,
                                          set<Index_t>& usedCells)
{

    const vector<FragData>& contigFrags = pContigMap->getFrags();
    const int n = pScoreMatrix->n_;
    double score = pScoreMatrix->d_[n*end_index.first + end_index.second].score_;

    // Build the trail through scoreMatrix
    #if MATCHBUILDER_DEBUG > 0
    {
    ScoreElement_t * pE = &(pScoreMatrix->d_[n*end_index.first + end_index.second]);
    std::cout << "Building trail from element: "
              << pE
              << std::endl;
    }
    #endif
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

        #if MATCHBUILDER_DEBUG > 0
        {
        ScoreElement_t * pE = &(pScoreMatrix->d_[n*ind.first + ind.second]);
        std::cout << "Adding element to trail: "
                  << pE
                  << std::endl;
        }
        #endif


        trail.push_back(ind);
    }
    reverse(trail.begin(), trail.end());

    const vector<Index_t>::iterator tb = trail.begin();
    const vector<Index_t>::iterator te = trail.end();

    //////////////////////////////////////
    // If any of the cells have already been used, do not build a Match.
    bool allCellsUnused = true;
    for (vector<Index_t>::iterator iter = tb;
         iter != te;
         iter++)
    {
        if (usedCells.count(*iter) > 0)
        {
            allCellsUnused = false;
            break;
        }
    }

    if (!allCellsUnused)
    {
        // The MatchResult for this trail should not be built since it uses cells that were
        // already used in a previous MatchResult.
        return NULL;
    }
    /////////////////////////////////////

    MatchResult * pMatch = new MatchResult(pContigMap->getId(), pOpticalMap->getId(), 
                                           pContigMap->getLength(), contigIsForward, score);

    #if MATCHBUILDER_DEBUG > 0
    std::cout << "Building match for: "
              << pContigMap->getId() << " "
              << pOpticalMap->getId() << " "
              << "forward: " << contigIsForward
              << " end_index: " << end_index.first << " , " << end_index.second
              << " trailSize: " << trail.size()
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

        MatchedChunk chunk = MatchedChunk(os, oe, opStartBp, opEndBp, pOpticalMap,
                                          cs, ce, cStartBp, cEndBp, pContigMap);
        matchedChunkList.push_back(chunk);
        pc = ce; po = oe;
    }

    pMatch->buildAlignmentAttributes();

    bool matchOK = filterFunction(pMatch);
    if (!matchOK)
    {
        delete pMatch;
        pMatch = NULL;
        return NULL;
    }

    // Mark the cells in the trail for the MatchResult as used.
    usedCells.insert(tb, te);
    return pMatch;
}

// Return true if MatchResult is acceptable
bool LocalMatchMaker::filterFunction(const MatchResult * pResult)
{

    //Chi2 Filter
    bool failChi2Filter = false;
    if (pResult->numAlignedInnerBlocks_ > 0)
        failChi2Filter = (pResult->chi2_/((double) pResult->numAlignedInnerBlocks_)) > avgChi2Threshold_;
       
    bool failLengthRatio = pResult->alignedLengthRatio_ < minLengthRatio_;
    bool failContigMissRate = pResult->contigMissRate_ > maxMissRateContig_;
    bool failContigHitsCheck = (size_t) pResult->contigHits_ < minContigHits_;

    bool fail = (failChi2Filter || 
                 failLengthRatio ||
                 failContigMissRate ||
                 failContigHitsCheck);

    #if FILTER_DEBUG > 0
    std::cout << "MatchResult: " << *pResult << "\nt"
              << "FailChi2: " << failChi2Filter
              << " FailLengthRatio: " << failLengthRatio
              << " FailContigMissRate: " << failContigMissRate
              << " FailContigHitsCheck: " << failContigHitsCheck
              << std::endl;
    #endif

    return (!fail);
}
