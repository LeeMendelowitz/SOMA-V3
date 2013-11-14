#include "StandardMatchMaker.h"
#include "MatchResult.h"
#include "MatchedChunk.h"
#include "debugUtils.h"

#include <algorithm>
#include <iostream>

#define BUILDMATCH_DEBUG 0
#define FILTER_DEBUG 0
#define MAKEMATCH_DEBUG 0

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



MatchResult * StandardMatchMaker::buildMatch(const Index_t& end_index, const ScoreMatrix_t * pScoreMatrix, const MapData * pOpticalMap,
                                             const ContigMapData * pContigMap, bool contigIsForward)
{

    const vector<FragData>& contigFrags = pContigMap->getFrags();
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

    // Use the id of forward contig map
    const ContigMapData * pContigMapForward = contigIsForward ? pContigMap : pContigMap->getTwin();
    MatchResult * pMatch = new MatchResult(pContigMapForward->getId(), pOpticalMap->getId(),
                                           pContigMap->getLength(), contigIsForward, score);

    #if BUILDMATCH_DEBUG > 0
    std::cout << "\n\n\nBuilding match for: "
              << pContigMap->getId() << " "
              << pOpticalMap->getId() << " "
              << " forward: " << contigIsForward
              << " end_index: " << end_index.first << " , " << end_index.second
              << std::endl;

    // Print the contig fragments
    std::cout << "Contig Frags: " ;
    for(vector<FragData>::const_iterator iter = contigFrags.begin();
                                         iter != contigFrags.end();
                                         iter++)
    { std::cout << iter->size_ << " "; }
    std::cout << std::endl;
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

        int cStartBp = pContigMap->getStartBp(cs);
        int cEndBp = pContigMap->getEndBp(ce-1); // substract 1 since end is exclusive
        bool gapAlignment = (oe == os);
        int opStartBp, opEndBp;

        #if BUILDMATCH_DEBUG > 0
        std::cout << "Contig: "
                  << "(" << cs << "," << ce << ")"
                  << " = " << "(" << cStartBp << "," << cEndBp << "=" << cEndBp - cStartBp << ")"
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

        #if BUILDMATCH_DEBUG > 0
        std::cout << "Optical: "
                  << "(" << os << "," << oe << ")"
                  << " = " << "(" << opStartBp << "," << opEndBp << "=" << opEndBp - opStartBp <<  ")"
                  << std::endl;
        #endif

        // For boundary cases where the matched chunk includes the first or last contig fragment,
        // the chunk is not bounded by two matched restriction sites. Adjust the optical start/end positions
        // accordingly.
        if (opt::useBoundaries)
        {
            if (cs == 0)
            {
                // Adjust the optical start position
                opStartBp = opEndBp - (cEndBp - cStartBp);
                #if BUILDMATCH_DEBUG > 0
                std::cout << "cs==0 adjustment: "
                          << "(" << os << "," << oe << ")"
                          << " = " << "(" << opStartBp << "," << opEndBp << ")"
                          << std::endl;
                #endif
            }
            else if ((size_t) ce == contigFrags.size())
            {
                opEndBp = opStartBp + (cEndBp - cStartBp);
                #if BUILDMATCH_DEBUG > 0
                std::cout << "ce==contigFrags.size() adjustment: "
                          << "(" << os << "," << oe << ")"
                          << " = " << "(" << opStartBp << "," << opEndBp << ")"
                          << std::endl;
                #endif
            }
        }

        bool isBoundaryChunk = opt::useBoundaries && ( (cs==0) || ( (size_t) ce == contigFrags.size()) );

        MatchedChunk chunk = MatchedChunk(os, oe, opStartBp, opEndBp, pOpticalMap,
                                          cs, ce, cStartBp, cEndBp, pContigMap, isBoundaryChunk);

        ////////////////////////////////////////////////////////////////////////////
        #if BUILDMATCH_DEBUG > 0
        {
            std::cout << "Built Chunk:\n"
                      << "\tcontig = (" << chunk.getContigStartIndex() << "," << chunk.getContigEndIndex() << "): ";

            // Add the contig frag lengths
            int cSize = 0;
            for(vector<FragData>::const_iterator iter = chunk.getContigFragB(); 
                                                 iter != chunk.getContigFragE();
                                                 iter++
                                                 )
            {
                std::cout  << iter->size_ << " ";
                if (iter != chunk.getContigFragE()-1)
                    std::cout << "+ ";
                cSize += iter->size_;
            }
            std::cout << " = " << cSize << "\n";

            std::cout << "\toptical = (" << chunk.getOpticalStartIndex() << "," << chunk.getOpticalEndIndex() << ")";

            // Add the optical frag lengths
            int oSize = 0;
            for(vector<FragData>::const_iterator iter = chunk.getOpticalFragB(); 
                                                 iter != chunk.getOpticalFragE();
                                                 iter++
                                                 )
            {
                std::cout  << iter->size_ << " ";
                if (iter != chunk.getOpticalFragE()-1)
                    std::cout << "+ ";
                oSize += iter->size_;
            }
            std::cout << " = " << oSize << std::endl;
        }
        #endif
        ////////////////////////////////////////////////////////////////////////////


        matchedChunkList.push_back(chunk);
        pc = ce; po = oe;
    }

    pMatch->buildAlignmentAttributes();

    return pMatch;
}

// Return true if MatchResult is acceptable
bool StandardMatchMaker::filterFunction(const MatchResult * pResult)
{

    //Chi2 Filter
    bool failChi2Filter = false;
    if (pResult->numAlignedInnerBlocks_ > 0)
        failChi2Filter = (pResult->chi2_/((double) pResult->numAlignedInnerBlocks_)) > avgChi2Threshold_;
       
    bool failLengthRatio = pResult->alignedLengthRatio_ < minLengthRatio_;
    bool failContigMissRate = pResult->contigMissRate_ > maxMissRateContig_;
    bool failContigHitsCheck = (size_t) pResult->contigHits_ <= minContigHits_;

    bool fail = (failChi2Filter || 
                 failLengthRatio ||
                 failContigMissRate ||
                 failContigHitsCheck);

    #if FILTER_DEBUG > 0
    std::cout << "FILTER: MatchResult=" << *pResult << "\n"
              << "FailChi2: " << failChi2Filter
              << " FailLengthRatio: " << failLengthRatio
              << " FailContigMissRate: " << failContigMissRate
              << " FailContigHitsCheck: " << failContigHitsCheck << " ";
    printAttributes(std::cout, *pResult);
    std::cout << "\n\n\n" << std::endl;
    #endif

    return (!fail);
}

// Build a match result from the scoreMatrix which ends at end_index
MatchResult * StandardMatchMaker::buildMatch(const Index_t& end_index, const seeded::ScoreMatrix& scoreMatrix, const MapData * pOpticalMap,
                                             const ContigMapData * pContigMap, bool contigIsForward)
{

    const vector<FragData>& contigFrags = pContigMap->getFrags();
    const int n = scoreMatrix.getNumCols() ;
    const ScoreCell * pCell = scoreMatrix.getCell(end_index.first, end_index.second);
    double score = pCell->score_;

    assert(pCell->backPointer_.getSource() == pCell);

    // Use the id of forward contig map
    const ContigMapData * pContigMapForward = contigIsForward ? pContigMap : pContigMap->getTwin();
    MatchResult * pMatch = new MatchResult(pContigMapForward->getId(), pOpticalMap->getId(),
                                           pContigMap->getLength(), contigIsForward, score);

    // Build the trail through scoreMatrix
    vector<const ScoreCell*> trail;
    trail.reserve(scoreMatrix.getNumRows());
    trail.push_back(pCell);
    while(true)
    {
        #if BUILDMATCH_DEBUG > 0
        std::cout << "pCell: " << *pCell << endl;
        #endif
        ScoreCell * pNext = pCell->backPointer_.getTarget();
        if(!pNext) break;
        assert(pCell->backPointer_.getSource() == pCell);
        trail.push_back(pNext);
        pCell = pNext;
    }
    reverse(trail.begin(), trail.end());

    const vector<const ScoreCell*>::iterator tb = trail.begin();
    const vector<const ScoreCell*>::iterator te = trail.end();

    #if BUILDMATCH_DEBUG > 0
    std::cout << "\n\n\nBuilding match for: "
              << pContigMap->getId() << " "
              << pOpticalMap->getId() << " "
              << " forward: " << contigIsForward
              << " end_index: " << end_index.first << " , " << end_index.second
              << std::endl;

    // Print the contig fragments
    std::cout << "Contig Frags: " ;
    for(vector<FragData>::const_iterator iter = contigFrags.begin();
                                         iter != contigFrags.end();
                                         iter++)
    { std::cout << iter->size_ << " "; }
    std::cout << std::endl;
    # endif

    // Build the vector of matched fragments. Add them directory to the match.
    vector<MatchedChunk>& matchedChunkList = pMatch->matchedChunkList_;
    matchedChunkList.reserve(trail.size());
    int pc = (*tb)->q_; // query end index (exclusive) of predecessor matched fragment
    int po = (*tb)->r_; // ref end index (exclusive) of predecessor matched fragment
    int os, oe, cs, ce; // optical and contig start index inclusive and end index (exclusive) for matched fragment
    for(auto ti = tb+1; ti != te; ti++)
    {
        // Add this matched fragment to the matchedChunkList
        // Reminder: matched fragment are given by the range of indexes
        // from opStart (inclusive) to opEnd (exclusive)

        // Here, indices refer to the Fragment index (not the site index)
        // ------|----------|-----------|-------- ...
        //   0         1           2         3    ...

        cs = pc; // contig start index (inclusive, 0 based)
        ce = (*ti)->q_; // contig end index (exclusive)
        os = po; // optical start index (inclusive, 0 based)
        oe = (*ti)->r_; // optical end index (exclusive, 0 based)

        int cStartBp = pContigMap->getStartBp(cs);
        int cEndBp = pContigMap->getEndBp(ce-1); // substract 1 since end is exclusive
        bool gapAlignment = (oe == os);
        int opStartBp, opEndBp;

        #if BUILDMATCH_DEBUG > 0
        std::cout << "Contig: "
                  << "(" << cs << "," << ce << ")"
                  << " = " << "(" << cStartBp << "," << cEndBp << "=" << cEndBp - cStartBp << ")"
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

        #if BUILDMATCH_DEBUG > 0
        std::cout << "Optical: "
                  << "(" << os << "," << oe << ")"
                  << " = " << "(" << opStartBp << "," << opEndBp << "=" << opEndBp - opStartBp <<  ")"
                  << std::endl;
        #endif

        // For boundary cases where the matched chunk includes the first or last contig fragment,
        // the chunk is not bounded by two matched restriction sites. Adjust the optical start/end positions
        // accordingly.
        if (opt::useBoundaries)
        {
            if (cs == 0)
            {
                // Adjust the optical start position
                opStartBp = opEndBp - (cEndBp - cStartBp);
                #if BUILDMATCH_DEBUG > 0
                std::cout << "cs==0 adjustment: "
                          << "(" << os << "," << oe << ")"
                          << " = " << "(" << opStartBp << "," << opEndBp << ")"
                          << std::endl;
                #endif
            }
            else if ((size_t) ce == contigFrags.size())
            {
                opEndBp = opStartBp + (cEndBp - cStartBp);
                #if BUILDMATCH_DEBUG > 0
                std::cout << "ce==contigFrags.size() adjustment: "
                          << "(" << os << "," << oe << ")"
                          << " = " << "(" << opStartBp << "," << opEndBp << ")"
                          << std::endl;
                #endif
            }
        }

        bool isBoundaryChunk = opt::useBoundaries && ( (cs==0) || ( (size_t) ce == contigFrags.size()) );

        MatchedChunk chunk = MatchedChunk(os, oe, opStartBp, opEndBp, pOpticalMap,
                                          cs, ce, cStartBp, cEndBp, pContigMap, isBoundaryChunk);

        ////////////////////////////////////////////////////////////////////////////
        #if BUILDMATCH_DEBUG > 0
        {
            std::cout << "Built Chunk:\n"
                      << "\tcontig = (" << chunk.getContigStartIndex() << "," << chunk.getContigEndIndex() << "): ";

            // Add the contig frag lengths
            int cSize = 0;
            for(vector<FragData>::const_iterator iter = chunk.getContigFragB(); 
                                                 iter != chunk.getContigFragE();
                                                 iter++
                                                 )
            {
                std::cout  << iter->size_ << " ";
                if (iter != chunk.getContigFragE()-1)
                    std::cout << "+ ";
                cSize += iter->size_;
            }
            std::cout << " = " << cSize << "\n";

            std::cout << "\toptical = (" << chunk.getOpticalStartIndex() << "," << chunk.getOpticalEndIndex() << ")";

            // Add the optical frag lengths
            int oSize = 0;
            for(vector<FragData>::const_iterator iter = chunk.getOpticalFragB(); 
                                                 iter != chunk.getOpticalFragE();
                                                 iter++
                                                 )
            {
                std::cout  << iter->size_ << " ";
                if (iter != chunk.getOpticalFragE()-1)
                    std::cout << "+ ";
                oSize += iter->size_;
            }
            std::cout << " = " << oSize << std::endl;
        }
        #endif
        ////////////////////////////////////////////////////////////////////////////


        matchedChunkList.push_back(chunk);
        pc = ce; po = oe;
    }

    pMatch->buildAlignmentAttributes();

    return pMatch;
}
