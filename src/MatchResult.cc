#include <algorithm>

#include "MatchResult.h"
#include "dp.h"
#include "ContigMapData.h"
#include "OpticalMapData.h"

#include <iostream>

#define MR_DEBUG 0

// This is the constructor for MatchResult
// The MatchResult data is determined from the scoreMatrix and the end_index
// and walking backwards through the scoreMatrix.
MatchResult::MatchResult(const Index_t& end_index, const ScoreMatrix_t * pScoreMatrix,
                         const ContigMapData * pContigMapData, const OpticalMapData * pOpticalMapData,
                         bool forward)
{
    reset();
    chromosomeId_ = pOpticalMapData->opticalId_;
    contigId_ = pContigMapData->contigId_;
    forward_ = forward;
    const int n = pScoreMatrix->n_;
    const ScoreElement_t * pLastElement = &(pScoreMatrix->d_[end_index.first*n + end_index.second]);

    score_ = pLastElement->score_; 

    // Use the ScoreMatrix to determine the alignment from end_index
    computeAlignment(end_index, pScoreMatrix, pContigMapData, pOpticalMapData);
}

// Build the matchedChunkList_  from a trail through the dynamic programming score table.
// NOTE: Instead of calling computeAlignment, a MatchResult could easily be handed a trail of MatchedChunks
// and the MatchResultAlignment could be computed from this. This would allow for much more flexible
// construction of MatchResults from the dynamic programming table, if needed.
void MatchResult::computeAlignment(const Index_t& end_index, const ScoreMatrix_t * pScoreMatrix,
                          const ContigMapData * pContigMapData,
                          const OpticalMapData * pOpticalMapData)
{

    const int m = pScoreMatrix->m_;
    const int n = pScoreMatrix->n_;
    ScoreElement_t * pE;
    vector<Index_t> trail; // Trail through dynamic programming score table for the best alignment
    vector<Index_t>::iterator ti;
    Index_t ind = end_index;
    const vector<FragData>& contigFrags = forward_ ? pContigMapData->frags_ : pContigMapData->reverseFrags_;
    const vector<FragData>& opticalFrags = pOpticalMapData->frags_;

    // Build the trail through scoreMatrix
    trail.reserve(m);
    trail.push_back(ind);
    while (true)
    {
        pE = &(pScoreMatrix->d_[n*ind.first + ind.second]);
        ind = make_pair(pE->pi_, pE->pj_);
        if ((ind.first < 0) || (ind.second < 0) )
            break;
        trail.push_back(ind);
    }
    reverse(trail.begin(), trail.end());

    #if MR_DEBUG > 0
    std::cout << "Trail size: " << trail.size() << std::endl;
    #endif

    const vector<Index_t>::iterator tb = trail.begin();
    const vector<Index_t>::iterator te = trail.end();
    cStartIndex_ = tb->first;
    opStartIndex_ = tb->second; // index of first aligned fragment in optical map
    // Subtract 1 from the end_index of the dynamic programming score matrix
    // because the dynamic programming score matrix has an inserted first row/column.
    cEndIndex_ = end_index.first-1;
    opEndIndex_ = end_index.second-1;

    assert (opStartIndex_ >= 0);
    assert (opEndIndex_ >= 0);

    // Build the vector of matched fragments
    matchedChunkList_.clear();
    matchedChunkList_.reserve(trail.size());
    int pc, po; // Contig and optical end index (exclusive) of predecessor matched fragment
    int os, oe, cs, ce; // optical and contig start index inclusive and end index (exclusive) for matched fragment
    pc = cStartIndex_;
    po = opStartIndex_;
    for(ti = tb+1; ti != te; ti++)
    {
        // Add this matched fragment to the matchedChunkList
        // Reminder: matched fragment are given by the range of indexes
        // from opStart inclusive to opEnd exclusive (ala python slicing)
        cs = pc; // contig start index (inclusive, 0 based)
        ce = ti->first; // contig end index (exclusive)
        os = po; // optical start index (inclusive, 0 based)
        oe = ti->second; // optical end index (exclusive, 0 based)

        int cStartBp = pContigMapData->getStartBp(cs, forward_);
        int cEndBp = pContigMapData->getEndBp(ce, forward_);
        int opStartBp = pOpticalMapData->getStartBp(os);
        int opEndBp = pOpticalMapData->getEndBp(os);

        // For boundary cases where the matched chunk includes the first or last contig fragment,
        // the chunk is not bounded by two matched restriction sites. Adjust the optical start/end positions
        // accordingly.
        if (cs == 0)
        {
            // Adjust the optical start position
            opStartBp = opEndBp - (cEndBp - cStartBp);
        }
        else if ((size_t) ce == contigFrags.size())
        {
            opEndBp = opStartBp + (cEndBp - cStartBp);
        }

        MatchedChunk chunk = MatchedChunk(os, oe, opStartBp, opEndBp, opticalFrags,
                                          cs, ce, cStartBp, cEndBp, contigFrags);
        matchedChunkList_.push_back(chunk);
        pc = ce; po = oe;
    }

    //////////////////////////////////////////////////////
    //#ifdef DEBUG
    // Check that there are at most two boundary chunks.
    int boundaryCount = 0;
    const vector<MatchedChunk>::const_iterator E = matchedChunkList_.end();
    for (vector<MatchedChunk>::const_iterator iter = matchedChunkList_.begin();
         iter != E;
         iter++)
    {
        if (iter->boundaryChunk_)
        {
            boundaryCount++;
            assert( iter == matchedChunkList_.begin() || iter == matchedChunkList_.end());
        }
    }
    assert(boundaryCount <= 2);
    //#endif
    ////////////////////////////////////////////////////////
}


ostream& operator<<(ostream& os, const MatchResult& mr)
{
    os << mr.contigId_ << " " << mr.contigSize_ << " " << mr.forward_ << " "  << mr.opStartIndex_ << " " << mr.opEndIndex_ << "\n"
                 << mr.totalMisses_ << " " << " " << mr.pval_ << "\n"
                 << mr.opticalMatchString_ << "\n" << mr.contigMatchString_ << "\n";
    return os;
}


void MatchResult::buildAlignmentAttributes()
{
    if (builtAlignmentAttributes_)
        return;


    if (matchedChunkList_.empty()) return;

    // Assign the starting and ending bp locations.
    const MatchedChunk& firstChunk = matchedChunkList_.front();
    opStartBp_ = firstChunk.opStartBp_;
    cStartBp_ = firstChunk.cStartBp_;

    const MatchedChunk& lastChunk = matchedChunkList_.back();
    opEndBp_ = lastChunk.opEndBp_;
    cEndBp_ = lastChunk.cEndBp_;

    vector<FragData>::iterator pFrag;

    // Loop over the vector of matched fragment chunks and
    // compute alignment statistics and descriptive strings.
    typedef vector<MatchedChunk>::const_iterator ChunkIter;
    const ChunkIter mb = matchedChunkList_.begin();
    const ChunkIter me = matchedChunkList_.end();
    for (ChunkIter mi = mb; mi != me; mi++)
    {
        bool firstAlignment = (mi == mb);

        contigSize_ += mi->cLength_;

        if (mi->contigGap_)
        {
            contigUnalignedFrags_++;
            contigUnalignedBases_ += mi->cLength_;
            continue;
        } 

        opticalMisses_ += mi->opMisses_;
        contigMisses_ += mi->cMisses_;
        opticalTotalAlignedBases_ += mi->opLength_;
        contigTotalAlignedBases_ += mi->cLength_;
        if (!mi->boundaryChunk_)
        {
            // Model assumes that Optical Frag ~ N(C, SIGMA^2*C) where C is contig frag size
            double var = Constants::SIGMA2 * mi->cLength_;
            double d = mi->opLength_ - mi->cLength_;
            chi2_ += d*d/var;
            contigInnerAlignedBases_ += mi->cLength_;
            opticalInnerAlignedBases_ += mi->opLength_;
            numAlignedInnerBlocks_++;
        }

        if (firstAlignment && !mi->boundaryChunk_)
        {
            // The first alignment is not a boundary chunk.
            // Therefore, it is bounded by two restriction sites.
            // Add the contribution of the left matched sites to the hit count
            contigHits_ += 1;
            opticalHits_ += 1;
        }

        contigHits_ += 1;
        opticalHits_ += 1;
    }
    
    totalHits_ = contigHits_ + opticalHits_;
    totalMisses_ = (contigMisses_ + opticalMisses_);
    totalMissRate_ = (totalMisses_ == 0) ? 0.0 : 1.0*totalMisses_/(totalHits_ + totalMisses_);
    contigMissRate_ = (contigMisses_== 0) ? 0.0 : 1.0*contigMisses_/(contigHits_ + contigMisses_);
    opticalMissRate_ = (opticalMisses_ == 0) ? 0.0 : 1.0*opticalMisses_/(opticalHits_ + opticalMisses_);
    alignedLengthRatio_ = 1.0*min(contigInnerAlignedBases_, opticalInnerAlignedBases_) /
                              max(contigInnerAlignedBases_, opticalInnerAlignedBases_);
    contigUnalignedBaseRatio_ = 1.0*contigUnalignedBases_/contigSize_;

    builtAlignmentAttributes_ = true;
}

// Populate the descriptive strings for the alignment
void MatchResult::annotate()
{
    bool lastAlignment;
    //bool firstAlignment;
    bool firstAlignedFrag = true;
    stringstream opt_ss, c_ss; // Strings representing the overall alignment
    stringstream opt_aligned, c_aligned; // Strings representing indices of aligned fragments
    stringstream c_lost; // String representing the lost fragment indices

    // Loop over the vector of matched fragment chunks and
    // compute alignment statistics and descriptive strings.
    vector<MatchedChunk>::const_iterator mi,mb,me;
    mb = matchedChunkList_.begin();
    me = matchedChunkList_.end();
    for (mi = mb; mi != me; mi++)
    {
        lastAlignment = (mi == me-1);
        //firstAlignment = (mi == mb);

        if (mi->contigGap_)
        {
            opt_ss <<  " --- "; // Print a gap in the optical sequence
            c_ss << " (" << mi->cLength_ << ") ";
            c_lost << (c_lost.str().size() > 0 ? "," : "") << mi->cStart_;

        } else {

            vector<FragData>::const_iterator it,itb,ite;

            if (firstAlignedFrag && !mi->boundaryChunk_)
            {
                // This chunk is bounded by a restriction site
                // Semi-colon denotes a matched restriction site
                opt_ss << "; ";
                opt_aligned << "; ";
                c_ss << "; ";
                c_aligned << "; ";
            }

            firstAlignedFrag = false;

            // Add to alignment descriptive strings
            for (it = mi->opFragB_; it != mi->opFragE_; it++)
                opt_ss << " " <<  it->size_ << "," << it->sd_ << " ";

            int lasti = mi->opEnd_-1;
            for (int i = mi->opStart_; i <= lasti; i++)
                opt_aligned << i << ( i==lasti  ? "" : ",");

            for (it = mi->cFragB_; it != mi->cFragE_; it++)
                c_ss << " " <<  it->size_ << "," << " ";

            lasti = mi->cEnd_-1;
            for (int i = mi->cStart_; i <= lasti; i++)
                c_aligned << i << ( i==lasti  ? "" : ",");

            if (!lastAlignment || !mi->boundaryChunk_)
            {
                // Semi-colon denotes a matched restriction site
                c_ss << "; ";
                c_aligned << ";";
                opt_ss << "; ";
                opt_aligned << ";";
            }
        }
    }

    contigAlignedIndexString_ = c_aligned.str();
    contigMatchString_ = c_ss.str();
    contigLostIndexString_ = c_lost.str();
    opticalAlignedIndexString_ = opt_aligned.str();
    opticalMatchString_ = opt_ss.str();
}


// Reset all attributes
void MatchResult::reset()
{
    builtAlignmentAttributes_ = false;

    // Contig information
    contigId_.clear();
    contigSize_ = 0;

    // Alignment location & orientation
    chromosomeId_.clear();
    opStartIndex_ = 0;
    opEndIndex_ = 0;
    opStartBp_ = 0;
    opEndBp_ = 0;
    forward_ = 0;
    matchedChunkList_.clear();

    // Alignment statistics
    score_ = 0;
    totalMisses_ = 0;
    totalHits_ = 0;
    totalMissRate_ = 0;
    numAlignedInnerBlocks_ = 0;
    contigMisses_ = 0;
    contigHits_ = 0;
    contigMissRate_ = 0;
    contigUnalignedFrags_ = 0;
    contigUnalignedBases_ = 0;
    contigTotalAlignedBases_ = 0;
    contigInnerAlignedBases_ = 0;
    opticalMisses_ = 0;
    opticalHits_ = 0;
    opticalMissRate_ = 0;
    opticalTotalAlignedBases_ = 0;
    opticalInnerAlignedBases_ = 0;
    alignedLengthRatio_ = 0;
    contigUnalignedBaseRatio_ = 0;
    chi2_ = 0.0;
    pval_ = 0;

    // Alignment Descriptive strings
    contigMatchString_.clear();
    contigAlignedIndexString_.clear();
    contigLostIndexString_.clear();
    opticalMatchString_.clear();
    opticalAlignedIndexString_.clear();
}
