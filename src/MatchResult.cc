#include <algorithm>

#include "MatchResult.h"
#include "dp.h"
#include "ContigMapData.h"
#include "OpticalMapData.h"

#include <iostream>

#define MR_DEBUG 0

void MatchResult::setMatchedChunks(const MatchedChunkVec& matchedChunkList)
{
    matchedChunkList_ = matchedChunkList;
    buildAlignmentAttributes();
}

ostream& operator<<(ostream& os, const MatchResult& mr)
{
    os << mr.contigId_ << " "
       << mr.chromosomeId_ << " "
       << mr.contigLength_ << " "
       << mr.forward_ << " "
       << mr.opStartIndex_ << " "
       << mr.opEndIndex_ << " "
       << mr.cStartIndex_ << " "
       << mr.cEndIndex_ << " "
       << mr.totalMisses_ << " ";
       //<< mr.pval_;
    return os;
}


void MatchResult::buildAlignmentAttributes()
{

    if (matchedChunkList_.empty()) return;

    contigUnalignedFrags_ = 0;
    opticalMisses_ = 0;
    contigMisses_ = 0;
    contigHits_  = 0;
    opticalHits_ = 0;
    opticalTotalAlignedBases_ = 0;
    contigTotalAlignedBases_ = 0;
    contigInnerAlignedBases_ = 0;
    opticalInnerAlignedBases_ = 0;
    contigUnalignedBases_ = 0;
    numAlignedInnerBlocks_ = 0;

    // Assign the starting and ending bp locations.
    const MatchedChunk& firstChunk = matchedChunkList_.front();
    opStartIndex_ = firstChunk.opStart_;
    cStartIndex_ = firstChunk.cStart_;
    opStartBp_ = firstChunk.opStartBp_;
    cStartBp_ = firstChunk.cStartBp_;


    const MatchedChunk& lastChunk = matchedChunkList_.back();
    opEndIndex_ = lastChunk.opEnd_ - 1; // inclusive
    cEndIndex_ = lastChunk.cEnd_ - 1; // inclusive
    opEndBp_ = lastChunk.opEndBp_;
    cEndBp_ = lastChunk.cEndBp_;

    // Loop over the vector of matched fragment chunks and
    // compute alignment statistics and descriptive strings.
    typedef vector<MatchedChunk>::const_iterator ChunkIter;
    const ChunkIter mb = matchedChunkList_.begin();
    const ChunkIter me = matchedChunkList_.end();


    for (ChunkIter mi = mb; mi != me; mi++)
    {
        bool firstAlignment = (mi == mb);

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
    contigUnalignedBaseRatio_ = 1.0*contigUnalignedBases_/contigLength_;

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
    stringstream score_ss;

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

    // Make a string depicting the components of the 
    // score for a matched chunk
    score_ss.precision(3);
    for (vector<Score>::const_iterator iter = scoreList_.begin();
         iter != scoreList_.end();
         iter++)
    {
        score_ss << "(" << iter->contig << "," << iter->optical << "," << iter->sizing << ")";
        if (iter != scoreList_.end() - 1)
            score_ss << " / ";
    }


    contigAlignedIndexString_ = c_aligned.str();
    contigMatchString_ = c_ss.str();
    contigLostIndexString_ = c_lost.str();
    opticalAlignedIndexString_ = opt_aligned.str();
    opticalMatchString_ = opt_ss.str();
    scoreString_ = score_ss.str();
}


// Reset all attributes
void MatchResult::reset()
{

    // Contig information
    contigId_.clear();

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
