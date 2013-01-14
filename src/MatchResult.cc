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
    opStartIndex_ = firstChunk.getOpticalStartIndex();
    cStartIndex_ = firstChunk.getContigStartIndex();
    opStartBp_ = firstChunk.getOpticalStartBp();
    cStartBp_ = firstChunk.getContigStartBp();


    const MatchedChunk& lastChunk = matchedChunkList_.back();
    opEndIndex_ = lastChunk.getOpticalEndIndex() - 1; // inclusive
    cEndIndex_ = lastChunk.getContigEndIndex() - 1; // inclusive
    opEndBp_ = lastChunk.getOpticalEndBp();
    cEndBp_ = lastChunk.getContigEndBp();

    // Loop over the vector of matched fragment chunks and
    // compute alignment statistics and descriptive strings.
    typedef vector<MatchedChunk>::const_iterator ChunkIter;
    const ChunkIter mb = matchedChunkList_.begin();
    const ChunkIter me = matchedChunkList_.end();


    for (ChunkIter mi = mb; mi != me; mi++)
    {
        bool firstAlignment = (mi == mb);
        int cl = mi->getContigMatchLengthBp();
        int ol = mi->getOpticalMatchLengthBp();

        if (mi->isContigGap())
        {
            contigUnalignedFrags_ += mi->getNumContigFrags();
            contigUnalignedBases_ += cl;
            continue;
        } 

        opticalMisses_ += mi->getNumOpticalMisses();
        contigMisses_ += mi->getNumContigMisses();
        opticalTotalAlignedBases_ += ol;
        contigTotalAlignedBases_ += cl;
        if (!mi->isBoundaryChunk())
        {
            // Model assumes that Optical Frag ~ N(C, SIGMA^2*C) where C is contig frag size
            double var = Constants::SIGMA2 * cl;
            double d = ol - cl;
            chi2_ += d*d/var;
            contigInnerAlignedBases_ += cl;
            opticalInnerAlignedBases_ += ol;
            numAlignedInnerBlocks_++;
        }

        if (firstAlignment && !mi->isBoundaryChunk())
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
    bool firstAlignedFrag = true;
    bool lastChunk = false;
    stringstream opt_ss, c_ss; // Strings representing the overall alignment
    stringstream opt_aligned, c_aligned; // Strings representing indices of aligned fragments
    stringstream c_lost; // String representing the lost fragment indices
    stringstream score_ss; score_ss.precision(3);

    // Loop over the vector of matched fragment chunks and
    // compute alignment statistics and descriptive strings.
    vector<MatchedChunk>::const_iterator mi,mb,me;
    mb = matchedChunkList_.begin();
    me = matchedChunkList_.end();
    for (mi = mb; mi != me; mi++)
    {
        lastChunk = (mi == me-1);

        if (mi->isContigGap())
        {
            opt_ss <<  " --- "; // Print a gap in the optical sequence
            c_ss << " (" << mi->getContigMatchLengthBp() << ") ";
            c_lost << (c_lost.str().size() > 0 ? "," : "") << mi->getContigStartIndex();

        } else {

            vector<FragData>::const_iterator it,itb,ite;

            if (firstAlignedFrag && !mi->isBoundaryChunk())
            {
                // This chunk is bounded by a restriction site
                // Semi-colon denotes a matched restriction site
                opt_ss << "; ";
                opt_aligned << "; ";
                c_ss << "; ";
                c_aligned << "; ";
            }


            // Add to alignment descriptive strings
            const vector<FragData>::const_iterator opE = mi->getOpticalFragE();
            const vector<FragData>::const_iterator cE = mi->getContigFragE();
            for (it = mi->getOpticalFragB();
                 it != opE;
                 it++)
                opt_ss << " " <<  it->size_ << "," << it->sd_ << " ";

            int lasti = mi->getOpticalEndIndex()-1;
            for (int i = mi->getOpticalStartIndex();
                 i <= lasti;
                 i++)
                opt_aligned << i << ( i==lasti  ? "" : ",");

            for (it = mi->getContigFragB();
                 it != cE;
                 it++)
                c_ss << " " <<  it->size_ << "," << " ";

            lasti = mi->getContigEndIndex()-1;
            for (int i = mi->getContigStartIndex(); i <= lasti; i++)
                c_aligned << i << ( i==lasti  ? "" : ",");

            // We post a trailing semi-colon denoting a matched site
            // if it is a any chunk other than the first
            if (firstAlignedFrag || !mi->isBoundaryChunk())
            {
                // Semi-colon denotes a matched restriction site
                c_ss << "; ";
                c_aligned << ";";
                opt_ss << "; ";
                opt_aligned << ";";
            }

            firstAlignedFrag = false;
        }

        const Score& score = mi->getScore();
        score_ss << "(" << score.contig << "," << score.optical << "," << score.sizing << ")";
        if (!lastChunk)
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
