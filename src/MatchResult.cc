#include <algorithm>

#include "MatchResult.h"
#include "dp.h"
#include "ContigMapData.h"
#include "OpticalMapData.h"

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
    const int m = pScoreMatrix->m_;
    const int n = pScoreMatrix->n_;
    const ScoreElement_t * pLastElement = &(pScoreMatrix->d_[end_index.first*n + end_index.second]);

    score_ = pLastElement->score_; 

    // Use the ScoreMatrix to determine the alignment from end_index
    computeAlignment(end_index, pScoreMatrix, pContigMapData, pOpticalMapData);
}

// Build the matchedFragList_  from a trail through the dynamic programming score table.
void MatchResult::computeAlignment(const Index_t& end_index, const ScoreMatrix_t * pScoreMatrix,
                          const ContigMapData * pContigMapData,
                          const OpticalMapData * pOpticalMapData)
{

    int optIndex, contigIndex;
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
    while (1)
    {
        pE = &(pScoreMatrix->d_[n*ind.first + ind.second]);
        ind = make_pair(pE->pi_, pE->pj_);
        if (ind.first >= 0 && ind.second >=0)
            trail.push_back(ind);
        else
            break;
    }
    reverse(trail.begin(), trail.end());
    const vector<Index_t>::iterator tb = trail.begin();
    const vector<Index_t>::iterator te = trail.end();
    assert(tb->first == 0); // Trail should start in first row, since we are aligning entire contig.
    opStartIndex_ = tb->second; // index of first aligned fragment in optical map
    opEndIndex_ = end_index.second-1; // index of last aligned fragment in optical map (inclusive)
    assert (opStartIndex_ >= 0);
    assert (opEndIndex_ >= 0);

    // Build the vector of matched fragments
    matchedFragList_.clear();
    matchedFragList_.reserve(trail.size());
    int pc, po; // Contig and optical end index (exclusive) of predecessor matched fragment
    int os, oe, cs, ce; // optical and contig start index (inclusive and end index (exclusive) for matched fragment
    MatchedFrag block;
    pc = 0;
    po = opStartIndex_;
    for(ti = tb+1; ti != te; ti++)
    {
        // Add this matched fragment to the matchedFragList
        // Reminder: matched fragment are given by the range of indexes
        // from opStart inclusive to opEnd exclusive (ala python slicing)
        cs = pc; // contig start index (inclusive, 0 based)
        ce = ti->first; // contig end index (exclusive)
        os = po; // optical start index (inclusive, 0 based)
        oe = ti->second; // optical end index (exclusive, 0 based)
        block = MatchedFrag(os, oe, opticalFrags, cs, ce, contigFrags);
        matchedFragList_.push_back(block);
        pc = ce; po = oe;
    }

    // Record the starting and ending position in the Optical Map.
    // This value is refined in buildAlignmentAttributes
    opStartBp_ = pOpticalMapData->getEndBp(opStartIndex_);
    opEndBp_ = pOpticalMapData->getStartBp(opEndIndex_);
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

    bool lastAlignment, firstAlignment;
    bool firstAlignmentIsGap = false;
    vector<FragData>::iterator pFrag;
    bool firstMatchedSite = true;
    int  cLeftOverhang=0; // contig size of first aligned chunk
    int  cRightOverhang=0; // contig size of last aligned chunk

    // Loop over the vector of matched fragment blocks and
    // compute alignment statistics and descriptive strings.
    vector<MatchedFrag>::const_iterator mi,mb,me;
    mb = matchedFragList_.begin();
    me = matchedFragList_.end();
    for (mi = mb; mi != me; mi++)
    {
        firstAlignment = (mi == mb);
        lastAlignment = (mi == me-1);

        contigSize_ += mi->cLength_;

        if (mi->contigGap_)
        {
            contigUnalignedFrags_++;
            contigUnalignedBases_ += mi->cLength_;
            if (firstAlignment)
                firstAlignmentIsGap = true;

        } else {

            opticalMisses_ += mi->opMisses_;
            contigMisses_ += mi->cMisses_;
            opticalTotalAlignedBases_ += mi->opLength_;
            contigTotalAlignedBases_ += mi->cLength_;
            if (!firstAlignment && !lastAlignment)
            {
                // Model assumes that Optical Frag ~ N(C, SIGMA^2*C) where C is contig frag size
                double var = Constants::SIGMA2 * mi->cLength_;
                double d = mi->opLength_ - mi->cLength_;
                chi2_ += d*d/var;
                contigInnerAlignedBases_ += mi->cLength_;
                opticalInnerAlignedBases_ += mi->opLength_;
                numAlignedInnerBlocks_++;
            }

            if (firstMatchedSite)
            {
                firstMatchedSite = false;
                cLeftOverhang = mi->cLength_;
                if (firstAlignmentIsGap)
                {
                    // The first truely aligned contig fragment is bounded by two restriction sites
                    // since one or more of the first contig fragments are aligned to gap.
                    // Account for the first restriction site here.
                    contigHits_ += 1;
                    opticalHits_ += 1;
                }
            }

            cRightOverhang = mi->cLength_;

            if (!lastAlignment)
            {
                contigHits_ += 1;
                opticalHits_ += 1;
            }
        }
    }
    
    totalHits_ = contigHits_ + opticalHits_;
    totalMisses_ = (contigMisses_ + opticalMisses_);
    totalMissRate_ = (totalMisses_ == 0) ? 0.0 : 1.0*totalMisses_/(totalHits_ + totalMisses_);
    contigMissRate_ = (contigMisses_== 0) ? 0.0 : 1.0*contigMisses_/(contigHits_ + contigMisses_);
    opticalMissRate_ = (opticalMisses_ == 0) ? 0.0 : 1.0*opticalMisses_/(opticalHits_ + opticalMisses_);
    alignedLengthRatio_ = 1.0*min(contigInnerAlignedBases_, opticalInnerAlignedBases_) /
                              max(contigInnerAlignedBases_, opticalInnerAlignedBases_);
    contigUnalignedBaseRatio_ = 1.0*contigUnalignedBases_/contigSize_;

    // Calculating the contig length before the first matched restriction site
    // and after the last matched restriction site.
    // Then calculate the approximate starting and ending bp position of the alignment.
    // from the beginning of the optical map
    opStartBp_ = opStartBp_-cLeftOverhang;
    opEndBp_ = opEndBp_+cRightOverhang;

    builtAlignmentAttributes_ = true;
}

// Populate the descriptive strings for the alignment
void MatchResult::annotate()
{
    bool lastAlignment, firstAlignment;
    bool firstAlignedFrag = true;
    vector<FragData>::const_iterator pFrag;
    stringstream opt_ss, c_ss; // Strings representing the overall alignment
    stringstream opt_aligned, c_aligned; // Strings representing indices of aligned fragments
    stringstream c_lost; // String representing the lost fragment indices

    // Loop over the vector of matched fragment blocks and
    // compute alignment statistics and descriptive strings.
    vector<MatchedFrag>::const_iterator mi,mb,me;
    mb = matchedFragList_.begin();
    me = matchedFragList_.end();
    for (mi = mb; mi != me; mi++)
    {
        lastAlignment = (mi == me-1);
        firstAlignment = (mi == mb);

        if (mi->contigGap_)
        {
            pFrag = mi->cFrags_.begin();
            opt_ss <<  " --- "; // Print a gap in the optical sequence
            c_ss << " (" << mi->cLength_ << ") ";
            c_lost << (c_lost.str().size() > 0 ? "," : "") << mi->cStart_;

        } else {

            vector<FragData>::const_iterator it,itb,ite;

            if (firstAlignedFrag && !firstAlignment)
            {
                // This means that one or more fragments at the beginning of contig
                // are in a gapped alignment. Therefore, this fragment is bounded by
                // two contig restriction sites.

                // Add the first restriction site to the annotation
                opt_ss << "; ";
                opt_aligned << "; ";
                c_ss << "; ";
                c_aligned << "; ";
            }

            firstAlignedFrag = false;

            // Add to alignment descriptive strings
            ite = mi->opFrags_.end();
            for (it = mi->opFrags_.begin(); it!=ite; it++)
                opt_ss << " " <<  it->size_ << "," << it->sd_ << " ";

            int lasti = mi->opEnd_-1;
            for (int i = mi->opStart_; i <= lasti; i++)
                opt_aligned << i << ( i==lasti  ? "" : ",");

            ite = mi->cFrags_.end();
            for (it = mi->cFrags_.begin(); it!=ite; it++)
                c_ss << " " <<  it->size_ << "," << " ";

            lasti = mi->cEnd_-1;
            for (int i = mi->cStart_; i <= lasti; i++)
                c_aligned << i << ( i==lasti  ? "" : ",");

            if (!lastAlignment)
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
    matchedFragList_.clear();

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
