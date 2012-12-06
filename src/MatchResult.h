#ifndef MATCHRESULT_H
#define MATCHRESULT_H

#include <numeric>

#include "ScoreMatrix.h"
#include "mapTypes.h"
#include "OpticalMapData.h"
#include "ContigMapData.h"

using namespace std;

typedef pair<int,int> Index_t;

// Represents an aligned fragment block ("chunk")
class MatchedFrag
{
    public:
        // Constructor
        MatchedFrag () {}; // Default constructor
        MatchedFrag(int opStart, int opEnd, const vector<FragData>& oData,
                    int cStart, int cEnd, const vector<FragData>& cData) :
        opStart_(opStart), opEnd_(opEnd), cStart_(cStart), cEnd_(cEnd),
        opLength_(0), cLength_(0), opMisses_(opEnd-opStart-1),
        cMisses_(cEnd-cStart-1)
        {
            vector<FragData>::const_iterator it, ite, itb;
            contigGap_ = (opStart_ == opEnd_);

            // Populate opticalFrags_
            opFrags_.resize(opMisses_+1);
            it = oData.begin();
            copy(it+opStart_, it+opEnd_, opFrags_.begin());

            // Populate contigFrags_;
            cFrags_.resize(cMisses_+1);
            it = cData.begin();
            copy(it+cStart_, it+cEnd_, cFrags_.begin());
    
            itb = oData.begin() + opStart;
            ite = oData.begin() + opEnd;
            for(it=itb; it!=ite; it++)
                opLength_ += it->size_;

            itb = cData.begin() + cStart;
            ite = cData.begin() + cEnd;
            for(it=itb; it!=ite; it++)
                cLength_ += it->size_;
        }

        // Data
        bool contigGap_;
        int opStart_; // optical start index
        int opEnd_; // optical end index (exclusive)
        int opLength_; // optical fragment length (bp)
        int cStart_; // contig start index
        int cEnd_; // contig end index (exclusive)
        int cLength_; // contig fragment length (bp)
        int opMisses_; // number of unaligned optical sites within the block
        int cMisses_; // number of unaligned contig sites within the block
        vector<FragData> cFrags_;
        vector<FragData> opFrags_;
};

//A struct for storing match result information about an alignment
class MatchResult {

    public:
    /////////////////////////////////////////////////////////////////
    // Functions Declarations

    // Default Constructor
    MatchResult()
    {
        reset();
    };

    MatchResult(const Index_t& end_index, const ScoreMatrix_t * pScoreMatrix,
                 const ContigMapData * pContigMapData, const OpticalMapData * pOpticalMapData,
                 bool forward);

    // Reset all attributes
    void reset();

    // Compute the alignment using the dynamic programming scoring table
    void computeAlignment(const Index_t& end_index, const ScoreMatrix_t * pScoreMatrix,
                          const ContigMapData * pContigMapData,
                          const OpticalMapData * pOpticalMapData);

    // Build additional alignment attributes
    void buildAlignmentAttributes();
    
    // Build descriptive strings of the alignment
    void annotate();

    bool operator==(const MatchResult& a) const {
        return score_ == a.score_;
    }

    // Return true if this MatchResult overlaps with
    // matchResult pMatch.
    bool overlaps(const MatchResult * pMatch) const
    {
        if (pMatch->chromosomeId_ != this->chromosomeId_)
            return false;
        if (pMatch->forward_ != this->forward_)
            return false;

        return !((this->opEndIndex_ <= pMatch->opStartIndex_) ||
                (this->opStartIndex_ >= pMatch->opEndIndex_));
    }

    static bool compareScore(const MatchResult * mr1, const MatchResult * mr2)
    {
        return (mr1->score_ < mr2->score_);
    }

    ////////////////////////////////////////////////////////////////////
    // Attributes

    // Contig information
    string contigId_;
    int contigSize_; // size of the contig (bp)

    // Alignment location & orientation
    string chromosomeId_;
    int opStartIndex_; // index of first aligned fragment from optical map
    int opEndIndex_; // index of last aligned fragment from optical map (inclusive)
    int opStartBp_; //bp position of the start of the match, relative to the beginning of the entire optical map
    int opEndBp_; // bp position of the end of the match, relative to the end of the entire optical map
    bool forward_; // contig oriented forward or backward
    vector<MatchedFrag> matchedFragList_;

    // Alignment statistics
    double score_;
    int totalMisses_; // number of missed restriction sites
    int totalHits_; // number of matched restriction sites
    double totalMissRate_; // fraction of contig sites that are missed in the contig + the aligned portion of the optical map
    int contigMisses_;
    int contigHits_; // number of matched contig sites
    int numAlignedInnerBlocks_; // number of aligned blocks
    double contigMissRate_;
    int contigUnalignedFrags_;
    int contigUnalignedBases_;
    int contigTotalAlignedBases_;
    int contigInnerAlignedBases_;
    int opticalMisses_;
    int opticalHits_; // number of matched optical sites
    double opticalMissRate_;
    int opticalTotalAlignedBases_;
    int opticalInnerAlignedBases_;
    double alignedLengthRatio_; // ratio of aligned contig & optical lengths
    double contigUnalignedBaseRatio_;
    double chi2_; // chi-squared score (sizing error)
    double pval_; // p-value from permutation test

    // Alignment Descriptive strings
    string contigMatchString_;
    string contigAlignedIndexString_;
    string contigLostIndexString_;
    string opticalMatchString_;
    string opticalAlignedIndexString_;


    private:
    bool builtAlignmentAttributes_;

};

ostream& operator<<(ostream& os, const MatchResult& mr);

#endif
