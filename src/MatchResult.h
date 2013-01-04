#ifndef MATCHRESULT_H
#define MATCHRESULT_H

#include <cstdlib>
#include <numeric>

#include "ScoreMatrix.h"
#include "mapTypes.h"
#include "OpticalMapData.h"
#include "ContigMapData.h"

using namespace std;


// Represents an aligned fragment block ("chunk")
class MatchedChunk
{
    public:

        typedef typename vector<FragData>::const_iterator FragDataConstIter;

        // Constructor
        //MatchedChunk () {}; // Default constructor
        MatchedChunk(int opStart, int opEnd, int opStartBp, int opEndBp, const vector<FragData>& oData,
                    int cStart, int cEnd, int cStartBp, int cEndBp, const vector<FragData>& cData) :
                    opStart_(opStart), opEnd_(opEnd), opLength_(opEndBp - opStartBp),
                    cStart_(cStart), cEnd_(cEnd), cLength_(abs(cEndBp - cStartBp)),
                    opStartBp_(opStartBp), opEndBp_(opEndBp),
                    cStartBp_(cStartBp), cEndBp_(cEndBp), 
                    cMisses_(cEnd-cStart-1)
        {

            contigGap_ = (opStart_ == opEnd_);
            opFragB_ = oData.begin() + opStart_;
            opFragE_ = oData.begin() + opEnd_;
            cFragB_ = cData.begin() + cStart_;
            cFragE_ = cData.begin() + cEnd_;

            opMisses_ = contigGap_ ? 0 : opEnd - opStart - 1;

            assert ( (opFragB_ <= opFragE_) && (opFragB_ >= oData.begin()) && (opFragE_ <= oData.end()) );
            assert ( (cFragB_ <= cFragE_) && (cFragB_ >= cData.begin()) && (cFragE_ <= cData.end()) );
            assert ( ( cStart_ >= 0 ) && ((size_t) cStart_ <= cData.size()) );
            assert ( ( cEnd_ >= 0) && ((size_t) cEnd_ <= cData.size()) );
            assert ( cStart_ <= cEnd_ );

            // This is a boundary chunk if it includes the first contig fragment
            // or the last contig fragment.
            boundaryChunk_ =  cStart_ == 0 || ((size_t) cEnd_ == cData.size());
        }

        // Data
        bool contigGap_;
        bool boundaryChunk_; // True if this is either the first or last matched chunk and is not bounded by two matched restriction sites. 
        int opStart_; // optical start index (zero based, inclusive)
        int opEnd_; // optical end index (exclusive)
        int opLength_; // optical fragment length (bp)
        int cStart_; // contig start index (zero based, inclusive)
        int cEnd_; // contig end index (exclusive)
        int cLength_; // contig fragment length (bp)

        // Note: The optical map starting/ending positions are only approximate for the case
        // of boundary chunks.
        int opStartBp_; // approximate optical start bp (zero based, inclusive)
        int opEndBp_; // approximate optical end bp (exclusive).
        int cStartBp_; // contig fragment start, bp. (zero based, inclusive)
        int cEndBp_; // contig end, bp. (exclusive)

        int opMisses_; // number of unaligned optical sites within the block
        int cMisses_; // number of unaligned contig sites within the block

        FragDataConstIter cFragB_; // pointer to first contig fragment in chunk
        FragDataConstIter cFragE_; // pointer to one beyond the last contig fragment in chunk
        FragDataConstIter opFragB_; // pointer to first optical fragment in chunk
        FragDataConstIter opFragE_; // pointer to one beyond the last optical fragment in chunk
};

typedef std::vector<MatchedChunk> MatchedChunkVec;

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

    MatchResult(const string& contigId, const string& opticalId, bool contigIsForward, double score)
    {
        reset();
        contigId_ = contigId;
        chromosomeId_ = opticalId;
        forward_ = contigIsForward;
        score_ = score;
    };

    // Reset all attributes
    void reset();

    // Set matched chunks, which defines the alignment.
    void setMatchedChunks(const MatchedChunkVec& chunks);

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
    int contigSize_; // size of the matched portion of the contig (bp)

    // Alignment location & orientation
    string chromosomeId_;
    int opStartIndex_; // index of first aligned fragment from optical map
    int opEndIndex_; // index of last aligned fragment from optical map (inclusive)
    int opStartBp_; //bp position of the start of the match, relative to the beginning of the entire optical map
    int opEndBp_; // bp position of the end of the match, relative to the end of the entire optical map
    int cStartIndex_;
    int cEndIndex_; // (inclusive)
    int cStartBp_;
    int cEndBp_;
    bool forward_; // contig oriented forward or backward
    vector<MatchedChunk> matchedChunkList_;

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
};

ostream& operator<<(ostream& os, const MatchResult& mr);

typedef pair<int,int> Index_t;

#endif
