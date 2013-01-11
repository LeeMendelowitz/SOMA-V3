#ifndef MATCHRESULT_H
#define MATCHRESULT_H

#include <cstdlib>
#include <numeric>

#include "ScoreMatrix.h"
#include "mapTypes.h"
#include "OpticalMapData.h"
#include "ContigMapData.h"
#include "Scorer.h"
#include "MatchedChunk.h"

using namespace std;

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

    MatchResult(const string& contigId, const string& opticalId, int contigSize, bool contigIsForward, double score)
    {
        reset();
        contigId_ = contigId;
        chromosomeId_ = opticalId;
        forward_ = contigIsForward;
        contigLength_ = contigSize;
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
    // Two MatchResults overlap if they align to the same strand of the same reference
    // map.
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
    int contigLength_; // total contig length (in bp)

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
    vector<Score> scoreList_; // scores for each chunk

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
    string scoreString_;
};

ostream& operator<<(ostream& os, const MatchResult& mr);

typedef pair<int,int> Index_t;

#endif
