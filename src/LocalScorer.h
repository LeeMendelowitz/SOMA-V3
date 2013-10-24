#ifndef LOCAL_SCORER_H
#define LOCAL_SCORER_H


// In this scoring scheme, the scores are negative (representing penalties).
// The more positive (less negative) the score is, the better.

class LocalScorer
{

    public:
    LocalScorer(const AlignmentParams& ap) :
        ap_(ap)
        {}
    ~LocalScorer() {}

    // Local scorer does not allow gaps.
    inline Score scoreGap(int gapSize) { return Score(-Constants::INF, -Constants::INF, -Constants::INF); }

    inline Score scoreAlignment(const std::vector<FragData>::const_iterator cB,
                         const std::vector<FragData>::const_iterator cE,
                         const std::vector<FragData>::const_iterator oB,
                         const std::vector<FragData>::const_iterator oE,
                         bool boundaryFrag);
    Score scoreMatchedChunk(MatchedChunk& chunk);
    void scoreMatchResult(MatchResult * pResult);

    private:
    inline Score scoringFunction(int nContigSites, int nOpticalSites,
                           int contigLength, int opticalLength,
                           bool boundaryFrag);
    AlignmentParams ap_;
};


inline Score LocalScorer::scoreAlignment(const std::vector<FragData>::const_iterator cB,
                     const std::vector<FragData>::const_iterator cE,
                     const std::vector<FragData>::const_iterator oB,
                     const std::vector<FragData>::const_iterator oE,
                     bool boundaryFrag)
{
    int nContigSites = cE - cB - 1;
    int nOpticalSites = oE - oB - 1;
    int opticalLength = 0;
    int contigLength = 0;

    for(std::vector<FragData>::const_iterator iter = cB;
        iter != cE; 
        iter++)
        contigLength += iter->size_;

    for(std::vector<FragData>::const_iterator iter = oB;
        iter != oE;
        iter++)
        opticalLength += iter->size_;

    return scoringFunction(nContigSites, nOpticalSites, contigLength, opticalLength, boundaryFrag);
}

// Scoring functions for local alignment
// The greater (i.e. more positive) the score, the better the alignment
// nContigSites: number of unaligned contig sites
// nOpticalSites: number of unaligned optical sites
// contigLength: length of the contig alignment chunk
// opticalLength: length of the optical alignment chunk
// Returns a positive number for good alignments, a negative
// number for bad alignments

// The sizing score is given by a quadratic:
// f(delta) = (H/T^2)*(delta^2 - T^2)
// where delta = (opticalLength - contigLength)/VAR[opticalLength]
// Assuming that the observed optical fragment corresponds to the contig fragment,
// VAR[opticalLength] depends on the length of the contig chunk.
// Here, we use the model: VAR[opticalLength] = sigma^2*contigLength

inline Score LocalScorer::scoringFunction(int nContigSites, int nOpticalSites,
                       int contigLength, int opticalLength,
                       bool boundaryFrag)
{
    // Model assumes that Optical Frag ~ N(C, SIGMA^2*C) where C is contig frag size 
    double delta = (opticalLength - contigLength);
    double var = contigLength*ap_.sigma2;
    double chi2 = delta*delta/var;

    // Check if this alignment block exceeds the maximum allowed sizing error
    if (chi2 > ap_.chi2Max)
        chi2 = -Constants::INF;

    double sizeScore = -ap_.A*(chi2 - ap_.T2); // This is positive for small delta
    double contigScore = -nContigSites*ap_.C_r_contig;
    double opticalScore = -nOpticalSites*ap_.C_r_optical;
    return Score(contigScore, opticalScore, sizeScore);
}

Score LocalScorer::scoreMatchedChunk(MatchedChunk& chunk)
{
    Score score;

    if (chunk.isContigGap())
    {
        score = scoreGap(chunk.getContigMatchLengthBp());
    }
    else
    {
        score = scoreAlignment(chunk.getContigFragB(), chunk.getContigFragE(), chunk.getOpticalFragB(),
                               chunk.getOpticalFragE(), chunk.isBoundaryChunk());
    }

    chunk.setScore(score);
    return score;
}

// Score the matched chunks in a MatchResult
void LocalScorer::scoreMatchResult(MatchResult * pResult)
{
    std::vector<MatchedChunk>& chunkList = pResult->matchedChunkList_;
    std::vector<MatchedChunk>::iterator iter = chunkList.begin();
    std::vector<MatchedChunk>::iterator E = chunkList.end();
    for(; iter != E; iter++)
        scoreMatchedChunk(*iter);
}

#endif
