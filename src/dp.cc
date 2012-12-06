#include <sstream>
#include <assert.h>
#include <stdlib.h>
#include <set>
#include "itoa.h"
#include <algorithm>
#include <numeric>

#include "ScoreMatrix.h"
#include "gamma-prob.c"
#include "globals.h"
#include "utils.h"
#include "match.h"
#include "dp.h"
#include "MatchResult.h"
#include "debugUtils.h"

using namespace std;

// Penalty for "losing" a small contig fragment.
// (gapped alignment)
double gapPenalty(int fragSize)
{
    if (fragSize < Options::smallFrag)
        return -fragSize * Options::smallFragSlope;
    else
        return -Constants::INF;
}

// Penalty for missing a contig restriction site
double contigMissedSitePenalty(int dToClosestSite)
{
    if (dToClosestSite > Options::smallFrag)
        return -Options::C_r_contig;
    else
    {   
        // Add a small amount to the penalty to break ties with the gapPenalty
        // in boundary fragments where there is no sizing error
        return (gapPenalty(dToClosestSite) + 1.0E-6);
    }
}


// nContigSites: number of unaligned contig sites
// nOpticalSites: number of unaligned optical sites
double scoringFunction(int nContigSites, int nOpticalSites,
                       int contigLength, int opticalLength,
                       bool boundaryFrag)
{
    // Model assumes that Optical Frag ~ N(C, SIGMA^2*C) where C is contig frag size 

    double chi2;
    
    if (boundaryFrag && contigLength < opticalLength)
        chi2 = 0.0; // Do not penalize boundary fragments for being too small
    else
    {
        double var = contigLength*Constants::SIGMA2;
        double dl = (contigLength-opticalLength);
        chi2 = dl*dl/var;
        assert (chi2 >= -1E-12);
        // Only allow Options::sdMax standard deviations in sizing error
        if (chi2 > Options::sdMax*Options::sdMax)
            return -Constants::INF;
    }
    return -(nContigSites*Options::C_r_contig +
           nOpticalSites*Options::C_r_optical +
           chi2);
}


// Scoring function which accounts for the distance to the closest
// aligned site for each unaligned interior site.
double scoringFunction2( const vector<FragData>::const_iterator& cB,
                         const vector<FragData>::const_iterator& cE,
                         const vector<FragData>::const_iterator& oB,
                         const vector<FragData>::const_iterator& oE,
                         bool boundaryFrag )
{
    vector<FragData>::const_iterator ci, oi;

    // Calculate total contig and optical block length
    int contigLength = 0;
    for (ci = cB; ci != cE; ci++)
        contigLength += ci->size_;

    int opticalLength = 0;
    for (oi = oB; oi != oE; oi++)
        opticalLength += oi->size_;

    // Compute sizing error
    double chi2 = 0.0;
    if (boundaryFrag && contigLength < opticalLength)
        chi2 = 0.0; // Do not penalize boundary fragments for being too small
    else
    {
        double var = contigLength*Constants::SIGMA2;
        double dl = (contigLength-opticalLength);
        chi2 = dl*dl/var;
        assert (chi2 >= -1E-12);
        // Only allow Options::sdMax standard deviations in sizing error
        if (chi2 > Options::sdMax*Options::sdMax)
            return -Constants::INF;
    }
    chi2 = -1.0*chi2;

    // Compute miss penalty
    double contigMissScore = 0.0;
    int numOpticalFrags = oE - oB;
    double opticalMissScore = -Options::C_r_optical*(numOpticalFrags-1);
    int pos = 0;
    for(ci = cB; ci != cE-1; ci++)
    {
        pos += ci->size_;
        // The miss score is a function of the distance to the first missed site
        contigMissScore += contigMissedSitePenalty(min(pos, contigLength-pos));
    }

    return (contigMissScore + opticalMissScore + chi2);
}




class IndexVectorCmp_t
{
    public:
    ScoreMatrix_t * pMat_;
    bool operator() (const Index_t& ind1, const Index_t& ind2)
    {
        ScoreElement_t * pE1 = &pMat_->d_[(pMat_->n_ - 1)*ind1.first + ind1.second];
        ScoreElement_t * pE2 = &pMat_->d_[(pMat_->n_ - 1)*ind2.first + ind2.second];
        return (pE1->score_ < pE2->score_);
    }
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Populate a ScoreMatrix_t
// Score all possible alignments
// Use the scoring function which take as input fragment lengths
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
ScoreMatrix_t * createScoreMatrix1(const vector<FragData>& contigFrags, const vector<FragData>& opticalFrags,
    const int numFragsInChromosome, const AlignmentParams& alignParams)
{
    const int m = contigFrags.size()+1; // number of rows
    const int n = opticalFrags.size()+1; // number of columns
    ScoreMatrix_t * pScoreMatrix = new ScoreMatrix_t(m ,n);

    // Initialize first row
    ScoreElement_t d(0.0, -1, -1);
    for (int j=0; j < n; j++)
    {
        pScoreMatrix->d_[j] = d;
    }

    // Initialize first column. This corresponds alignments which start with "lost" fragments from contig
    ScoreElement_t * pCur,  * pPrev;
    for (int i=1; i < m; i++)
    {
        pPrev = &(pScoreMatrix->d_[(i-1)*n]);
        pCur = &(pScoreMatrix->d_[i*n]);
        pCur -> pi_ = i-1;
        pCur -> pj_ = 0;
        // Calc score for gapped alignment
        pCur -> score_ = pPrev -> score_ + gapPenalty(contigFrags[i-1].size_);
    }

    // Initialze rest of table
    d = ScoreElement_t(-Constants::INF,-1,-1);
    for (int i = 1; i < m; i++)
        for (int j = 1; j < n; j++)
            pScoreMatrix->d_[i*n+j] = d;

    // Dynamic programming
    int i1, j1, nSitesContig, nSitesOptical, cFrag, oFrag;
    int cFragLength, oFragLength;
    double newScore;
    bool boundaryFrag; // True if the fragment is not bounded by two restriction sites
    for (int i=1; i < m; i++)
    {
        cFrag = contigFrags[i-1].size_;
        for (int j=1; j < n; j++)
        {
            // Do not start an alignment in the second half of a circular chromosome
            if (i==1 && j > numFragsInChromosome)
                continue;
            oFrag = opticalFrags[j-1].size_;
            pCur = &pScoreMatrix->d_[i*n+j];

            // Determine eligible extensions
            i1 = max(0, i-alignParams.delta-1);
            j1 = max(0, j-alignParams.delta-1); 
            cFragLength = 0;
            for (int k = i-1; k >= i1; k--) // Loop over contig frags for alignment block
            {
                oFragLength = 0;
                cFragLength += contigFrags[k].size_;
                for (int l = j-1; l >= j1; l--) // Loop over optical frags for alignment block
                {
                    oFragLength += opticalFrags[l].size_;
                    boundaryFrag = (k==0) || (i==m-1);
                    pPrev = &pScoreMatrix->d_[k*n+l];
                    if (pPrev->score_ > -Constants::INF)
                    {
                        nSitesContig = i-k-1;
                        nSitesOptical = j-l-1; 
                        newScore = pPrev->score_ +  scoringFunction(nSitesContig, nSitesOptical, cFragLength, oFragLength, boundaryFrag);
                        if (newScore > pCur->score_)
                        {
                            pCur->score_ = newScore;
                            pCur->pi_ = k;
                            pCur->pj_ = l;
                        }
                    }
                } // end k loop
            } // end l loop

            // Try a gap open/extension
            pPrev = &pScoreMatrix->d_[(i-1)*n+j];
            if (pPrev->score_ > -Constants::INF)
            {
                newScore = pPrev->score_ + gapPenalty(cFrag);
                if (newScore > pCur->score_)
                {
                    pCur->score_ = newScore;
                    pCur->pi_ = i-1;
                    pCur->pj_ = j;
                }
            }
        } // end j loop
    } // end i loop

    return pScoreMatrix;
} // end create score matrix


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Populate a ScoreMatrix_t
// Score all possible alignments
// Use the scoring function which take as input the restriction fragment pattern for a block
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
ScoreMatrix_t * createScoreMatrix2(const vector<FragData>& contigFrags, const vector<FragData>& opticalFrags,
    const int numFragsInChromosome, const AlignmentParams& alignParams)
{
    const int m = contigFrags.size()+1; // number of rows
    const int n = opticalFrags.size()+1; // number of columns
    ScoreMatrix_t * pScoreMatrix = new ScoreMatrix_t(m ,n);

    const vector<FragData>::const_iterator cB = contigFrags.begin();
    const vector<FragData>::const_iterator oB = opticalFrags.begin();

    // Initialize first row
    ScoreElement_t d(0.0, -1, -1);
    for (int j=0; j < n; j++)
    {
        pScoreMatrix->d_[j] = d;
    }

    // Initialize first column. This corresponds alignments which start with "lost" fragments from contig
    ScoreElement_t * pCur,  * pPrev;
    for (int i=1; i < m; i++)
    {
        pPrev = &(pScoreMatrix->d_[(i-1)*n]);
        pCur = &(pScoreMatrix->d_[i*n]);
        pCur -> pi_ = i-1;
        pCur -> pj_ = 0;
        // Calc score for gapped alignment
        pCur -> score_ = pPrev -> score_ + gapPenalty(contigFrags[i-1].size_);
    }

    // Initialze rest of table
    d = ScoreElement_t(-Constants::INF,-1,-1);
    for (int i = 1; i < m; i++)
        for (int j = 1; j < n; j++)
            pScoreMatrix->d_[i*n+j] = d;

    // Dynamic programming
    int i1, j1;
    double newScore;
    vector<FragData>::const_iterator cbB, cbE, obB, obE; // Beginning and ending iterators for contig and optical block
    bool boundaryFrag; // True if the block is not bounded by two restriction sites
    int cFragSize;
    for (int i=1; i < m; i++)
    {
        cbE = cB + i; // one past the last contig frag in the alignment block
        cFragSize = (cbE-1)->size_;
        for (int j=1; j < n; j++)
        {
            // Do not start an alignment in the second half of a circular chromosome
            if (i==1 && j > numFragsInChromosome)
                continue;

            obE = oB + j; // one past the last optical frag in the alignment block
            pCur = &pScoreMatrix->d_[i*n+j];

            // Determine eligible extensions
            i1 = max(0, i-alignParams.delta);
            j1 = max(0, j-alignParams.delta); 
            for (int k = i-1; k >= i1; k--) // Loop over contig frags for alignment block
            {
                cbB = cB + k;
                for (int l = j-1; l >= j1; l--) // Loop over optical frags for alignment block
                {
                    obB = oB + l;
                    boundaryFrag = (k==0) || (i==m-1);
                    pPrev = &pScoreMatrix->d_[k*n+l];
                    if (pPrev->score_ > -Constants::INF)
                    {
                        newScore = pPrev->score_ + scoringFunction2(cbB, cbE, obB, obE, boundaryFrag);
                        if (newScore > pCur->score_)
                        {
                            pCur->score_ = newScore;
                            pCur->pi_ = k;
                            pCur->pj_ = l;
                        }
                    }
                } // end k loop
            } // end l loop

            // Try a gap open/extension
            pPrev = &pScoreMatrix->d_[(i-1)*n+j];
            if (pPrev->score_ > -Constants::INF)
            {
                newScore = pPrev->score_ + gapPenalty(cFragSize);
                if (newScore > pCur->score_)
                {
                    pCur->score_ = newScore;
                    pCur->pi_ = i-1;
                    pCur->pj_ = j;
                }
            }
        } // end j loop
    } // end i loop

    return pScoreMatrix;
} // end create score matrix

// Match the fragments from a single contig to the entire optical map
// return the best match as a new MatchResult
// if pAllMatches is not NULL, add all matches (inlcuding the best one) to pAllMatches vector.
MatchResult *  match(const ContigMapData * pContigMapData, const OpticalMapData * pOpticalMapData,
                     vector<MatchResult *> * pAllMatches, bool forward, const AlignmentParams& alignParams)
{
    MatchResult * bestMatch = 0;
    double best_score = -Constants::INF;
    Index_t best_index = Index_t(-1,-1);
    vector<FragData> opticalFrags = pOpticalMapData->frags_;
    const vector<FragData>& contigFrags = forward ? pContigMapData->frags_: pContigMapData->reverseFrags_;
    const int m = contigFrags.size() + 1; // num rows
    const int n = opticalFrags.size() + 1; // num cols
    const int lr = m-1;
    const int lro = lr*n;
    const ScoreElement_t * pE;
    bool foundMatch = false;

    // Create the scoreMatrix
    ScoreMatrix_t * pScoreMatrix = createScoreMatrix2(contigFrags, opticalFrags, pOpticalMapData->numFrags_, alignParams);

    /*
    string debugFileName;
    if (forward)
        debugFileName = "debugMatrix_forward.txt";
    else
        debugFileName = "debugMatrix_reverse.txt";
    writeMatrixToFile(pScoreMatrix, debugFileName);
    */

    assert (pScoreMatrix->m_ == m);
    assert (pScoreMatrix->n_ == n);

    // In first pass over last row of the pScoreMatrix, find the best match
    if (Options::oneToOneMatch)
    {
        // We are matching the entire contig to the entire optical map
        pE = &pScoreMatrix->d_[lro + n-1];
        if (pE->score_ > -Constants::INF)
        {
            best_index = Index_t(lr,n-1); // Last element in scoreMatrix
            foundMatch = true;
        }
    }
    else
    {
        for(int j=1; j < n; j++)
        {
            pE = &pScoreMatrix->d_[lro + j];
            if (pE->score_ > best_score)
            {
                best_score = pE->score_;
                best_index = Index_t(lr, j);
                foundMatch = true;
            }
        }
    }

    // Create a MatchResult for the best scoring alignment, and make sure the alignment is valid.
    if (foundMatch)
    {
        bestMatch = new MatchResult(best_index, pScoreMatrix, pContigMapData, pOpticalMapData, forward);
        assert(bestMatch!=0);
        bestMatch->buildAlignmentAttributes();
        if (bestMatch->contigHits_ <= 0)
        {
            // The best scoring match is not valid, since all contigs fragments
            // are in a gap alignment!
            foundMatch = false;
            delete bestMatch;
            bestMatch = 0;
        } 
        else
        {
            if (pAllMatches)
                pAllMatches->push_back(bestMatch);
        }
    }

    // If necessary, find other alignments
    if(!Options::oneToOneMatch && pAllMatches && foundMatch)
    {
        assert (bestMatch != 0);
        assert (bestMatch->contigHits_ > 0);

        vector<Index_t> inds; //indices of candidate alignments other than the best match
        inds.reserve(n);
        for(int j = 1; j < n; j++)
        {
            if(j == best_index.second)
                continue;
            pE = &pScoreMatrix->d_[lro + j];
            if(pE->score_ > -Constants::INF)
            {
                inds.push_back(Index_t(lr, j));
            }
        }

        // Add all other Match Results
        pAllMatches->reserve(pAllMatches->size() + inds.size());
        vector<Index_t>::iterator it,ite;
        it = inds.begin();
        ite = inds.end();
        for (; it != ite; it++)
        {
            MatchResult * other = new MatchResult(*it, pScoreMatrix, pContigMapData, pOpticalMapData, forward);
            assert(other!=0);
            other->buildAlignmentAttributes();
            if (other->contigHits_ > 0)
                pAllMatches->push_back(other);
            else
                delete other; // This alignment is not valid
        }
    }

    delete pScoreMatrix;
    return bestMatch;
}

// Match the fragments from a single contig to the entire optical map. Return the best score.
// Note: This is a light-weight alignment, since it avoids constructing MatchResults.
// Note: The permutation test WILL NOTE USE filterFunction to discard bad matches. Therefore,
// if the match penalties or scoring functions are not set reasonably, it may be possible that the
// best MatchResult (after filtering) will have a low score, which will easily be beat by the permutation alignments,
// which are not filtered.
double matchPermutationTest(const ContigMapData * pContigMapData, const OpticalMapData * pOpticalMapData,
                     bool forward, const AlignmentParams& alignParams)
{
    MatchResult * bestMatch = 0;
    double best_score = -Constants::INF;
    Index_t best_index = Index_t(-1,-1);
    vector<FragData> opticalFrags = pOpticalMapData->frags_;
    const vector<FragData>& contigFrags = forward ? pContigMapData->frags_: pContigMapData->reverseFrags_;
    const int m = contigFrags.size() + 1; // num rows
    const int n = opticalFrags.size() + 1; // num cols
    const int lr = m-1;
    const int lro = lr*n;
    const ScoreElement_t * pE;
    bool foundMatch = false;
    bool contigHit = false;

    // Create the scoreMatrix
    ScoreMatrix_t * pScoreMatrix = createScoreMatrix2(contigFrags, opticalFrags, pOpticalMapData->numFrags_, alignParams);

    assert (pScoreMatrix->m_ == m);
    assert (pScoreMatrix->n_ == n);

    // In first pass over last row of the pScoreMatrix, find the best match
    if (Options::oneToOneMatch)
    {
        // We are matching the entire contig to the entire optical map
        pE = &pScoreMatrix->d_[lro + n-1];
        if (pE->score_ > -Constants::INF)
        {
            best_index = Index_t(lr,n-1); // Last element in scoreMatrix
            foundMatch = true;
        }
    }
    else
    {
        for(int j=1; j < n; j++)
        {
            pE = &pScoreMatrix->d_[lro + j];
            if (pE->score_ > best_score)
            {
                best_score = pE->score_;
                best_index = Index_t(lr, j);
                foundMatch = true;
            }
        }
    }

    // Create a MatchResult for the best scoring alignment, and make sure the alignment is valid.
    if (foundMatch)
    {
        // Make sure that the number of contig hits is greater than 1. Then return the score.
        // This avoids returning a score of 0 for an alignment where every fragment is in a 
        // gapped alignment.
        Index_t curInd, prevInd;
        ScoreElement_t * pPrev;
        prevInd = best_index;
        pPrev = &(pScoreMatrix->d_[n*prevInd.first + prevInd.second]);
        while(true)
        {
            curInd = Index_t(pPrev->pi_, pPrev->pj_); // index to next element in trail backwards
            if (curInd.first == 0)
                break;
            if (curInd.first > 0 && curInd.first < prevInd.first)
            {
                contigHit = true;
                break;
            }
            pPrev = &(pScoreMatrix->d_[n*curInd.first + curInd.second]);
            prevInd = curInd;
        }
    }

    if (!contigHit)
        best_score = -Constants::INF;

    delete pScoreMatrix;
    return best_score;
}
