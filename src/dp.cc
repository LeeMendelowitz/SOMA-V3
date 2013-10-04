#include <sstream>
#include <cassert>
#include <cstdlib>
#include <set>
#include "itoa.h"
#include <algorithm>
#include <numeric>

#include "ScoreMatrix.h"
//#include "gamma-prob.c"
#include "globals.h"
#include "utils.h"
#include "match.h"
#include "dp.h"
#include "MatchResult.h"
#include "debugUtils.h"
#include "scoringFunctions.h"
#include "StandardMatchMaker.h"
#include "LocalMatchMaker.h"
#include "LocalScorer.h"
#include "GlobalScorer.h"

#define DEBUG_MATRIX 0
#define DEBUG_MATCH_LOCAL 0
#define DEBUG_MATCH_GLOBAL 0

using namespace std;

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


ScoreMatrix_t * createScoreMatrix1(const vector<FragData>& contigFrags, const vector<FragData>& opticalFrags,
    const int numFragsInChromosome, const AlignmentParams& alignParams);
ScoreMatrix_t * createScoreMatrix2(const vector<FragData>& contigFrags, const vector<FragData>& opticalFrags,
    const int numFragsInChromosome, const AlignmentParams& alignParams);

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
        pCur -> score_ = pPrev -> score_ + gapPenalty(contigFrags[i-1].size_, alignParams);
    }

    // Initialze rest of table
    d = ScoreElement_t(-Constants::INF,-1,-1);
    for (int i = 1; i < m; i++)
        for (int j = 1; j < n; j++)
            pScoreMatrix->d_[i*n+j] = d;

    // Dynamic programming
    int i1, j1, nSitesContig, nSitesOptical, cFrag;
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
                    boundaryFrag = ((k==0) || (i==m-1)) && opt::useBoundaries;
                    pPrev = &pScoreMatrix->d_[k*n+l];
                    if (pPrev->score_ > -Constants::INF)
                    {
                        nSitesContig = i-k-1;
                        nSitesOptical = j-l-1; 
                        newScore = pPrev->score_ +  scoringFunction(nSitesContig, nSitesOptical, cFragLength, oFragLength, boundaryFrag, alignParams);
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
                newScore = pPrev->score_ + gapPenalty(cFrag, alignParams);
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
        pCur -> score_ = pPrev -> score_ + gapPenalty(contigFrags[i-1].size_, alignParams);
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
                    boundaryFrag = ((k==0) || (i==m-1)) && opt::useBoundaries;
                    pPrev = &pScoreMatrix->d_[k*n+l];
                    if (pPrev->score_ > -Constants::INF)
                    {
                        newScore = pPrev->score_ + scoringFunction2(cbB, cbE, obB, obE, boundaryFrag, alignParams);
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
                newScore = pPrev->score_ + gapPenalty(cFragSize, alignParams);
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
ScoreMatrix_t * createLocalScoreMatrix(const vector<FragData>& contigFrags, const vector<FragData>& opticalFrags,
    const int numFragsInChromosome, const AlignmentParams& alignParams)
{
    const int m = contigFrags.size()+1; // number of rows
    const int n = opticalFrags.size()+1; // number of columns
    ScoreMatrix_t * pScoreMatrix = new ScoreMatrix_t(m ,n);

    // Initialize Matrix
    ScoreElement_t d(0.0, -1, -1);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            pScoreMatrix->d_[i*n+j] = d;

    // Dynamic programming
    int i1, j1, nSitesContig, nSitesOptical;
    int cFragLength, oFragLength;
    double newScore;
    bool boundaryFrag; // True if the fragment is not bounded by two restriction sites
    ScoreElement_t * pCur,  * pPrev;
    for (int i=1; i < m; i++)
    {
        for (int j=1; j < n; j++)
        {
            // Do not start an alignment in the second half of a circular chromosome
            if (i==1 && j > numFragsInChromosome)
                continue;
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
                    boundaryFrag = ((k==0) || (i==m-1)) && opt::useBoundaries;
                    if (boundaryFrag) continue;
                    pPrev = &pScoreMatrix->d_[k*n+l];
                    nSitesContig = i-k-1;
                    nSitesOptical = j-l-1; 
                    newScore = pPrev->score_ +  localScoringFunction(nSitesContig, nSitesOptical, cFragLength, oFragLength, boundaryFrag, alignParams);
                    if (newScore > pCur->score_)
                    {
                        pCur->score_ = newScore;
                        pCur->pi_ = k;
                        pCur->pj_ = l;
                    }
                } // end k loop
            } // end l loop
        } // end j loop
    } // end i loop

    return pScoreMatrix;
} // end create score matrix


// Match the fragments from a single contig to the entire optical map
// return the best match as a new MatchResult
// if pAllMatches is not NULL, add all matches (inlcuding the best one) to pAllMatches vector.
MatchResult *  match(const ContigMapData * pContigMap, const OpticalMapData * pOpticalMap,
                     vector<MatchResult *> * pAllMatches, const AlignmentParams& alignParams)
{
    MatchResult * bestMatch = 0;
    bool contigIsForward = pContigMap->isForward();
    const vector<FragData>& opticalFrags = pOpticalMap->getFrags();
    const vector<FragData>& contigFrags = pContigMap->getFrags();
    //const int m = contigFrags.size() + 1; // num rows
    //const int n = opticalFrags.size() + 1; // num cols
    //const int lr = m-1;
    //const int lro = lr*n;

    // Create the scoreMatrix
    ScoreMatrix_t * pScoreMatrix = createScoreMatrix2(contigFrags, opticalFrags, pOpticalMap->getNumFrags(), alignParams);

    // Build matches from the score matrix
    GlobalScorer scorer(alignParams);
    StandardMatchMaker matchMaker(&scorer, opt::maxMatchesPerContig, opt::minContigHits,
                                  opt::minLengthRatio, opt::maxMissRateContig, opt::avgChi2Threshold);
    MatchResultPtrVec matches;
    matchMaker.makeMatches(pScoreMatrix, matches, pOpticalMap, pContigMap, contigIsForward);

    #if DEBUG_MATCH_GLOBAL > 0
    std::cout << "Found " << matches.size()
              << " matches for contig " << pContigMap->getId()
              << " to " << pOpticalMap->getId()
              << " in " << (contigIsForward ? " forward " : " reverse ") 
              << " direction." << std::endl;
    #endif

    // Set the matches to be returned
    if (!matches.empty())
    {
        // The first match is the best match
        bestMatch = matches.front();

        if (pAllMatches)
        {
            pAllMatches->insert(pAllMatches->end(), matches.begin(), matches.end());
        }
        else
        {
            // Delete all other matches other than the best match
            const MatchResultPtrVec::iterator E = matches.end();
            for(MatchResultPtrVec::iterator iter = matches.begin() + 1;
                iter != E;
                iter++)
            {
                delete *iter;
            }
            matches.resize(1);
        }
    }

    #if DEBUG_MATCH_GLOBAL > 0
    std::cout << "pAllMatches->size(): " << pAllMatches->size() << std::endl;
    #endif

    delete pScoreMatrix;
    return bestMatch;
}

// Match the fragments from a single contig to the entire optical map
// return the best match as a new MatchResult
// if pAllMatches is not NULL, add all matches (inlcuding the best one) to pAllMatches vector.
// This is a copy of the match() function, but to test the 
MatchResult *  match2(const ContigMapData * pContigMap, const OpticalMapData * pOpticalMap,
                     vector<MatchResult *> * pAllMatches, const AlignmentParams& alignParams)
{
    MatchResult * bestMatch = 0;
    bool contigIsForward = pContigMap->isForward();
    const vector<FragData>& opticalFrags = pOpticalMap->getFrags();
    const vector<FragData>& contigFrags = pContigMap->getFrags();
    //const int m = contigFrags.size() + 1; // num rows
    //const int n = opticalFrags.size() + 1; // num cols
    //const int lr = m-1;
    //const int lro = lr*n;

    // Create the scoreMatrix
    ScoreMatrix_t * pScoreMatrix = createScoreMatrix1(contigFrags, opticalFrags, pOpticalMap->getNumFrags(), alignParams);

    // Build matches from the score matrix
    GlobalScorer scorer(alignParams);
    StandardMatchMaker matchMaker(&scorer, opt::maxMatchesPerContig, opt::minContigHits,
                                  opt::minLengthRatio, opt::maxMissRateContig, opt::avgChi2Threshold);
    MatchResultPtrVec matches;
    matchMaker.makeMatches(pScoreMatrix, matches, pOpticalMap, pContigMap, contigIsForward);

    #if DEBUG_MATCH_GLOBAL > 0
    std::cout << "Found " << matches.size()
              << " matches for contig " << pContigMap->getId()
              << " to " << pOpticalMap->getId()
              << " in " << (contigIsForward ? " forward " : " reverse ") 
              << " direction." << std::endl;
    #endif

    // Set the matches to be returned
    if (!matches.empty())
    {
        // The first match is the best match
        bestMatch = matches.front();

        if (pAllMatches)
        {
            pAllMatches->insert(pAllMatches->end(), matches.begin(), matches.end());
        }
        else
        {
            // Delete all other matches other than the best match
            const MatchResultPtrVec::iterator E = matches.end();
            for(MatchResultPtrVec::iterator iter = matches.begin() + 1;
                iter != E;
                iter++)
            {
                delete *iter;
            }
            matches.resize(1);
        }
    }

    #if DEBUG_MATCH_GLOBAL > 0
    std::cout << "pAllMatches->size(): " << pAllMatches->size() << std::endl;
    #endif

    delete pScoreMatrix;
    return bestMatch;
}

// Match the fragments from a single contig to the entire optical map
// return the best match as a new MatchResult
// if pAllMatches is not NULL, add all matches (inlcuding the best one) to pAllMatches vector.
MatchResult *  matchLocal(const ContigMapData * pContigMap, const OpticalMapData * pOpticalMap,
                     vector<MatchResult *> * pAllMatches, const AlignmentParams& alignParams)
{
    MatchResult * bestMatch = 0;
    bool contigIsForward = pContigMap->isForward();
    const vector<FragData>& opticalFrags = pOpticalMap->getFrags();
    const vector<FragData>& contigFrags = pContigMap->getFrags();
    const int m = contigFrags.size() + 1; // num rows
    const int n = opticalFrags.size() + 1; // num cols

    // Create the scoreMatrix
    ScoreMatrix_t * pScoreMatrix = createLocalScoreMatrix(contigFrags, opticalFrags, pOpticalMap->getNumFrags(), alignParams);

    #if DEBUG_MATRIX > 0 
    string debugFileName;
    debugFileName = pContigMap->getId() + "_" + pOpticalMap->getId();
    if (contigIsForward)
        debugFileName += ".forward.txt";
    else
        debugFileName += ".reverse.txt";
    writeMatrixToFile(pScoreMatrix, debugFileName);
    #endif

    assert (pScoreMatrix->m_ == m);
    assert (pScoreMatrix->n_ == n);

    // Build matches from the score matrix
    LocalScorer localScorer(alignParams);
    LocalMatchMaker matchMaker(&localScorer, opt::maxMatchesPerContig, opt::minContigHits,
                               opt::minLengthRatio, opt::maxMissRateContig,
                               opt::avgChi2Threshold);
    MatchResultPtrVec matches;
    matchMaker.makeMatches(pScoreMatrix, matches, pOpticalMap, pContigMap, contigIsForward);

    #if DEBUG_MATCH_LOCAL > 0
    std::cout << "Found " << matches.size()
              << " local alignment for contig " << pContigMap->getId()
              << " to " << pOpticalMap->getId()
              << " in " << (contigIsForward ? " forward " : " reverse ") 
              << " direction." << std::endl;
    #endif

    // Set the matches to be returned
    if (!matches.empty())
    {
        // The first match is the best match
        bestMatch = matches.front();

        if (pAllMatches)
        {
            pAllMatches->insert(pAllMatches->end(), matches.begin(), matches.end());
        }
        else
        {
            // Delete all other matches other than the best match
            const MatchResultPtrVec::iterator E = matches.end();
            for(MatchResultPtrVec::iterator iter = matches.begin() + 1;
                iter != E;
                iter++)
            {
                delete *iter;
            }
            matches.resize(1);
        }
    }

    #if DEBUG_MATCH_GLOBAL > 0
    std::cout << "pAllMatches->size(): " << pAllMatches->size() << std::endl;
    #endif

    delete pScoreMatrix;
    return bestMatch;
}

// Match the fragments from a single contig to the entire optical map. Return the best score.
// Note: This is a light-weight alignment, since it avoids constructing MatchResults.
// Note: The permutation test WILL NOTE USE filterFunction to discard bad matches. Therefore,
// if the match penalties or scoring functions are not set reasonably, it may be possible that the
// best MatchResult (after filtering) will have a low score, which will easily be beat by the permutation alignments,
// which are not filtered.
double matchPermutationTest(const ContigMapData * pContigMap, const OpticalMapData * pOpticalMap,
                     const AlignmentParams& alignParams)
{
    double best_score = -Constants::INF;
    Index_t best_index = Index_t(-1,-1);
    const vector<FragData>& opticalFrags = pOpticalMap->getFrags();
    const vector<FragData>& contigFrags = pContigMap->getFrags();
    const int m = contigFrags.size() + 1; // num rows
    const int n = opticalFrags.size() + 1; // num cols
    const int lr = m-1;
    const int lro = lr*n;
    const ScoreElement_t * pE;
    bool foundMatch = false;
    bool contigHit = false;

    // Create the scoreMatrix
    ScoreMatrix_t * pScoreMatrix = createScoreMatrix2(contigFrags, opticalFrags, pOpticalMap->getNumFrags(), alignParams);

    assert (pScoreMatrix->m_ == m);
    assert (pScoreMatrix->n_ == n);

    // In first pass over last row of the pScoreMatrix, find the best match
    if (opt::oneToOneMatch)
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
