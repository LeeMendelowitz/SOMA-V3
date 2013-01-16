// match.cc

#include <vector>
#include <cstdlib>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <set>
#include <assert.h>
#include <omp.h>
#include <cstdlib>
#include <cassert>

#include "globals.h"
#include "parseArgs.h"
#include "utils.h"
#include "dp.h"
#include "MatchResult.h"
#include "match.h"
#include "xmlWriter.h"
#include "exception.h"

//#include <boost/math/distributions.hpp> 
//#include "gamma-prob.c"

#define DEBUG_MATCH 0
#define DEBUG_FILTER 0
#define DEBUG_READMAPS 1

using namespace std;
using namespace Globals; // from globals.h

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) > (b) ? (b) : (a))


float betai(float a, float b, float x);

typedef vector<MatchResult *> MatchResultPtrVec;
typedef map<ContigMapData *, MatchResultPtrVec> MatchResultMap;

ostream& operator<<(ostream& ss, const SiteData& siteData)
{
    return (ss << siteData.loc_);
}

// Sort Algorithm Comparisons
bool contigDataVec_comp(const ContigMapData * pContig1, const ContigMapData * pContig2)
    { return pContig1->getNumFrags() > pContig2->getNumFrags(); }
bool contigSiteDataVec_lt(const SiteData& sd1, const SiteData& sd2) {return sd1.loc_ < sd2.loc_;}
bool contigFragDataVec_lt(const FragData& fd1, const FragData& fd2) { return fd1.size_ < fd2.size_; }
bool resultScoreComp(const MatchResult * r1, const MatchResult * r2) { return r1->score_ > r2->score_; }


void matchContigToOpticalMaps(const ContigMapData * pContigMap, const vector<OpticalMapData *>& opticalMapList,
                              vector<MatchResult *> * const pResults, bool localAlignment = false);

void runPermutationTests(MatchResultMap * pMatchResultMap, const vector<OpticalMapData *> opticalMapList, int numTrials);
void runPermutationTests2(MatchResultMap * pMatchResultMap, const vector<OpticalMapData *> opticalMapList, const vector<ContigMapData *>& contigVec, int numTrials);

// Create a new vector replaced from orig where elements at positions
// are replaced by elements from replacements
void replace(const vector<FragData>& orig, const vector<int>& positions,
                    const vector<FragData>& replacements, vector<FragData>& replaced) {

    replaced = orig;
    vector<int>::const_iterator iter = positions.begin();
    const vector<int>::const_iterator end = positions.end();
    for(; iter != end; iter++)
        replaced[*iter] = replacements[*iter];
    return;
}


// Computer the P-value by seeing what fraction of those entries in scoreVec
// score higher than myScore
double computePValue(double myScore, const vector<double>& scoreVec)
{
    int tailCount = 0;
    int numSamples = scoreVec.size();
    assert(numSamples > 0);
    for(int i=0; i<numSamples; i++)
        if (scoreVec[i] > myScore) tailCount++;
    return ((double) tailCount)/numSamples;
}

void readMaps(vector<OpticalMapData *>& opMapVec, vector<ContigMapData *>& contigVec)
{
    size_t numOpticalMaps = opt::opticalMapList.size();
    for(size_t i=0; i<numOpticalMaps; i++)
        readMaps(opt::opticalMapList[i], opMapVec);

    // Read Silico file with contig map data    
    readMaps(opt::silicoMap, contigVec);


    // Sort contigs in contigVec by the number of sites
    sort(contigVec.begin(), contigVec.end(), contigDataVec_comp);

    #if DEBUG_READMAPS > 0
    const size_t opN = opMapVec.size();
    const size_t cN = contigVec.size();
    cout << "Read " << opN << " optical maps: " << endl;
    for(size_t i = 0; i < opN; i++)
    {
        const OpticalMapData * pMap = opMapVec[i];
        cout << "\t" << pMap->getId() << ": "
                     << pMap->getNumFrags() << " frags, "
                     << pMap->getLength() << " bp."
                     << endl;
    }

    cout << "Read " << cN << " contigs: " << endl;
    for(size_t i = 0; i < cN; i++)
    {
        const ContigMapData * pMap = contigVec[i];
        cout << "\t" << pMap->getId() << ": "
                     << pMap->getNumFrags() << " frags, "
                     << pMap->getLength() << " bp."
                     << endl;
    }
    #endif
}

int main(int argc, char ** argv)
{
    // Parse command line arguments
    Parser::ArgParser * ap = Parser::ArgParser::instance();
    ap->parseArgs(argc, argv);
    ap->printArgs();
    
    Globals::initialize(); // Initialize dynamically created global variables

    //////////////////////////////////////////////////////////////////////////////
    // READ INPUT MAPS

    // Construct Optical Maps from file
    vector<OpticalMapData *> opticalMaps;
    vector<ContigMapData *> contigMaps;
    readMaps(opticalMaps, contigMaps);
    const int numContigs = contigMaps.size();
    //////////////////////////////////////////////////////////////////////////////

    //Setup output file streams
    XMLWriter xmlWriterAll, xmlWriterSig;
    xmlWriterAll.open((opt::outputPrefix + "_AllMatches.xml").c_str());
    xmlWriterSig.open((opt::outputPrefix + "_SigMatches.xml").c_str());

    ///////////////////////////////////////////////////////////////////////
    // START align
    //
    // Loop over contigs, starting with those with the most restriction sites
    MatchResultPtrVec allResults;
    MatchResultMap  matchResultMap;
    omp_set_num_threads(opt::numThreads);
    int reportPeriod = numContigs/100;

    int i;
    #pragma omp parallel for schedule(dynamic) private(i) shared(allResults, matchResultMap)
    for(i=0; i<numContigs; i++)
    {
        if (reportPeriod==0 || i%reportPeriod == 0)
            cout << "Aligning contig " << i << " of " << numContigs << "\n";

        ContigMapData * pContigMap = contigMaps[i];
        int numFrags = pContigMap->getNumFrags();

        if(numFrags <= 2 && !opt::oneToOneMatch)
            continue;

        MatchResultPtrVec resultList;
        matchContigToOpticalMaps(pContigMap, opticalMaps, &resultList, opt::localAlignment);
        sort(resultList.begin(), resultList.end(), resultScoreComp);
      
        if (resultList.size() > 0)
        {
            # pragma omp critical
            {
            allResults.insert(allResults.end(), resultList.begin(), resultList.end());
            matchResultMap.insert(make_pair(pContigMap, resultList));
            }
        }

        resultList.clear();
    } //end for loop over contigs

    cout << "Done Aligning Contigs\n";

    // If the permutation test is on, run the permutation test for those
    // contigs that have a match
    //runPermutationTests(&matchResultMap, opticalMaps, opt::numPermutationTrials);
    runPermutationTests2(&matchResultMap, opticalMaps, contigMaps, opt::numPermutationTrials);

    // DONE align
    ///////////////////////////////////////////////////////////////////////

    vector<MatchResult *>::iterator resultListIter = allResults.begin();
    vector<MatchResult *>::iterator resultListEnd = allResults.end();
    MatchResult * pResult;
    cout << "Writing alignments to file..." << endl;
    for(; resultListIter != resultListEnd; resultListIter++)
    {
        pResult = *resultListIter;
        pResult->annotate();
        xmlWriterAll.writeMatchResult(*pResult);
        if( pResult->pval_ <= opt::pThreshold)
            xmlWriterSig.writeMatchResult(*pResult);
        delete pResult; pResult = 0;
    }

    xmlWriterAll.close();
    xmlWriterSig.close();

    Globals::free(); // Free dynamically created global variables
   
    // Free ContigMapData memory
    {
        vector<ContigMapData *>::iterator vecIter, vecIterEnd;
        vecIter = contigMaps.begin();
        vecIterEnd = contigMaps.end();
        for(; vecIter!=vecIterEnd; vecIter++)
            delete *vecIter;
        contigMaps.clear();
    }
    // Free OpticalMapData memory
    {
        vector<OpticalMapData *>::iterator vecIter, vecIterEnd;
        vecIter = opticalMaps.begin();
        vecIterEnd = opticalMaps.end();
        for(; vecIter!=vecIterEnd; vecIter++)
            delete *vecIter;
        opticalMaps.clear();
    }

    return EXIT_SUCCESS;
}


// Run permutation test by permuting the optical map numTrials times, and then
// aligning each contig to the permuted maps
void runPermutationTests(MatchResultMap * pMatchResultMap, const vector<OpticalMapData *> opticalMapList,
                         int numTrials)
{
    if (numTrials < 1) return;

    // Create a pool of permuted optical maps
    vector<FragData> opticalFrags;
    vector<OpticalMapData *>::const_iterator oIter = opticalMapList.begin();
    vector<OpticalMapData *>::const_iterator oIterEnd = opticalMapList.end();
    for(; oIter!=oIterEnd ; oIter++)
    {
        const OpticalMapData * pOpticalMap = *oIter;
        const vector<FragData> * pFrags = &pOpticalMap->getFrags();
        opticalFrags.insert(opticalFrags.end(), pFrags->begin(), pFrags->end());
    }

    int numFrags = opticalFrags.size();
    vector< OpticalMapData * > permutedOpticalMaps;
    permutedOpticalMaps.reserve(numTrials);
    for(int i=0; i<numTrials; i++)
    {
        permute(opticalFrags);
        OpticalMapData * pPermutedMap = new OpticalMapData(numFrags, string(), false, opticalFrags);
        assert(pPermutedMap);
        permutedOpticalMaps.push_back(pPermutedMap);
    }

    AlignmentParams alignParams(opt::C_r_contig, opt::C_r_optical,
                                opt::sdMax, Constants::SIGMA2, opt::delta,
                                opt::smallFrag, opt::smallFragSlope,
                                opt::H, opt::T);

    // Loop over the contigs. For each contig, perform trials for the permutation test. 
    MatchResultMap::iterator matchMapIter = pMatchResultMap->begin();
    MatchResultMap::iterator matchMapEnd = pMatchResultMap->end();
    int numPermutationTests = pMatchResultMap->size();
    int reportPeriod = ((int) 0.01*numPermutationTests) + 1;
    int numProcessed = 0;
    for (; matchMapIter != matchMapEnd; matchMapIter++)
    {
        ContigMapData * pContigMap = matchMapIter->first;
        if (numProcessed%reportPeriod == 0)
            cout << "Processed permutation test " << numProcessed << " of " << numPermutationTests << "\n";
        vector<double> permutedScores = vector<double>(numTrials, -Constants::INF);

        #pragma omp parallel for schedule(dynamic) private(i)
        for (int i=0; i<numTrials; i++)
        {
            const OpticalMapData * pPermutedOpticalMap = permutedOpticalMaps[i];
            double scoreForward =  matchPermutationTest(pContigMap, pPermutedOpticalMap, true, alignParams);
            double scoreReverse =  matchPermutationTest(pContigMap, pPermutedOpticalMap, false, alignParams);
            permutedScores[i] = max(scoreForward, scoreReverse); // This appears to be thread safe.
        }

        // Sort permuted scores
        sort(permutedScores.begin(), permutedScores.end());
        const vector<double>::const_iterator psBegin = permutedScores.begin();
        const vector<double>::const_iterator psEnd = permutedScores.end();

        // Set pvalue for each of the results for this contig
        int numPermutedResultsBetter;
        vector<double>::const_iterator pLowerBound;
        vector<MatchResult *> * pResultList = &matchMapIter->second;
        vector<MatchResult *>::iterator it,itb,ite;
        MatchResult * pMatchResult;
        itb = pResultList->begin();
        ite = pResultList->end();
        for (it = itb; it != ite; it++)
        {
            pMatchResult = *it;
            // Get pointer to the first element which does not score less than this
            // MatchResult's score. Ties should give a larger pvalue!
            pLowerBound = lower_bound(psBegin, psEnd, pMatchResult->score_);
            numPermutedResultsBetter = psEnd - pLowerBound;
            pMatchResult->pval_ = ((double) numPermutedResultsBetter)/numTrials;
        }
        numProcessed++;
    }

    // Free memory
    for(int i=0; i<numTrials; i++) delete permutedOpticalMaps[i];
    permutedOpticalMaps.clear();
}

// Run permutation tests by selecting random contigs from the distribution of contig fragments,
// and aligning numTrials times to the original optical maps
void runPermutationTests2(MatchResultMap * pMatchResultMap, const vector<OpticalMapData *> opticalMapList,
                          const vector<ContigMapData *>& contigVec, int numTrials)
{

    if (numTrials < 1) return;

    // Populate vector of contig fragments
    vector<FragData> contigFrags;
    {
    vector<ContigMapData *>::const_iterator i = contigVec.begin();
    const vector<ContigMapData *>::const_iterator e = contigVec.end();
    vector<FragData>::const_iterator fragIter, fragEnd;
    for(; i != e; i++)
    {
       const ContigMapData * pMap = *i;
       const vector<FragData>& frags = pMap->getFrags();
       vector<FragData>::const_iterator fragEnd = frags.end();
       vector<FragData>::const_iterator fragStart = frags.begin();
       contigFrags.insert(contigFrags.end(), fragStart, fragEnd);
    }
    }

    const int N = contigFrags.size();

    cout << "Read " << N << " insilico contig fragments." << endl;
    const int numMaps = opticalMapList.size();

    set<int> fragNums; // Set of number of fragments in contigs that have a match
    {
    // Loop over the contigs. For each contig, perform trials for the permutation test. 
    MatchResultMap::iterator matchMapIter = pMatchResultMap->begin();
    MatchResultMap::iterator matchMapEnd = pMatchResultMap->end();
    for(; matchMapIter != matchMapEnd; matchMapIter++)
    {
        const ContigMapData * pMap = matchMapIter->first;
        fragNums.insert(pMap->getNumFrags());
    }
    }

    AlignmentParams alignParams(opt::C_r_contig, opt::C_r_optical,
                                opt::sdMax, Constants::SIGMA2, opt::delta,
                                opt::smallFrag, opt::smallFragSlope,
                                opt::H, opt::T);

    typedef vector<double> FloatVec; 
    typedef map<int, FloatVec> PermutationScoreMap;
    PermutationScoreMap permutationScoreMap;
    {
    set<int>::const_iterator i = fragNums.begin();
    set<int>::const_iterator e = fragNums.end();
    for(; i != e; i++)
        permutationScoreMap.insert(PermutationScoreMap::value_type(*i, FloatVec(numTrials, -Constants::INF)));
    }

    // For each fragment number, perform a permutation test
    {
    PermutationScoreMap::iterator i = permutationScoreMap.begin();
    PermutationScoreMap::iterator e = permutationScoreMap.end();
    int numRounds = permutationScoreMap.size();
    int reportPeriod = ((int) 0.01*numRounds) + 1;
    int numProcessed = 0;
    for(; i != e; i++)
    {
        if (numProcessed%reportPeriod == 0)
            cout << "Processed permutation test " << numProcessed << " of " << numRounds << endl;

        // Generate random contigs
        int numFrags = i->first;
        FloatVec& scores = i->second;
        vector<ContigMapData> permutedContigs(numTrials);
        for(int j = 0; j < numTrials; j++)
        {
            vector<FragData> frags(numFrags);
            for(int k=0; k < numFrags; k++)
                frags[k] =  contigFrags[rand() % N]; // random frag
            permutedContigs[j].setFrags(frags);
        }

        #pragma omp parallel for schedule(dynamic) private(j)
        for (int j=0; j < numTrials; j++)
        {
            // Match permuted contig to all optical maps
            const ContigMapData * pContigMap = &permutedContigs[j];
            double maxScore = -Constants::INF;
            for (int k=0; k<numMaps; k++)
            {
                const OpticalMapData * pOpMap = opticalMapList[k];
                double scoreForward =  matchPermutationTest(pContigMap, pOpMap, true, alignParams);
                double scoreReverse =  matchPermutationTest(pContigMap, pOpMap, false, alignParams);
                if (scoreForward > maxScore) maxScore = scoreForward;
                if (scoreReverse > maxScore) maxScore = scoreReverse;
            }
            scores[j] = maxScore;
        }
        sort(scores.begin(), scores.end()); // sort the scores
        numProcessed++;
    }
    }

    cout << "Assigning pvals...";
    MatchResultMap::iterator matchMapIter = pMatchResultMap->begin();
    MatchResultMap::iterator matchMapEnd = pMatchResultMap->end();
    for (; matchMapIter != matchMapEnd; matchMapIter++)
    {

        const ContigMapData * pMap = matchMapIter->first;
        int numFrags = pMap->getNumFrags();

        PermutationScoreMap::iterator pIter = permutationScoreMap.find(numFrags);
        assert(pIter != permutationScoreMap.end());
        const vector<double>& permutedScores = pIter->second;
        const vector<double>::const_iterator psBegin = permutedScores.begin();
        const vector<double>::const_iterator psEnd = permutedScores.end();

        // Set pvalue for each of the results for this contig
        int numPermutedResultsBetter;
        vector<double>::const_iterator pLowerBound;
        vector<MatchResult *> * pResultList = &matchMapIter->second;
        vector<MatchResult *>::iterator it,itb,ite;
        MatchResult * pMatchResult;
        itb = pResultList->begin();
        ite = pResultList->end();
        for (it = itb; it != ite; it++)
        {
            pMatchResult = *it;
            // Get pointer to the first element which does not score less than this
            // MatchResult's score. Ties should give a larger pvalue!
            pLowerBound = lower_bound(psBegin, psEnd, pMatchResult->score_);
            numPermutedResultsBetter = psEnd - pLowerBound;
            pMatchResult->pval_ = ((double) numPermutedResultsBetter)/numTrials;
        }
    }
    cout << "DONE." << endl;
}
// Run the match pipeline for a single contig to all of the optical maps
void matchContigToOpticalMaps(const ContigMapData * pContigMap, const vector<OpticalMapData *>& opticalMapList,
                              vector<MatchResult *> * const pResults, bool localAlignment)
{
    pResults->clear();

    int numContigFrags = pContigMap->getNumFrags();
    if (numContigFrags <= 2 && !opt::oneToOneMatch) return;
    int numOpticalMaps = opticalMapList.size();

    // Maximum number of missed restriction sites allowed for the contig
    double C_sigma = opt::sdMin;
   
    // Call dynampic programming algorithm to make match
    vector<MatchResult *> matches; // vector of all match results (including the best match)
    matches.reserve(1024);
    vector<MatchResult *> * pMatches = &matches;
    MatchResult * pResult = 0;
    MatchResult * pRevResult = 0;

    typedef MatchResult* (*MatchFunc)(const ContigMapData * pContigMap, const OpticalMapData * pOpticalMap,
                     vector<MatchResult *> * pMatches, bool forward, const AlignmentParams& alignParams);

    MatchFunc matchFunc = localAlignment ? matchLocal : match;

    // Align to all optical maps
    for (int i=0; i<numOpticalMaps; i++)
    {

        OpticalMapData * pOpticalMap = opticalMapList[i];

        // Course-grained search for alignment, relaxing the standard devation restriction
        C_sigma = opt::sdMin;
        while (C_sigma <= opt::sdMax)
        {
            AlignmentParams alignParams(opt::C_r_contig, opt::C_r_optical,
                                        C_sigma, Constants::SIGMA2, opt::delta,
                                        opt::smallFrag, opt::smallFragSlope,
                                        opt::H, opt::T);

            // Match in the forward direction
            bool forward = true;
            pResult = matchFunc(pContigMap, pOpticalMap, pMatches, forward, alignParams);

            if (!opt::noReverse)
            {
                // Match in the reverse direction
                forward = false;
                pRevResult = matchFunc(pContigMap, pOpticalMap, pMatches, forward, alignParams);
            }

            if( pResult || pRevResult)
                break;

            C_sigma += 1;
        }
    } // end for loop over optical maps

    if (matches.empty()) return;

    // sort MatchResults by score in descending order
    sort(matches.begin(), matches.end(), MatchResult::compareScore);
    reverse(matches.begin(), matches.end());

    // Only take the top opt::maxMatchesPerContig results;
    MatchResultPtrVec::iterator B = matches.begin();
    MatchResultPtrVec::iterator E = matches.end();
    MatchResultPtrVec::iterator E1;
    if (opt::maxMatchesPerContig > 0 && matches.size() > (size_t) opt::maxMatchesPerContig)
        E1 = matches.begin() + opt::maxMatchesPerContig;
    else
        E1 = matches.end();

    pResults->insert(pResults->end(), B, E1);

    for (MatchResultPtrVec::iterator it = E1; it != E; it++)
        delete *it;
    matches.clear();
}
