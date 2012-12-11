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
#include "globals.h"
#include "parseArgs.h"
#include "utils.h"
#include "dp.h"
#include "match.h"
#include "xmlWriter.h"
#include "exception.h"
#include <omp.h>
#include <cstdlib>
#include <cassert>

//#include <boost/math/distributions.hpp> 
//#include "gamma-prob.c"

using namespace std;
using namespace Globals; // from globals.h

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) > (b) ? (b) : (a))


float betai(float a, float b, float x);
void filterResults(vector<MatchResult *>& results);
bool filterFunction(MatchResult * pMatchResult);

typedef vector<MatchResult *> MatchResultList;
typedef map<ContigMapData *, MatchResultList> MatchResultMap;

ostream& operator<<(ostream& ss, const SiteData& siteData)
{
    return (ss << siteData.loc_);
}

// Sort Algorithm Comparisons
bool contigDataVec_comp(const ContigMapData * pContig1, const ContigMapData * pContig2)
    { return pContig1->numFrags_ > pContig2->numFrags_; }
bool contigSiteDataVec_lt(const SiteData& sd1, const SiteData& sd2) {return sd1.loc_ < sd2.loc_;}
bool contigFragDataVec_lt(const FragData& fd1, const FragData& fd2) { return fd1.size_ < fd2.size_; }
bool resultScoreComp(const MatchResult * r1, const MatchResult * r2) { return r1->score_ > r2->score_; }


//MatchResult * runMatch(const ContigMapData * pContigMapData, const OpticalMapData * pOpticalMapData,
//              vector<MatchResult *> * const pOthers, const AlignmentParams& alignParams);

void matchContigToOpticalMaps(const ContigMapData * pContigMapData, const vector<OpticalMapData *>& opticalMapDataList,
                              vector<MatchResult *> * const pResults, bool getOthers = true, bool localAlignment = false);
void runPermutationTests(MatchResultMap * pMatchResultMap, const vector<OpticalMapData *> opticalMapDataList, int numTrials);
void runPermutationTests2(MatchResultMap * pMatchResultMap, const vector<OpticalMapData *> opticalMapDataList, const vector<ContigMapData *>& contigVec, int numTrials);

// Create a new vector replaced from orig where elements at positions
// are replaced by eleemnts from replacements
void replace(const vector<FragData>& orig, const vector<int>& positions,
                    const vector<FragData>& replacements, vector<FragData>& replaced) {

    replaced = orig;
    int n = positions.size();
    vector<int>::const_iterator iter = positions.begin();
    vector<int>::const_iterator end = positions.end();
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

// Read the contig silico file
// Create the ContigMapData instances, return pointers in a vector
void readSilicoFile(const string& silicoFileName, vector<ContigMapData *>& retVec)
{
    retVec.clear();
    ifstream silicoFile(silicoFileName.c_str(), ios_base::in);
    if (silicoFile.fail()) {
        ostringstream msg;
        msg << "ERROR: Could not open input file " << silicoFileName;
        throw(Exception(msg.str()));
    }

    string line, siteDataString;
    int siteLoc;
    set<string> contigIdSet; // set of contigId's
    istringstream issLine, issContigData;
    while(getline(silicoFile, line))
    {
        issLine.clear(); issLine.str(line); // clear the state flag bits
        string contigId;
        int length, numSites;
        vector<SiteData> contigSites;
        issLine >> contigId >> length >> numSites;
        assert(contigIdSet.count(contigId)==0); // Check that contigId is unique
        contigIdSet.insert(contigId);

        //Read the list of restriction site locations
        // The line that specifies the in-silico restriction site locations should separate each site by a semicolon
        getline(silicoFile, line);
        assert(!silicoFile.eof());
        issLine.clear(); issLine.str(line);
        while(getline(issLine, siteDataString, ';'))
        {
            siteLoc = stringToInt(siteDataString);
            contigSites.push_back(SiteData(siteLoc));
        }

        // Check that the correct number of sites were read
        assert(numSites == contigSites.size());

        // Sort the sites based on location
        sort(contigSites.begin(), contigSites.end(), contigSiteDataVec_lt);

        // Construct contigMapData
        ContigMapData * pContigMapData = new ContigMapData(length, contigId, contigSites);
        retVec.push_back(pContigMapData);
    }
    silicoFile.close();
    return;
}

int main(int argc, char ** argv)
{
    // Parse command line arguments
    Parser::ArgParser * ap = Parser::ArgParser::instance();
    ap->parseArgs(argc, argv);
    ap->printArgs();
    
    Globals::initialize(); // Initialize dynamically created global variables
    double size, sd;

    // Construct Optical Maps from file
    vector<OpticalMapData *> opticalMapDataList;
    int numOpticalMaps = Options::opticalMapList.size();
    for(int i=0; i<numOpticalMaps; i++)
    {
        string opticalMapFile = Options::opticalMapList[i];
        OpticalMapData * pOpticalMapData = 0;
        try {
            pOpticalMapData = new OpticalMapData(opticalMapFile, Options::circular);
        }
        catch (Exception& e)
        {
            cerr << e.what() << "\n";
            return(1);
        }
        opticalMapDataList.push_back(pOpticalMapData);
        cout << "Read Optical Map " << opticalMapFile << " with " << pOpticalMapData->numFrags_ << " fragments." << "\n";
    }

    // Read Silico file with contig map data    
    vector<ContigMapData *> contigVec;
    try {
        readSilicoFile(Options::silicoMap, contigVec);
    } 
    catch(Exception& e)
    {
        cerr << e.what()  << "\n";
        return(1);
    }
    const int numContigs = contigVec.size();
    cout << "Read Contig Silico file with " << numContigs << " contigs." << "\n";

    // Sort contigs in contigVec by the number of sites
    sort(contigVec.begin(), contigVec.end(), contigDataVec_comp);

    //Setup output file streams
    XMLWriter xmlWriterAll, xmlWriterSig;
    try {
        xmlWriterAll.open((Options::outputPrefix + "_AllMatches.xml").c_str());
        xmlWriterSig.open((Options::outputPrefix + "_SigMatches.xml").c_str());
    }
    catch (Exception& e)
    {
        cerr << e.what() << "\n";
        return(1);
    }

    // Loop over contigs, starting with those with the most restriction sites
    double C_sigma = Options::sdMin;
    const string emptyString = string("");
    
    MatchResultList allResults;
    MatchResultMap  matchResultMap;
    omp_set_num_threads(Options::numThreads);
    int reportPeriod = numContigs/100;


    #pragma omp parallel for schedule(dynamic) private(i) shared(allResults, matchResultMap)
    for(int i=0; i<numContigs; i++)
    {
        if (reportPeriod==0 || i%reportPeriod == 0)
            cout << "Aligning contig " << i << " of " << numContigs << "\n";

        ContigMapData * pContigMapData = contigVec[i];
        int numFrags = pContigMapData->numFrags_;

        if(numFrags <= 2 && !Options::oneToOneMatch)
            continue;

        MatchResultList resultList;
        bool getOthers = true;
        matchContigToOpticalMaps(pContigMapData, opticalMapDataList, &resultList, getOthers, Options::localAlignment);
        //cout << "Matching contig " << pContigMapData->contigId_ << ". Found " << resultList.size() << " results" << "\n";
        sort(resultList.begin(), resultList.end(), resultScoreComp);
      
        if (resultList.size() > 0)
        {
            # pragma omp critical
            {
            allResults.insert(allResults.end(), resultList.begin(), resultList.end());
            matchResultMap.insert(make_pair(pContigMapData, resultList));
            }
        }

        resultList.clear();
    } //end for loop over contigs

    cout << "Done Aligning Contigs\n";

    // If the permutation test is on, run the permutation test for those
    // contigs that have a match
    //runPermutationTests(&matchResultMap, opticalMapDataList, Options::numPermutationTrials);
    runPermutationTests2(&matchResultMap, opticalMapDataList, contigVec, Options::numPermutationTrials);

    vector<MatchResult *>::iterator resultListIter = allResults.begin();
    vector<MatchResult *>::iterator resultListEnd = allResults.end();
    MatchResult * pResult;
    cout << "Writing alignments to file..." << endl;
    for(; resultListIter != resultListEnd; resultListIter++)
    {
        pResult = *resultListIter;
        pResult->annotate();
        xmlWriterAll.writeAlignment(*pResult);
        if( pResult->pval_ <= Options::pThreshold)
            xmlWriterSig.writeAlignment(*pResult);
        delete pResult; pResult = 0;
    }

    xmlWriterAll.close();
    xmlWriterSig.close();

    Globals::free(); // Free dynamically created global variables
   
    // Free ContigMapData memory
    {
        vector<ContigMapData *>::iterator vecIter, vecIterEnd;
        vecIter = contigVec.begin();
        vecIterEnd = contigVec.end();
        for(; vecIter!=vecIterEnd; vecIter++)
            delete *vecIter;
        contigVec.clear();
    }
    // Free OpticalMapData memory
    {
        vector<OpticalMapData *>::iterator vecIter, vecIterEnd;
        vecIter = opticalMapDataList.begin();
        vecIterEnd = opticalMapDataList.end();
        for(; vecIter!=vecIterEnd; vecIter++)
            delete *vecIter;
        opticalMapDataList.clear();
    }

    return EXIT_SUCCESS;
}


// Run permutation test by permuting the optical map numTrials times, and then
// aligning each contig to the permuted maps
void runPermutationTests(MatchResultMap * pMatchResultMap, const vector<OpticalMapData *> opticalMapDataList,
                         int numTrials)
{
    if (numTrials < 1) return;

    // Create a pool of permuted optical maps
    vector<FragData> opticalFrags;
    vector<OpticalMapData *>::const_iterator oIter = opticalMapDataList.begin();
    vector<OpticalMapData *>::const_iterator oIterEnd = opticalMapDataList.end();
    for(; oIter!=oIterEnd ; oIter++)
    {
        const OpticalMapData * pOpticalMap = *oIter;
        const vector<FragData> * pFrags = &pOpticalMap->frags_;
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

    AlignmentParams alignParams(0, Options::C_r_contig, Options::C_r_optical,
                                Options::sdMax, Options::delta);  

    // Loop over the contigs. For each contig, perform trials for the permutation test. 
    MatchResultMap::iterator matchMapIter = pMatchResultMap->begin();
    MatchResultMap::iterator matchMapEnd = pMatchResultMap->end();
    int numPermutationTests = pMatchResultMap->size();
    int reportPeriod = ((int) 0.01*numPermutationTests) + 1;
    int numProcessed = 0;
    for (; matchMapIter != matchMapEnd; matchMapIter++)
    {
        ContigMapData * pContigMapData = matchMapIter->first;
        //cout << "Running permutation test for " << pContigMapData->contigId_ << "\n";
        if (numProcessed%reportPeriod == 0)
            cout << "Processed permutation test " << numProcessed << " of " << numPermutationTests << "\n";
        vector<double> permutedScores = vector<double>(numTrials, -Constants::INF);

        #pragma omp parallel for schedule(dynamic) private(i)
        for (int i=0; i<numTrials; i++)
        {
            const OpticalMapData * pPermutedOpticalMap = permutedOpticalMaps[i];
            double scoreForward =  matchPermutationTest(pContigMapData, pPermutedOpticalMap, true, alignParams);
            double scoreReverse =  matchPermutationTest(pContigMapData, pPermutedOpticalMap, false, alignParams);
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
void runPermutationTests2(MatchResultMap * pMatchResultMap, const vector<OpticalMapData *> opticalMapDataList,
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
       vector<FragData>::const_iterator fragEnd = pMap->frags_.end();
       vector<FragData>::const_iterator fragStart = pMap->frags_.begin();
       contigFrags.insert(contigFrags.end(), fragStart, fragEnd);
    }
    }

    const int N = contigFrags.size();

    cout << "Read " << N << " insilico contig fragments." << endl;
    const int numMaps = opticalMapDataList.size();

    set<int> fragNums; // Set of number of fragments in contigs that have a match
    {
    // Loop over the contigs. For each contig, perform trials for the permutation test. 
    MatchResultMap::iterator matchMapIter = pMatchResultMap->begin();
    MatchResultMap::iterator matchMapEnd = pMatchResultMap->end();
    for(; matchMapIter != matchMapEnd; matchMapIter++)
    {
        const ContigMapData * pMap = matchMapIter->first;
        fragNums.insert(pMap->frags_.size());
    }
    }

    AlignmentParams alignParams(0, Options::C_r_contig, Options::C_r_optical,
                                Options::sdMax, Options::delta);  

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
                const OpticalMapData * pOpMap = opticalMapDataList[k];
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
        int numFrags = pMap->frags_.size();

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
// If getAllMatches is true, add all matches to pResults. Otherwise, only add the best match to pResults.
void matchContigToOpticalMaps(const ContigMapData * pContigMapData, const vector<OpticalMapData *>& opticalMapDataList,
                              vector<MatchResult *> * const pResults, bool getAllMatches, bool localAlignment)
{
    pResults->clear();

    int numContigFrags = pContigMapData->frags_.size();
    if (numContigFrags <= 2 && !Options::oneToOneMatch) return;
    int numOpticalMaps = opticalMapDataList.size();

    // Maximum number of missed restriction sites allowed for the contig
    double C_sigma = Options::sdMin;
   
    // Call dynampic programming algorithm to make match
    vector<MatchResult *> allMatches; // vector of all match results (including the best match)
    allMatches.reserve(1024);
    vector<MatchResult *> * pAllMatches = getAllMatches ? &allMatches : 0; // if 0, then match() will only return the best match
    MatchResult * pResult = 0;
    MatchResult * pRevResult = 0;

    typedef MatchResult* (*MatchFunc)(const ContigMapData * pContigMapData, const OpticalMapData * pOpticalMapData,
                     vector<MatchResult *> * pAllMatches, bool forward, const AlignmentParams& alignParams);

    MatchFunc matchFunc = localAlignment ? matchLocal : match;

    if (matchFunc == matchLocal)
        cout << "Using match local" << std::endl;
    else if (matchFunc == match)
    {
        cout << "Using match" << std::endl;
    }
    else
    {
        cout << "Using unknown match function" << std::endl;
    }


    // Align to all optical maps
    for (int i=0; i<numOpticalMaps; i++)
    {

        OpticalMapData * pOpticalMapData = opticalMapDataList[i];
        int opticalMapSize = pOpticalMapData->numFrags_;

        // Course-grained search for alignment, relaxing the standard devation restriction
        C_sigma = Options::sdMin;
        while (C_sigma <= Options::sdMax)
        {
            AlignmentParams alignParams(0, Options::C_r_contig, Options::C_r_optical,
                                        C_sigma, Options::delta);  

            // Match in the forward direction
            bool forward = true;
            pResult = matchFunc(pContigMapData, pOpticalMapData, pAllMatches, forward, alignParams);
            if (pResult && !pAllMatches)
                allMatches.push_back(pResult);

            if (!Options::noReverse)
            {
                // Match in the reverse direction
                forward = false;
                pRevResult = matchFunc(pContigMapData, pOpticalMapData, pAllMatches, forward, alignParams);
                if (pRevResult && !pAllMatches)
                    allMatches.push_back(pRevResult);
            }

            if( pResult || pRevResult)
                break;

            C_sigma += 1;
        }

    } // end for loop over optical maps

    // sort MatchResults by score in descending order
    sort(allMatches.begin(), allMatches.end(), MatchResult::compareScore);
    reverse(allMatches.begin(), allMatches.end());

    if (getAllMatches)
    {
        // Place the filtered results into pResults vector
        filterResults(allMatches);
        pResults->insert(pResults->end(), allMatches.begin(), allMatches.end());
        allMatches.clear();
    }
    else
    {
        // Place only the best result into pResults vector if it passes the filter.
        vector<MatchResult *>::iterator it, ite;
        it = allMatches.begin();
        ite = allMatches.end();
        if (it != ite)
        {
            bool passFilter = filterFunction(*it);
            if (passFilter)
            {
                pResults->push_back(*it); // Add the first (best scoring) MatchResult
                it++;
            }
        }
        // Delete the rest
        for (; it != ite; it++)
            delete *it;
        allMatches.clear();
    }
}


// Return true if MatchResult is acceptable
bool filterFunction(MatchResult * pResult)
{
    // Build additional attributes for this alignment candidate, and see if it
    // passes the thresholds
    pResult->buildAlignmentAttributes();
   
    //Chi2 Filter
    bool failChi2Filter = false;
    if (pResult->numAlignedInnerBlocks_ > 0)
        failChi2Filter = (pResult->chi2_/((double) pResult->numAlignedInnerBlocks_)) > Options::avgChi2Threshold;
       
    bool failLengthRatio = pResult->alignedLengthRatio_ < Options::minLengthRatio;
    bool failContigMissRate = pResult->contigMissRate_ > Options::maxMissRateContig;
    bool failContigHitsCheck = pResult->contigHits_ <= 0;

    bool fail = (failChi2Filter || 
                 failLengthRatio ||
                 failContigMissRate ||
                 failContigHitsCheck);

    return (!fail);
}

// Filter results by:
// - thresholding on the alignment length ratio
// - selecting non-overlapping matches
// Delete any MatchResults's that are filtered out
// NOTE: the results should be sorted by score in descending order BEFORE
//       calling this function
void filterResults(vector<MatchResult *>& results)
{
    if (results.size()==0)
        return;

    vector<MatchResult *> filteredResults;
    vector<MatchResult *>::iterator it, itb, ite;
    vector<MatchResult *>::const_iterator frit, frb, fre;
    const MatchResult * pfr;
    MatchResult * pResult;
    bool hasOverlap;


    // Add results to filteredResults in descending order of score, provided
    // that each MatchResult satisfies the length ratio and miss rate criteria
    itb = results.begin();
    ite = results.end();
    int numAdded = 0;
    double lastScore = (*itb)->score_+1;
    for(it = itb;
       (it != ite) && numAdded < Options::maxMatchesPerContig;
       it++)
    {
        pResult = *it;

        // Check that scores are sorted properly
        assert (lastScore >= pResult->score_);
        lastScore = pResult->score_;

        // Check for overlap with a higher scoring MatchResult
        hasOverlap = false;
        frb = filteredResults.begin();
        fre = filteredResults.end();
        for (frit = frb; frit != fre; frit++)
        {
            pfr = *frit;
            if (pResult->overlaps(pfr))
            {
                hasOverlap = true;
                break;
            }
        }

        if (hasOverlap)
        {
            delete pResult;
            continue;
        }
        else
        {
            // Build additional attributes for this alignment candidate, and see if it
            // passes the thresholds
            pResult->buildAlignmentAttributes();
            bool passFilter = filterFunction(pResult); 
            if (!passFilter)
            {
                delete pResult;
                continue;
            }
            else
            {
                filteredResults.push_back(pResult);
                numAdded++;
            }
        }
    }

    // Cleanup the rest
    for(; it != ite; it++)
        delete *it;

    results = filteredResults;
}
