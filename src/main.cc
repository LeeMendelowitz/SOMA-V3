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
#include <ctime>

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
#define DEBUG_READMAPS 0

using namespace std;
using namespace Globals; // from globals.h

typedef vector<MatchResult *> MatchResultPtrVec;
typedef map<ContigMapData *, MatchResultPtrVec> MatchResultMap;

int main(int argc, char ** argv)
{

    clock_t start_time = std::clock();

    // Parse command line arguments
    Parser::ArgParser * ap = Parser::ArgParser::instance();
    ap->parseArgs(argc, argv);
    ap->printArgs();
    
    Globals::initialize(); // Initialize dynamically created global variables

    //////////////////////////////////////////////////////////////////////////////
    // READ INPUT MAPS

    // Construct Optical Maps from file
    vector<OpticalMapData *> opticalMaps;
    vector<ContigMapData *> allContigMaps;
    vector<ContigMapData *> forwardContigMaps;
    readMaps(opticalMaps, allContigMaps, forwardContigMaps);
    const int numContigs = forwardContigMaps.size();
    //////////////////////////////////////////////////////////////////////////////

    //Setup output file streams
    XMLWriter xmlWriterAll, xmlWriterSig;
    xmlWriterAll.open((opt::outputPrefix + "_AllMatches.xml").c_str());
    //xmlWriterSig.open((opt::outputPrefix + "_SigMatches.xml").c_str());

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

        ContigMapData * pContigMap = forwardContigMaps[i];
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
    runPermutationTests2(&matchResultMap, opticalMaps, forwardContigMaps, opt::numPermutationTrials);

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
        //if( pResult->pval_ <= opt::pThreshold)
            //xmlWriterSig.writeMatchResult(*pResult);
        delete pResult; pResult = 0;
    }

    xmlWriterAll.close();
    //xmlWriterSig.close();

    Globals::free(); // Free dynamically created global variables
   
    // Free ContigMapData memory
    {
        vector<ContigMapData *>::iterator vecIter, vecIterEnd;
        vecIter = allContigMaps.begin();
        vecIterEnd = allContigMaps.end();
        for(; vecIter!=vecIterEnd; vecIter++)
            delete *vecIter;
        allContigMaps.clear();
        forwardContigMaps.clear();
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

    clock_t end_time = std::clock();
    float delta = ((float) (end_time - start_time))/CLOCKS_PER_SEC;
    cout << delta << " seconds elapsed.\n";

    return EXIT_SUCCESS;
}
