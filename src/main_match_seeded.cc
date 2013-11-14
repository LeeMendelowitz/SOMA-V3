// match.cc

#include <vector>
#include <cstdlib>
#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
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
#include "ChunkDatabase.h"
#include "ScoreMatrixSeeded.h"
#include "Clock.h"

#define DEBUG_MATCH 0
#define DEBUG_FILTER 0
#define DEBUG_READMAPS 0

using namespace std;
using namespace Globals; // from globals.h

typedef vector<MatchResult *> MatchResultPtrVec;
typedef map<ContigMapData *, MatchResultPtrVec> MatchResultMap;

int main(int argc, char ** argv)
{

    Clock clock;
    std::cout.precision(5);
    std::cout << std::fixed;

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

    //////////////////////////////////////////
    // Set the map chunks for the references,
    // populate the ChunkDatabase, and 
    // sort the database.
    std::cout << "Building database on reference..." << std::flush;
    clock.lap();
    ChunkDatabase chunkDB;
    for(auto pRefMap : opticalMaps)
    {
        setMapChunks(pRefMap, opt::maxChunkMissesReference);
        chunkDB.addChunks(pRefMap->getChunks());
    }
    std::cout << "sorting..." << std::flush;
    chunkDB.sortFrags();

    std::cout << "done. "
              << clock.lap() << " seconds elapsed." << std::endl;


    //Setup output file streams
    XMLWriter xmlWriterAll;
    xmlWriterAll.open((opt::outputPrefix + ".matches.xml").c_str());

    ///////////////////////////////////////////////////////////////////////
    // START align
    // Loop over contigs, starting with those with the most restriction sites
    omp_set_num_threads(opt::numThreads);
    int reportPeriod = min(numContigs/100, 100);

    // Create the ScoreMatrix. Each thread will get its own private copy.
    seeded::ScoreMatrix scoreMatrix {};
    int i = 0;
    #pragma omp parallel for schedule(dynamic) private(i, scoreMatrix) shared(clock)
    for(i = 0; i<numContigs; i++)
    {
        if (reportPeriod > 0 && i%reportPeriod == 0)
        {
            # pragma omp critical
            {
            std::cout << "Aligning contig " << i << " of " << numContigs
                 << ". Elapsed: " << clock.lap() << " seconds\n";
            }
        }

        ContigMapData * pContigMap = forwardContigMaps[i];
        MatchResultPtrVec resultList;
        matchContigToOpticalMapsSeeded(pContigMap, chunkDB, scoreMatrix, &resultList);
      
        if (resultList.size() > 0)
        {
            // Annotate the matches
            for (auto pMatch : resultList)
                pMatch->annotate();

            // Critical block: Write the output to the XML file
            # pragma omp critical
            {
                clock.lap();
                // Write alignments to output file
                for (auto pMatch : resultList)
                    xmlWriterAll.writeMatchResult(*pMatch);
                std::cout << "Writing to xml. " << clock.lap() << " seconds.\n";
            }

            // Free the matches
            for (auto pMatch : resultList)
            {
                delete pMatch;
                pMatch = 0;
            }
            resultList.clear();
        }
    } //end for loop over contigs

    cout << "Done Aligning Contigs\n";

    // DONE align
    ///////////////////////////////////////////////////////////////////////

    xmlWriterAll.close();

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

    cout << clock.elapsed() << " seconds elapsed.\n";

    return EXIT_SUCCESS;
}
