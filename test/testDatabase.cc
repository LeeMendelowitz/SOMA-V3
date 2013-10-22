#include <vector>
#include <iostream>

#include "ContigMapData.h"
#include "OpticalMapData.h"
#include "ChunkDatabase.h"
#include "MapChunkUtils.h"

int main()
{
    using namespace std;

    vector<ContigMapData *> contigMaps;
    vector<OpticalMapData *> opticalMaps;

    // Read input maps
    readOpticalMaps("map100.opt", opticalMaps);
    readContigMaps("contigs.silico", contigMaps);
    cout << "Read " << contigMaps.size() << " contig maps.\n";
    cout << "Read " << opticalMaps.size() << " optical maps.\n" << flush;

    // Build fragment database
    const size_t missedInteriorSites = 2;
    for (auto omap : opticalMaps)
    {
        cout << "setting map chunks for map " << omap->getId() << endl;
        setMapChunks(omap, missedInteriorSites);
    }
    cout << "done setting maps!" << endl;
    ChunkDatabase chunkDB;

    cout << "Adding map..." << endl;
    chunkDB.addChunks(opticalMaps.front()->getChunks());
    cout << "sorting..." << endl;
    chunkDB.sortFrags();

    cout << "Done sorting optical map fragments in database.\n";

    // Make a query
    IntPairVec bounds;
    // 40041   33760   49438
    bounds.push_back(IntPair(30000,50000));
    bounds.push_back(IntPair(30000,50000));
    bounds.push_back(IntPair(30000,50000));

    vector<MapChunk*> hits;
    chunkDB.getMapChunkHits(bounds, hits);

    cout << "Found " << hits.size() << " hits.\n";
    for(size_t i = 0; i < hits.size(); i++)
    {
        cout << *(hits[i]) << "\n";
    }

    bounds.resize(1);
    chunkDB.getMapChunkHits(bounds, hits);

    cout << "Found " << hits.size() << " hits.\n";
    for(size_t i = 0; i < hits.size(); i++)
    {
        cout << *(hits[i]) << "\n";
    }

    return 0;
}
