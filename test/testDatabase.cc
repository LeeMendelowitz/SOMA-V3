#include <vector>
#include <iostream>

#include "ContigMapData.h"
#include "OpticalMapData.h"
#include "FragDatabase.h"

int main()
{
    using namespace std;

    vector<ContigMapData *> contigMaps;
    vector<OpticalMapData *> opticalMaps;

    // Read input maps
    readOpticalMaps("map100.opt", opticalMaps);
    readContigMaps("contigs.silico", contigMaps);
    cout << "Read " << contigMaps.size() << " contig maps.\n";
    cout << "Read " << opticalMaps.size() << " optical maps.\n";

    // Build fragment database
    FragDatabase fragDB;
    fragDB.addMap(opticalMaps.front());
    fragDB.sortFrags();

    cout << "Done sorting optical map fragments in database.\n";

    // Make a query
    IntPairVec bounds;
    // 40041   33760   49438
    bounds.push_back(IntPair(30000,50000));
    bounds.push_back(IntPair(30000,50000));
    bounds.push_back(IntPair(30000,50000));

    vector<FragPtr*> hits;
    fragDB.getFragPtrHits(bounds, hits);

    cout << "Found " << hits.size() << " hits.\n";
    for(size_t i = 0; i < hits.size(); i++)
    {
        cout << *(hits[i]) << "\n";
    }

    return 0;
}
