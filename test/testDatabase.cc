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

    return 0;
}
