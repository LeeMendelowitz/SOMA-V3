#ifndef CONTIGMAPDATA_H
#define CONTIGMAPDATA_H

#include <string>
#include <vector>

#include "mapTypes.h"

using namespace std;

class ContigMapData
{
    public:
    int numFrags_;
    int length_; // The number of bases in contig
    string contigId_;
    vector<FragData> frags_;
    vector<FragData> reverseFrags_; // reverse of frags_

    // Default Constructor
    ContigMapData();

    // Constructor
    ContigMapData(int length, const string& contigId, const vector<SiteData>& sites);

    void getFrags(vector<FragData>& frags) const { frags = frags_;}
    void getReverseFrags(vector<FragData>& frags) const { frags = reverseFrags_;}
    void setFrags(const vector<FragData>& frags);

    private:

    // Filter out any silico fragments of zero size
    void computeFragsFromSites(const vector<SiteData>& sites);
};

void reverseContigFragDataVec(const vector<FragData>& orig, vector<FragData>& reversed);
#endif
