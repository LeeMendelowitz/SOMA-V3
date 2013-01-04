#ifndef CONTIGMAPDATA_H
#define CONTIGMAPDATA_H

#include <string>
#include <vector>

#include "mapTypes.h"
#include "MapData.h"

using namespace std;

class ContigMapData : public MapData
{
    public:

    // Default Constructor
    ContigMapData();

    // Constructor
    ContigMapData(int length, const string& contigId, const vector<SiteData>& sites);

    const vector<FragData>& getFrags(bool forward=true) const
    {
        if(forward) return frags_;
        else return reverseFrags_;
    }

    void setFrags(const vector<FragData>& frags);

    // Get the starting and ending location of a restriction fragment
    int getStartBp(int ind, bool forward=true) const;
    int getEndBp(int ind, bool forward=true) const;

    int numFrags_;
    int length_; // The number of bases in contig
    vector<FragData> frags_;
    vector<FragData> reverseFrags_; // reverse of frags_

    private:
    vector<int> fragBoundaryBp_; // indices of the boundaries for each contig fragment

    // Filter out any silico fragments of zero size
    void computeFragsFromSites(const vector<SiteData>& sites);
    void calcFragBoundaries();

};

void reverseContigFragDataVec(const vector<FragData>& orig, vector<FragData>& reversed);
#endif
