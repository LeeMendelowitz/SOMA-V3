#ifndef CONTIGMAPDATA_H
#define CONTIGMAPDATA_H

#include <string>
#include <vector>
#include <cassert>

#include "mapTypes.h"
#include "MapData.h"

using namespace std;

class ContigMapData : public MapData
{
    public:

    // Default Constructor
    ContigMapData();

    // Constructor
    ContigMapData(int length, const string& contigId, bool isForward);

    void setFrags(const vector<FragData>& frags);
    const vector<FragData>& getFrags() const
    {
        return frags_;
    }


    // Get the starting location of a restriction fragment
    int getStartBp(int ind, bool forward=true) const
    {
        assert((size_t) ind >= 0);
        assert((size_t) ind < frags_.size());

        if (forward)
            return fragBoundaryBp_[ind];

        // Contig is reverse, get the boundary position counting from the right
        size_t last = fragBoundaryBp_.size() - 1; size_t myInd = last - ind;
        return fragBoundaryBp_[myInd];
    }

    // Get the ending location of a restriction fragment
    int getEndBp(int ind, bool forward=true) const
    {
        assert((size_t) ind >= 0);
        assert((size_t) ind < frags_.size());

        if (forward)
            return fragBoundaryBp_[ind+1]; 

        // Contig is reverse, get the boundary position counting from the right
        size_t last = fragBoundaryBp_.size() - 1;
        size_t myInd = last - (ind + 1);
        return fragBoundaryBp_[myInd];
    }

    size_t getNumFrags() const { return frags_.size(); }

    int getLength() const { return length_; }

    bool isForward() const { return isForward_;}

    void setTwin(ContigMapData * pTwin) {pTwin_ = pTwin;}
    ContigMapData * getTwin() const { return pTwin_; }

    private:
    vector<int> fragBoundaryBp_; // indices of the boundaries for each contig fragment
    FragDataVec frags_;
    ContigMapData * pTwin_; // pointer to the reverse of this map
    int length_; // length in bp
    bool isForward_;

    void calcFragBoundaries();
};

void reverseContigFragDataVec(const vector<FragData>& orig, vector<FragData>& reversed);
bool readContigMaps(const std::string& fileName, vector<ContigMapData *>& contigMapVec);

#endif
