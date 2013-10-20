#ifndef OPTICALMAPDATA_H
#define OPTICALMAPDATA_H

#include <string>
#include <vector>
#include <cassert>

#include "MapData.h"
#include "mapTypes.h"

using namespace std;

class OpticalMapData : public MapData
{
    public:

    OpticalMapData();
    OpticalMapData(int numFrags, const string& opticalId, bool isCircular, const vector<FragData>& frags);
    OpticalMapData(const string& mapFile, bool isCircular = 0);

    int getStartBp(int ind, bool forward=true) const
    {
        assert((size_t) ind >= 0);
        assert((size_t) ind < frags_.size());

        if (forward)
            return fragBoundaryBp_[ind];

        // Contig is reverse, get the boundary position counting from the right
        size_t last = fragBoundaryBp_.size() - 1;
        size_t myInd = last - ind;
        return fragBoundaryBp_[myInd];
    }

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
   
    // Get the number of fragments in the Optical map
    size_t getNumFrags() const { return numFrags_; }


    const vector<FragData>& getFrags() const
    {
        return frags_;
    }

    bool isCircular() const {return isCircular_;}

    int getLength() const { return length_; }

    FragDataVec::const_iterator getFragsB() const { return frags_.begin(); }
    FragDataVec::const_iterator getFragsE() const { return frags_.end(); }

    private:
    vector<int> fragBoundaryBp_; // The boundaries (in bp) of the restriction fragments in map
    FragDataVec frags_; // vector of fragments (doubled if circular)
    int numFrags_; // Number of fragments in map (without circular trick!)
    bool isCircular_;
    int length_; // Length of the map, in bp

    // Calculate the starting and ending position of each optical fragment (in bp)
    // within the optical map
    void calcFragBoundaries();
    void makeCircular();
    void readFile(const string&);
};

bool readOpticalMaps(const std::string fileName, std::vector<OpticalMapData *>& opMapVec);

#endif
