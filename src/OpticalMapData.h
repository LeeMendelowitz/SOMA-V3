#ifndef OPTICALMAPDATA_H
#define OPTICALMAPDATA_H

#include <string>
#include <vector>
#include <cassert>

#include "mapTypes.h"

using namespace std;

class OpticalMapData
{
    public:

    OpticalMapData();
    OpticalMapData(int numFrags, const string& opticalId, bool isCircular, const vector<FragData>& frags);
    OpticalMapData(const string& mapFile, bool isCircular = 0);

    int getStartBp(int ind) const {
        assert( ind < numFrags_ && ind >= 0);
        return fragStartBp_[ind];
    }

    int getEndBp(int ind) const {
        assert( ind < numFrags_ && ind >= 0);
        return fragEndBp_[ind];
    }

    int numFrags_; // Number of fragments in map (without circular trick!)
    string opticalId_;
    bool isCircular_;
    vector<FragData> frags_; // vector of fragments (doubled if circular)

    private:
    vector<int> fragStartBp_; // starting location of fragment in bp
    vector<int> fragEndBp_; // ending location of fragment in bp

    // Calculate the starting and ending position of each optical fragment (in bp)
    // within the optical map
    void calcFragStartEnd();
    void makeCircular();
    void readFile(const string&);
};

#endif
