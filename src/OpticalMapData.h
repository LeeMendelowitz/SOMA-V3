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

    int getStartBp(int ind, bool forward=true) const {
        assert( (size_t) ind < frags_.size() && ind >= 0);
        if (forward)
        {
            return fragStartBp_[ind];
        }
        return fragEndBp_[frags_.size() - 1 - ind];
    }

    int getEndBp(int ind, bool forward=true) const {
        assert( (size_t) ind < frags_.size() && ind >= 0);
        if (forward)
        {
            return fragEndBp_[ind];
        }
        return fragStartBp_[frags_.size() -1 - ind];
    }

    const vector<FragData>& getFrags(bool forward=true) const
    {
        if (forward) return frags_;
        else return reverseFrags_;
    }

    int numFrags_; // Number of fragments in map (without circular trick!)
    bool isCircular_;
    vector<FragData> frags_; // vector of fragments (doubled if circular)
    vector<FragData> reverseFrags_;

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
