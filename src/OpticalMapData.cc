#include <fstream>
#include <sstream>
#include <algorithm>

#include "OpticalMapData.h"
#include "exception.h"

using namespace std;

// Default Constructor
OpticalMapData::OpticalMapData() :
    MapData(""),
    numFrags_(0),
    isCircular_(false)
    { };

// Constructor
OpticalMapData::OpticalMapData(int numFrags, const string& opticalId, bool isCircular, const vector<FragData>& frags) :
    MapData(opticalId),
    frags_(frags),
    numFrags_(numFrags),
    isCircular_(isCircular)
{
    calcFragBoundaries();
    if (isCircular_) makeCircular();
    reverseFrags_ = frags_;
    reverse(reverseFrags_.begin(), reverseFrags_.end());
};

//Constructor (from file)
OpticalMapData::OpticalMapData(const string& mapFile, bool isCircular) :
    MapData(mapFile),
    numFrags_(0),
    isCircular_(isCircular) 
{
    readFile(mapFile);
    numFrags_ = frags_.size(); //Number of original fragments (before made circular)
    frags_[0].firstOrLastFrag_ = true;
    frags_[numFrags_-1].firstOrLastFrag_ = true;
    calcFragBoundaries();
    if (isCircular_) makeCircular();
    reverseFrags_ = frags_;
    reverse(reverseFrags_.begin(), reverseFrags_.end());
}

void OpticalMapData::readFile(const string& mapFile)
{
    double sizef, sdf;
    int sizei, sdi;
    string line;
    istringstream iss;
    int pos = 0;

    frags_.clear();

    // Reading optical map data.
    ifstream mapFileStream(mapFile.c_str(), ios_base::in); // optical map file
    if(mapFileStream.fail())
    {
        ostringstream msg;
        msg << "Error: Could not open file " << mapFile << "!";
        throw(Exception(msg.str()));
        return;
    }

    while(getline(mapFileStream, line))
    {
        // parse line for size, sd
        sdf = 0;
        iss.clear(); iss.str(line); // must clear state flags
        iss >> sizef >> sdf;
        
        // convert size and standard deviation from kbp to bp
        sizei = (int) (1000*sizef + 0.5); // round to integer
        sdi = (int) (1000*sdf+0.5); // round to integer
        FragData myFrag = FragData(sizei, sdi, pos++);
        frags_.push_back(myFrag);
    }

    mapFileStream.close();
}
// Calculate the starting and ending position of each optical fragment (in bp)
// within the optical map
void OpticalMapData::calcFragBoundaries()
{
    const size_t N = frags_.size();
    fragBoundaryBp_ = vector<int>(N+1, 0);
    int curPos = 0;
    for (size_t i = 0; i < N; i++)
    {
        curPos += frags_[i].size_;
        fragBoundaryBp_[i+1] = curPos;
    }
}

void OpticalMapData::makeCircular()
{
    if (!isCircular_) return;

    // Double frags_
    frags_.insert(frags_.end(), frags_.begin(), frags_.end());
}
