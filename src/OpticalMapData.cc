#include <fstream>
#include <sstream>
#include <algorithm>

#include "OpticalMapData.h"
#include "MapReader.h"
#include "exception.h"
#include "globals.h"

using namespace std;

// Default Constructor
OpticalMapData::OpticalMapData() :
    MapData(""),
    numFrags_(0),
    isCircular_(false),
    length_(0)
    { }

// Constructor
OpticalMapData::OpticalMapData(int numFrags, const string& opticalId, bool isCircular, const vector<FragData>& frags) :
    MapData(opticalId),
    frags_(frags),
    numFrags_(numFrags),
    isCircular_(isCircular)
{
    calcFragBoundaries();
    if (isCircular_) makeCircular();
}

//Constructor (from file)
OpticalMapData::OpticalMapData(const string& mapFile, bool isCircular) :
    MapData(mapFile),
    numFrags_(0),
    isCircular_(isCircular),
    length_(0)
{
    readFile(mapFile);
    numFrags_ = frags_.size(); //Number of original fragments (before made circular)
    frags_[0].firstOrLastFrag_ = true;
    frags_[numFrags_-1].firstOrLastFrag_ = true;
    calcFragBoundaries();
    if (isCircular_) makeCircular();
}

void OpticalMapData::readFile(const string& mapFile)
{
    double sizef, sdf;
    int sizei;
    string line;
    istringstream iss;

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
        FragData myFrag = FragData(sizei);
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
    length_ = curPos;
}

void OpticalMapData::makeCircular()
{
    if (!isCircular_) return;

    // Double frags_
    frags_.insert(frags_.end(), frags_.begin(), frags_.end());
}


bool readOpticalMaps(const std::string fileName, std::vector<OpticalMapData *>& opMapVec)
{
    MapReader reader(fileName);
    MapInput mapData;
    bool success = false;
    while(reader.next(mapData))
    {
        OpticalMapData * pMap = new OpticalMapData(mapData.frags_.size(), mapData.mapId_, opt::circular, mapData.frags_);
        opMapVec.push_back(pMap);
        success = true;
    }
    return success;
}
