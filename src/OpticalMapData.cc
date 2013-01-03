#include <fstream>
#include <sstream>

#include "OpticalMapData.h"
#include "exception.h"

using namespace std;

// Default Constructor
OpticalMapData::OpticalMapData() :
    numFrags_(0), opticalId_(""), isCircular_(false), frags_(vector<FragData>()) {};

// Constructor
OpticalMapData::OpticalMapData(int numFrags, const string& opticalId, bool isCircular, const vector<FragData>& frags) :
    numFrags_(numFrags), opticalId_(opticalId), isCircular_(isCircular), frags_(frags)
{
    calcFragStartEnd();
    if (isCircular_) makeCircular();
};

//Constructor (from file)
OpticalMapData::OpticalMapData(const string& mapFile, bool isCircular) :
     numFrags_(0), opticalId_(mapFile), isCircular_(isCircular) 
{
    readFile(mapFile);
    numFrags_ = frags_.size(); //Number of original fragments (before made circular)
    frags_[0].firstOrLastFrag_ = true;
    frags_[numFrags_-1].firstOrLastFrag_ = true;
    calcFragStartEnd();
    if (isCircular_) makeCircular();
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
void OpticalMapData::calcFragStartEnd()
{
    fragStartBp_ = vector<int>(numFrags_, 0);
    fragEndBp_ = vector<int>(numFrags_, 0);
    int curPos = 0;
    // Set the start for the first fragment
    fragStartBp_[0] = 0;
    for(int i = 1; i < numFrags_-1; i++)
    {
        curPos += frags_[i-1].size_;
        fragEndBp_[i-1] = curPos;
        fragStartBp_[i] = curPos;
    }
    // Set the end for the last fragment
    curPos += frags_[numFrags_-1].size_;
    fragEndBp_[numFrags_-1] = curPos;
}

void OpticalMapData::makeCircular()
{
    if (!isCircular_) return;

    // Double frags_
    frags_.insert(frags_.end(), frags_.begin(), frags_.end());
}
