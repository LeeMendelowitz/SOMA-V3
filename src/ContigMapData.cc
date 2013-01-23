#include <algorithm>
#include <iostream>

#include "ContigMapData.h"
#include "MapReader.h"
#include "utils.h"

#define READMAPS_DEBUG 0

void printMap(std::ostream& os, const ContigMapData * pMap);

// Default Constructor
ContigMapData::ContigMapData() :
    MapData(""),
    pTwin_(NULL),
    length_(0)
    {};

// Constructor
ContigMapData::ContigMapData(int length, const string& contigId, bool isForward) :
    MapData(contigId),
    pTwin_(NULL),
    length_(length),
    isForward_(isForward)
{ }

void ContigMapData::setFrags(const vector<FragData>& frags)
{
    frags_ = frags;
    calcFragBoundaries();

    if (!frags_.empty())
    {
        FragData * pFragData = &(*frags_.begin());
        pFragData->firstOrLastFrag_ = true;

        pFragData = &*(frags_.end() - 1);
        pFragData->firstOrLastFrag_ = true;
    }
}

// Reverse the orig vector<FragData>.
void reverseContigFragDataVec(const vector<FragData>& orig, vector<FragData>& reversed)
{
    reversed = orig;
    reverse(reversed.begin(), reversed.end());
}

// Calculate the starting and ending position of each restriction fragment (in bp)
// within the contig
void ContigMapData::calcFragBoundaries()
{
    const size_t N = frags_.size();
    fragBoundaryBp_ = vector<int>(N+1);
    int curPos = 0;
    for (size_t i = 0; i < N; i++)
    {
        curPos += frags_[i].size_;
        fragBoundaryBp_[i+1] = curPos;
    }
}

// Read contig maps from input file.
// Create ContigMaps both in the forward and reverse direction.
bool readMaps(const string& fileName, vector<ContigMapData *>& contigVec)
{
    MapReader reader(fileName);
    MapInput mapData;
    bool success = false;
    while(reader.next(mapData))
    {
        // Construct new maps
        ContigMapData * pForward = new ContigMapData(mapData.length_, mapData.mapId_, true);
        ContigMapData * pReverse = new ContigMapData(mapData.length_, mapData.mapId_ + ".r", false);

        pForward->setTwin(pReverse);
        pReverse->setTwin(pForward);

        // Set map fragments
        pForward->setFrags(mapData.frags_);
        reverse(mapData.frags_.begin(), mapData.frags_.end());
        pReverse->setFrags(mapData.frags_);
        contigVec.push_back(pForward);
        contigVec.push_back(pReverse);
        success = true;

        #if READMAPS_DEBUG > 0
        printMap(std::cout, pForward);
        printMap(std::cout, pReverse);
        #endif
    }
    return success;
}


void printMap(std::ostream& os, const ContigMapData * pMap)
{
    os << "Created ContigMap: " 
              << pMap->getId()
              << " frags: ";

    const vector<FragData>& frags = pMap->getFrags();
    for(vector<FragData>::const_iterator iter = frags.begin();
                                         iter != frags.end();
                                         iter++)
    {
        os << iter->size_ << " ";
    }
    os << std::endl;
}


