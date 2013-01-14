#include <algorithm>

#include "ContigMapData.h"
#include "utils.h"


// Default Constructor
ContigMapData::ContigMapData() :
    MapData(""),
    pTwin_(NULL),
    length_(0)
    {};

// Constructor
// Create the frags_ from the sites vector provided
ContigMapData::ContigMapData(int length, const string& contigId, bool isForward, const vector<SiteData>& sites) :
    MapData(contigId),
    pTwin_(NULL),
    length_(length),
    isForward_(isForward)
{
    computeFragsFromSites(sites);
    calcFragBoundaries();
}

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
}

void ContigMapData::computeFragsFromSites(const vector<SiteData>& sites)
{
    frags_.clear();
    vector<SiteData>::const_iterator sIter = sites.begin();
    const vector<SiteData>::const_iterator sIterEnd = sites.end();

    // Create FragData from SiteData for first fragment
    int pos = 0;
    FragData fragData = FragData();
    fragData.size_ = sIter->loc_; fragData.sd_ = 0;
    fragData.pos_ = pos++;
    if (fragData.size_ > 0) frags_.push_back(fragData);
    sIter++;

    // Add the middle fragments
    for( ; sIter != sIterEnd; sIter++)
    {
        fragData = FragData(*(sIter-1), *sIter, pos++);
        if (fragData.size_>0) frags_.push_back(fragData);
    }
    
    // Add the last fragment
    vector<SiteData>::const_iterator pLastSite = sIterEnd-1;
    int lastFragLength = length_ - pLastSite->loc_;
    fragData = FragData();
    fragData.size_ = lastFragLength; fragData.sd_ = 0;
    fragData.pos_ = pos++;
    if (fragData.size_>0) frags_.push_back(fragData);
    size_t numFrags = frags_.size();
    frags_[0].firstOrLastFrag_ = true;
    frags_[numFrags-1].firstOrLastFrag_ = true;
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
