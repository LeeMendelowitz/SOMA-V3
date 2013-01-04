#include <algorithm>

#include "ContigMapData.h"
#include "utils.h"


// Default Constructor
ContigMapData::ContigMapData() :
    MapData(""),
    numFrags_(0), length_(0) 
    {};

// Constructor
// Create the frags_ from the sites vector provided
ContigMapData::ContigMapData(int length, const string& contigId, const vector<SiteData>& sites) :
    MapData(contigId),
    length_(length)
{
    computeFragsFromSites(sites);
    calcFragBoundaries();
}

void ContigMapData::setFrags(const vector<FragData>& frags)
{
    frags_ = frags;
    reverseContigFragDataVec(frags_, reverseFrags_);
    numFrags_ = frags.size();
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
    //fragData.missCost_ = missedFragmentCostFunc(fragData.size_);
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
    //fragData.missCost_ = missedFragmentCostFunc(fragData.size_);
    fragData.pos_ = pos++;
    if (fragData.size_>0) frags_.push_back(fragData);
    numFrags_ = frags_.size();
    frags_[0].firstOrLastFrag_ = true;
    frags_[numFrags_-1].firstOrLastFrag_ = true;

    // Compute reverseFrags_
    reverseContigFragDataVec(frags_, reverseFrags_);
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

int ContigMapData::getStartBp(int ind, bool forward) const
{
    if (forward)
    {
        assert((size_t) ind < fragBoundaryBp_.size());
        assert(ind >= 0);
        return fragBoundaryBp_[ind];
    }
    else
    {
        size_t last = fragBoundaryBp_.size() - 1;
        size_t myInd = last - ind;
        assert(myInd < fragBoundaryBp_.size());
        assert(myInd >= 0);
        return fragBoundaryBp_[myInd];
    }
}

int ContigMapData::getEndBp(int ind, bool forward) const
{
    if (forward)
    {
        size_t myInd = ind+1;
        assert(myInd < fragBoundaryBp_.size());
        assert(myInd >= 0);
        return fragBoundaryBp_[myInd];
    }
    else
    {
        size_t last = fragBoundaryBp_.size() - 1;
        size_t myInd = last - (ind+1);
        assert(myInd < fragBoundaryBp_.size());
        assert(myInd >= 0);
        return fragBoundaryBp_[myInd];
    }
}
