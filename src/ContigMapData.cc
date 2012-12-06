#include <algorithm>

#include "ContigMapData.h"
#include "utils.h"


// Default Constructor
ContigMapData::ContigMapData() :
    numFrags_(0), length_(0), contigId_(""), frags_(vector<FragData>()) {};

// Constructor
// Create the frags_ from the sites vector provided
ContigMapData::ContigMapData(int length, const string& contigId, const vector<SiteData>& sites) :
    length_(length), contigId_(contigId)
{
    computeFragsFromSites(sites);
}

void ContigMapData::setFrags(const vector<FragData>& frags)
{
    frags_ = frags;
    reverseContigFragDataVec(frags_, reverseFrags_);
    numFrags_ = frags.size();
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
