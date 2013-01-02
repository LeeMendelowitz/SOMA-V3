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
    calcFragStartEnd();
}

void ContigMapData::setFrags(const vector<FragData>& frags)
{
    frags_ = frags;
    reverseContigFragDataVec(frags_, reverseFrags_);
    numFrags_ = frags.size();
    calcFragStartEnd();
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
void ContigMapData::calcFragStartEnd()
{
    fragStartBp_ = vector<int>(frags_.size(), 0);
    fragEndBp_ = vector<int>(frags_.size(), 0);
    int curPos = 0;
    // Set the start for the first fragment
    fragStartBp_[0] = 0;
    for(int i = 1; i < frags_.size()-1; i++)
    {
        curPos += frags_[i-1].size_;
        fragEndBp_[i-1] = curPos;
        fragStartBp_[i] = curPos;
    }
    // Set the end for the last fragment
    curPos += frags_[frags_.size()-1].size_;
    fragEndBp_[frags_.size()-1] = curPos;
}

int ContigMapData::getStartBp(int ind, bool forward) const
{
    if (forward)
    {
        return fragStartBp_[ind];
    }
    else
    {
        return fragEndBp_[frags_.size() - 1 - ind];
    }
}

int ContigMapData::getEndBp(int ind, bool forward) const
{
    if (forward)
    {
        return fragEndBp_[ind];
    }
    else
    {
        return fragStartBp_[frags_.size() - 1 - ind];
    }
}
