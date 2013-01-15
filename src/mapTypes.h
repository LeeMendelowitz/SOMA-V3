#ifndef MAPTYPES_H
#define MAPTYPES_H

#include <string>

using namespace std;

// Data for an in-silico (synthetic) restriction site
class SiteData
{
    public:
        int loc_; // restriction site location
        SiteData(int loc) :
            loc_(loc) {};
        SiteData() :
            loc_(0) { };
        SiteData(const SiteData& sd) :
            loc_(sd.loc_) {};
};

// Fragment data
class FragData
{
    public: 
        int size_; // Size in bp
        bool firstOrLastFrag_; // True if this fragment is the first or last fragment in a sequence

        FragData() :
            size_(0),
            firstOrLastFrag_(false)
        { };

        FragData(int size) :
            size_(size),
            firstOrLastFrag_(false)
        { };

        // Construct the fragment from two consecutive sites.
        // This is for constructing a contig restriction pattern
        FragData(const SiteData& sd1, const SiteData& sd2) : 
        firstOrLastFrag_(false)
        {
            size_ = sd2.loc_ - sd1.loc_;       
        }

        bool operator==(const FragData& fd) const
        {
            return (size_ == fd.size_);
        }
};

#endif
