#ifndef MAPTYPES_H
#define MAPTYPES_H

#include <string>

using namespace std;

// Data for an in-silico (synthetic) restriction site
class SiteData
{
    public:
        int loc_; // restriction site location
        double missCost_; // Cost for allowing this site to go unaligned
        SiteData(int loc) :
            loc_(loc), missCost_(0) {};
        SiteData(int loc, double missCost) :
            loc_(loc), missCost_(missCost) {};
        SiteData() :
            loc_(0), missCost_(0) {};
        SiteData(const SiteData& sd) :
            loc_(sd.loc_), missCost_(sd.missCost_) {};
};

// Fragment data
class FragData
{
    public: 
        int size_; // Size in bp
        int sd_;   // standard deviation (bp). (Not applicable to insilico fragment)
        bool used_; // true if this fragment has been used in an alignment  
        double missCost_; // Cost of missing this fragment in an alignment. (Applicable only to in silico fragment)
        int pos_; // Position in a sequence of fragments.
        bool firstOrLastFrag_; // True if this fragment is the first or last fragment in a sequence

        FragData() :
            size_(0), sd_(0), used_(false),
            missCost_(0.0), pos_(0), firstOrLastFrag_(false){};

        FragData(int size, int sd, int pos = -1) :
            size_(size), sd_(sd), used_(false),
            missCost_(0.0), pos_(pos), firstOrLastFrag_(false)
        {
            //missCost_ = missedFragmentCostFunc(size_);
        };

        // Construct the fragment from two consecutive sites.
        // This is for constructing a contig restriction pattern
        FragData(const SiteData& sd1, const SiteData& sd2, int pos = -1) : 
        pos_(pos), firstOrLastFrag_(false)
        {
            size_ = sd2.loc_ - sd1.loc_;       
            sd_ = 0; // SD is not applicable to an insilico fragment
            used_ = false;
            //missCost_ = missedFragmentCostFunc(size_);
        }

        bool operator==(const FragData& fd) const
        {
            return (size_ == fd.size_ &&
                    sd_ == fd.sd_ );
        }
};

#endif
