// Lee Mendelowitz
// 10/6/2012
// Code to run and store permutation test results

#ifndef PERMUTATION_H
#define PERMUTATION_H

#include <vector>
#include <map>
#include <ostream>
#include <istream>
#include <set>

#include "mapTypes.h"

typedef unsigned int uint;
typedef std::vector<float> FloatVec;
typedef std::map<uint, FloatVec> UInt2FloatVec;

// Forward Declarations
class AlignmentParams;
//class FragData;
class ContigMapData;
class OpticalMapData;

// Class for storing permutation test results and getting a 
// p-value
class PermutationDB
{
    public:
        PermutationDB(unsigned int numTrials) : numTrials_(numTrials) {};
        PermutationDB(std::istream& is) {read(is);};

        // Clear the containers
        void reset() { scoreMap_.clear(); contigFrags_.clear(); N_=0;}

        // Run permutation tests
        void runTests(unsigned int maxN, const AlignmentParams& alignParams, const std::vector<OpticalMapData *>& opMaps);
        void runTests(const set<unsigned int>& numFrags, const AlignmentParams& alignParams, const vector<OpticalMapData *>& opMaps);

        // Get p-value
        float getPval(unsigned int numFrags, float score) const;

        // Populate vector of contig fragments. The FragData pointers are still
        // owned by the ContigMapData's.
        void setContigFrags(const std::vector<ContigMapData *>& contigVec);

        // Write out the distribution of scores 
        void write(std::ostream& os) const;

        // Read the distribution of scores
        void read(std::istream& is);

    private:
        unsigned int numTrials_;
        UInt2FloatVec scoreMap_;
        std::vector<FragData *> contigFrags_;
        size_t N_; // Number of fragments in contigFrags_

        std::vector<FragData *> generateRandomContig(unsigned int numFrags) const;
};
















#endif
