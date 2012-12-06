#ifndef SCOREMATRIX_H
#define SCOREMATRIX_H

#include <cstdlib>
#include <sstream>

#include "globals.h"
#include "exception.h"

using namespace std;

// A cell inside the dynamic programming scoring table
class ScoreElement_t
{
    public:
    double score_;
    int pi_; // pointer to predecessor row
    int pj_; // pointer to predecessor col

    ScoreElement_t(double score, int pi, int pj) :
        score_(score),
        pi_(pi),
        pj_(pj) {};

    ScoreElement_t() :
        score_(-Constants::INF),
        pi_(-1),
        pj_(-1) {};
};

// The dynamic programming scoring table.
// The table is stored as a one dimensional array of ScoreElement_t's.
class ScoreMatrix_t
{
    public:

    ScoreElement_t * d_;
    const int m_; // numRows
    const int n_; // numCols

    ScoreMatrix_t(int m, int n) : 
        d_(0),
        m_(m),
        n_(n)
    {
        int numBytes = m*n*sizeof(ScoreElement_t);
        d_ = (ScoreElement_t * ) malloc (numBytes); 
        if (d_ == 0)
        {
            ostringstream msg;
            msg << "Could not allocate ScoreMatrix_t: "
                << " m=" << m
                << " n=" << n
                << " numBytes=" << numBytes << endl;
            throw Exception(msg.str());
        }
    };

    ~ScoreMatrix_t()
    {
        if (d_ != 0)
            free (d_);
    }
};

#endif
