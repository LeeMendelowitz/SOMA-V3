#include "ScoreMatrixSeeded.h"
#include "globals.h"

#include <iostream>

using namespace seeded;

size_t ScoreMatrix::countFilledCells() {
    size_t count = 0;
    for (size_t i = 0; i < data_.size(); i++)
    {
        if (data_[i].score_ > -Constants::INF)
            count++;
        /*
         cout << "Score: " << data_[i].score_ << "\tBackPtrs: " << data_[i].backPointers_.size()
                                            << "\tForwardPtrs: " << data_[i].forwardPointers_.size()
                                            << endl;*/
    }
    return count;
}

double ScoreMatrix::getMaxScore() {
    double maxScore = -Constants::INF;
    for(size_t i = 0; i < numCols_; i++)
    {
        ScoreCell * pCell = getCell(numRows_-1, i);
        if (pCell->score_ > maxScore) 
            maxScore = pCell->score_;
    }
    cout << "Max Score: " << maxScore << "\n";
    return maxScore;
}

