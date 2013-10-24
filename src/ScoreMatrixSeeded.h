#ifndef SCOREMATRIXSEEDED_H
#define SCOREMATRIXSEEDED_H
/*
This ScoreMatrix will maintain a vector to hold ScoreCell
objects. Ideally a single ScoreMatrix should be used per thread
for the lifetime of the thread.
    
A ScoreMatrix should beinitialized up-front with ScoreCells, and
ideally given a capacity large enough to handle any alignment problem
it will see.

The ScoreMatrix will also track which cells are in use by row.
*/

#include "ScoreCell.h"
#include "ScorePathStep.h"


namespace seeded {

class ScoreMatrix
{
    public:
    ScoreMatrix(size_t numRows, size_t numCols);
    void resize(size_t numRows, size_t numCols);
    size_t getSize() const { return size_; }
    size_t reset();
    size_t getCapacity() const { return data_.size(); }
    size_t getNumRows() { return numRows_; }
    size_t getNumCols() { return numCols_; }
    ScoreCell * getCell(size_t row, size_t col);
    ScoreCell * getCell(const IntPair& coord);
    ScoreCell * getCell(size_t coord);
    void addScorePathStep(ScorePathStep& step);
    void resetCells();

    private:
    size_t numRows_;
    size_t numCols_;
    size_t size_;
    ScoreCellVec data_;
    vector<ScoreCellSet> rowToCellsInPlay_;
    ScoreCellSet cellsInPlay_;
};

}



inline seeded::ScoreMatrix::ScoreMatrix(size_t numRows, size_t numCols) :
    numRows_(numRows),
    numCols_(numCols),
    size_(numRows_*numCols_)
{

}

inline void seeded::ScoreMatrix::resize(size_t numRows, size_t numCols) {
    numRows_ = numRows;
    numCols_ = numCols;
    size_ = numRows*numCols;
    if (size_ > data_.size())
    {
        data_.resize(size_);
    }

    if (rowToCellsInPlay_.size() < numRows_)
        rowToCellsInPlay_.resize(numRows_);

    resetCells();
}

inline void seeded::ScoreMatrix::resetCells() {
    for(size_t i = 0; i < numRows_; i++)
    {
        for(size_t j = 0; j < numCols_; j++)
        {
            ScoreCell * pCell = getCell(i, j);
            pCell->q_ = i;
            pCell->r_ = j;
            pCell->inPlay_ = false;
        }
    }

    // Mark the rest of the cells as invalid.
    for(size_t i = getSize(); i < data_.size(); i++)
    {
            ScoreCell * pCell = &data_[i];
            pCell->q_ = -1;
            pCell->r_ = -1;
            pCell->inPlay_ = false;
    }

    cellsInPlay_.clear();

    for(auto& cellSet : rowToCellsInPlay_)
        cellSet.clear();
}

inline ScoreCell * seeded::ScoreMatrix::getCell(size_t row, size_t col) {
    // Check that coords within bounds?
    return &data_[row*numCols_ + col];
}

inline ScoreCell * seeded::ScoreMatrix::getCell(const IntPair& coords)
{
    // Check that coords within bounds?
    return getCell(coords.first, coords.second);
}

inline ScoreCell * seeded::ScoreMatrix::getCell(size_t n) {
    // Check that coord within bounds?
    return &data_[n];
}

inline void seeded::ScoreMatrix::addScorePathStep(ScorePathStep& step)
{
    // Determine the query start/end coords
    ScoreCell * tgt = getCell(step.startCoord());
    tgt->forwardPointers_.push_back(step);

    ScoreCell * src = getCell(step.endCoord());
    src->backPointers_.push_back(step);

    step.setTarget(tgt);
    step.setSource(src);
    // Should we add step.src and step.tgt to the rowToCellsInPlay_?
}


#endif
