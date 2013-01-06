#ifndef DEBUGUTIL_H
#define DEBUGUTIL_H

#include "ScoreMatrix.h"
#include <fstream>
#include <sstream>
#include <ostream>

using namespace std;

std::ostream& operator<<(std::ostream& os, const ScoreElement_t * pE);

// Write a ScoreMatrix to file
void writeMatrixToFile(const ScoreMatrix_t * pScoreMatrix, const string& fileName);

#endif
