#ifndef DEBUGUTIL_H
#define DEBUGUTIL_H

#include "ScoreMatrix.h"
#include "MatchResult.h"
#include <fstream>
#include <sstream>
#include <ostream>

using namespace std;

std::ostream& operator<<(std::ostream& os, const ScoreElement_t * pE);
std::ostream& operator <<(std::ostream& os, const Index_t& index);

// Write a ScoreMatrix to file
void writeMatrixToFile(const ScoreMatrix_t * pScoreMatrix, const string& fileName);



#endif
