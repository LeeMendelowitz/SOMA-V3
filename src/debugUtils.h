#ifndef DEBUGUTIL_H
#define DEBUGUTIL_H
#include "ScoreMatrix.h"
#include <fstream>
#include <sstream>
#include <ostream>

using namespace std;

std::ostream& operator<<(std::ostream& os, const ScoreElement_t * pE);
// Write a ScoreMatrix to file
void writeMatrixToFile(const ScoreMatrix_t * pScoreMatrix, const string& fileName)
{
    ofstream outFile;
    outFile.open(fileName.c_str());

    const int m = pScoreMatrix->m_;
    const int n = pScoreMatrix->n_;
    const ScoreElement_t * pE;
   
    // Write column headers
    outFile.width(5);
    outFile << " ";
    for (int j = 0; j < n; j++)
    {
        outFile.width(16); outFile << std::left << j;
    }
    outFile << std::endl;

    for (int i = 0; i < m; i++)
    {
        outFile.width(5);
        outFile << i;
        for (int j=0; j < n; j++)
        {
            pE = &(pScoreMatrix->d_[i*n + j]);
            /*
            outFile << "i: " << i
                    << " j: " << j 
                    << " score: " << pE->score_
                    << " pi: " << pE->pi_
                    << " pj: " << pE->pj_
                    << endl;
            */
            outFile << pE;
        }
        outFile << std::endl;
    }
    outFile.close();
}

std::ostream& operator<<(std::ostream& os, const ScoreElement_t * pE)
{
    std::ostringstream oss;
    oss.precision(2);
    oss << "(" << std::fixed << pE->score_
        << "," << pE->pi_
        << "," << pE->pj_
        << ")";
    os.width(16);
    os << std::left << oss.str();
    return os;
}

#endif
