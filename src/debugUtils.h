#include "ScoreMatrix.h"
#include <fstream>

using namespace std;

// Write a ScoreMatrix to file
void writeMatrixToFile(const ScoreMatrix_t * pScoreMatrix, const string& fileName)
{
    ofstream outFile;
    outFile.open(fileName.c_str());

    const int m = pScoreMatrix->m_;
    const int n = pScoreMatrix->n_;
    const ScoreElement_t * pE;

    for (int i = 0; i < m; i++)
        for (int j=0; j < n; j++)
        {
            pE = &(pScoreMatrix->d_[i*n + j]);
            outFile << "i: " << i
                    << " j: " << j 
                    << " score: " << pE->score_
                    << " pi: " << pE->pi_
                    << " pj: " << pE->pj_
                    << endl;
        }
    outFile.close();
}
