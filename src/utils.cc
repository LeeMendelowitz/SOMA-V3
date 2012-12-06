#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <assert.h>
#include "globals.h"
#include "utils.h"
#include "MersenneTwister.h"

using namespace std;

typedef pair<string, string> StringPair;
typedef map< const StringPair, int> DistanceMap;

// Given an input string, split the string using delimiter
// return vector of fields
void splitString(const string& input, char delim, vector<string> * pVecOutput)
{
    istringstream iss(input);
    string field;
    while(iss.good())
    {
        getline(iss, field, delim);
        pVecOutput->push_back(field);
    }
}

// Return the integer represented by the string
int stringToInt(const string& input)
{
    istringstream iss(input);
    int retVal;
    iss >> retVal;
    return retVal;
}
