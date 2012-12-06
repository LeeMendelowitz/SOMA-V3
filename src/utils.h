#ifndef UTILS_H
#define UTILS_H
#include <vector>
#include <string>
#include <map>
#include <assert.h>
#include "globals.h"
#include "utils.h"
#include "MersenneTwister.h"
#include <sstream>

using namespace std;

// Split a string using the given delimiter
void splitString(const string& input, char delim, vector<string> * pVecOutput);

// Return the integer represented by the string
int stringToInt(const string& input);

template <class T>
string toStr(const T& data)
{
    stringstream ss;
    ss << data;
    return ss.str();
}

template <class T>
const char * toCStr(const T&data)
{
    stringstream ss;
    ss << data;
    return ss.str().c_str();
}

// Permute a vector of types data_t
template<class data_t>
void permute(vector<data_t>& data) {
    MTRand generator;
    int s = data.size(); 
    for(int i = s-2; i >= 0; i--)
        swap(data[i], data[i+generator.randInt(s-i-1)]);
}

#endif
