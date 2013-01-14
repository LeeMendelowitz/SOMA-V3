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


// Split a string using the given delimiter
void splitString(const std::string& input, char delim, std::vector<std::string> * pVecOutput);

// Return the integer represented by the string
int stringToInt(const std::string& input);

template <class T>
std::string toStr(const T& data)
{
    std::stringstream ss;
    ss << data;
    return ss.str();
}

// Convert a datatype to a const char array
// Note: The char array must be copied
// before a subsequent call to this function.
template <class T>
const char * toCStr(const T& data)
{
    static string s;
    std::stringstream ss;
    ss << data;
    s = ss.str();
    return s.c_str();
}

// Permute a vector of types data_t
template<class data_t>
void permute(std::vector<data_t>& data) {
    MTRand generator;
    int s = data.size(); 
    for(int i = s-2; i >= 0; i--)
        std::swap(data[i], data[i+generator.randInt(s-i-1)]);
}

#endif
