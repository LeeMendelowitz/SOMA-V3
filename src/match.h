#ifndef MATCH_H
#define MATCH_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>

#include "exception.h"
#include "mapTypes.h"
#include "OpticalMapData.h"
#include "ContigMapData.h"
#include "MatchResult.h"

using namespace std;

void readSilicoFile(const string& silicoFileName, vector<ContigMapData *>& retVec);
double missedSiteCostFunc(int d);
double missedFragmentCostFunc(int l);

ostream& operator<<(ostream& ss, const SiteData& siteData);

typedef map<string, ContigMapData *> ContigMap;
typedef map< string, vector<SiteData> > ContigSiteMap;
typedef map<string, int> ContigSizeMap;
typedef vector<ContigSiteMap::iterator> ContigSiteMapIterVec;
typedef vector<MatchResult *> MatchResultPtrVec;
typedef map<ContigMapData *, MatchResultPtrVec> MatchResultMap;

// Cost function for missing a restriction site
// whose closest neighboring site is d bp away.
double missedSiteCostFunc(int d);

// Cost function for missing a fragment of length l
double missedFragmentCostFunc(int l);

bool resultScoreComp(const MatchResult * r1, const MatchResult * r2);

void matchContigToOpticalMaps(const ContigMapData * pContigMap, const vector<OpticalMapData *>& opticalMapList, vector<MatchResult *> * const pResults, bool localAlignment = false);

void runPermutationTests(MatchResultMap * pMatchResultMap, const vector<OpticalMapData *> opticalMapList, int numTrials);
void runPermutationTests2(MatchResultMap * pMatchResultMap, const vector<OpticalMapData *> opticalMapList, const vector<ContigMapData *>& contigVec, int numTrials);

// TO DO: Move this to another header & cc file
void readMaps(vector<OpticalMapData *>& opMapVec, vector<ContigMapData *>& allContigMaps, vector<ContigMapData *>& forwardContigMaps);


#endif
