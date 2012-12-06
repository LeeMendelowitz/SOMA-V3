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

using namespace std;

void readSilicoFile(const string& silicoFileName, vector<ContigMapData *>& retVec);
double missedSiteCostFunc(int d);
double missedFragmentCostFunc(int l);


ostream& operator<<(ostream& ss, const SiteData& siteData);

typedef map<string, ContigMapData *> ContigMap;
typedef map< string, vector<SiteData> > ContigSiteMap;
typedef map<string, int> ContigSizeMap;
typedef vector<ContigSiteMap::iterator> ContigSiteMapIterVec;

// Cost function for missing a restriction site
// whose closest neighboring site is d bp away.
double missedSiteCostFunc(int d);

// Cost function for missing a fragment of length l
double missedFragmentCostFunc(int l);

#endif
