#############################################################################################
# File: make_opt.py
# Author: Lee Mendelowitz
# Date: 4/7/2012

# Description:
# Script to read optical maps in either the XML or Schwartz format
# The optical map data can be written to separate files for each chromosome (as SOMA expects)

# The soma optical map format has one record per line, tab delimited:
# "MAPID LENGTH_BP NUM_FRAGS FRAG1 FRAG2 FRAG3 ..."
#############################################################################################

import sys
import re

#############################################################################################
class OpticalMapData:
    def __init__(self, name, frags, fragSD, enzyme):
        self.name = name
        self.frags = frags
        self.length = sum(frags)
        self.fragSD = fragSD
        self.enzyme = enzyme
        if self.fragSD is not None:
            assert(len(self.fragSD)==len(self.frags))

    def writeMap(self, handle):
        f = handle

        fields = [self.name,
                  str(self.length),
                  str(len(self.frags))]
        fields.extend(str(frag) for frag in self.frags)

        outS = '\t'.join(fields) + '\n'
        f.write(outS)


#############################################################################################
# Read an optical map in the Schwartz format
# Return a list of OpticalMapData instances
def readMapDataSchwartz(filename):
    omaps = []
    fin = open(filename)
    for line in fin:

        # Skip if line is empty
        if not line.strip():
            continue

        line2 = fin.next()
        # Read two consecutive lines
        # Line 1: Chromosome name (this is line)
        # Line 2: Enzyme, 'B', List of fragment lengths (this is line2)
        chromName = line.strip()
        fields = line2.strip().split()
        enzyme = fields[0]
        fragLengths = [int(1000.0*float(field)) for field in fields[2:]] # Convert to bp
        omaps.append(OpticalMapData(chromName, fragLengths, None, enzyme))
        print 'Read map from chromosome %s with enzyme %s and %i fragments'%(chromName, enzyme, len(fragLengths))
    fin.close()
    return omaps

#############################################################################################
# Reads (in a crude fashion) an optical map in xml format
# Return a list of OpticalMapData instances
def readMapDataXML(fileName):
    fin = open(fileName)
    opticalMapList = []

    # Define regular expression for a line with map data
    # (frag_number) (size) (std. dev) (pfx) (sfx)
    # Note: pfx and sfx are optional!
    fragmentPattern = re.compile(r'S="(\d+)" STDDEV="([\d\.]+)"') # Finds fragment length and standard deviation
    startMapPattern = re.compile(r'<RESTRICTION_MAP ID="([^"]+)" ENZYME="([^"]+)"') # Start of a restriction map
    endMapPattern = re.compile(r'</RESTRICTION_MAP>') # End of a restriction map

    fragLengths = []
    fragSDs = []
    enzyme = ''
    mapName = ''
    for line in fin:
        startMapMatch = startMapPattern.search(line)
        endMapMatch = endMapPattern.search(line)
        fragmentMatch = fragmentPattern.search(line)
        # Assert that only at most one of these patterns has a match
        assert sum(map(int, [m != None for m in [startMapMatch, endMapMatch, fragmentMatch]])) <= 1

        if startMapMatch:
            # Assert that the previous map was ended properly
            assert(len(fragLength)==0)
            assert(len(fragSD)==0)
            assert(len(enzyme)==0)
            assert(len(mapName)==0)
            mapName = startMapMatch.group(1)
            enzyme = startMapMatch.group(2)
        elif fragmentMatch:
            fragLength = float(fragmentMatch.group(1))/1000 # in kbp
            fragSD = float(fragmentMatch.group(2)) # in kbp
            if fragSD==0.000:
                print 'Warning: Ignoring fragment with SD=0.000 in optical map %s'%mapName
            else:
                fragLengths.append(fragLength)
                fragSDs.append(fragSD)
        elif endMapMatch:
            # Create the OpticalMapData and add to the list
            print 'Creating optical map with name %s enzyme %s and %i fragments'%(mapName, enzyme, len(fragLengths))
            opticalMapList.append(OpticalMapData(mapName, fragLengths, fragSD, enzyme))
            fragLengths = []
            fragSDs = []
            enzyme = ''
            mapName = ''
    fin.close()

    # Assert that the previous map was ended properly
    assert(len(fragLengths)==0)
    assert(len(fragSDs)==0)
    assert(len(enzyme)==0)
    assert(len(mapName)==0)
    return opticalMapList

#############################################################################################
# Write a list of optical maps
def writeMaps(opticalMapDataList, outFileName):
    fout = open(outFileName, 'w')
    for opMap in opticalMapDataList:
        opMap.writeMap(fout)
    fout.close()
