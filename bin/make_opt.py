#############################################################################################
# File: make_opt.py
# Author: Lee Mendelowitz
# Date: 4/7/2012

# Description:
# Script to read optical maps in either the XML or Schwartz format
# The optical map data can be written to separate files for each chromosome (as SOMA expects)
#############################################################################################

import sys
import re



#############################################################################################
class OpticalMapData:
    def __init__(self, name, fragLength, fragSD, enzyme):
        self.name = name
        self.fragLength = fragLength
        self.fragSD = fragSD
        self.enzyme = enzyme
        if self.fragSD is not None:
            assert(len(self.fragSD)==len(self.fragLength))

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
        fragLengths = map(float, fields[2:])
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
def writeMapDataSingleMap(opticalMapData, outFileName):
    fout = open(outFileName, 'w')
    hasSDData = (opticalMapData.fragSD != None)
    numFrags = len(opticalMapData.fragLength)
    for i in range(numFrags):
        sdString = '%.3f'%(opticalMapData.fragSD[i]) if hasSDData else ''
        fout.write('%.3f\t%s\n'%(opticalMapData.fragLength[i],sdString))
    fout.close()

#############################################################################################
# Write a list of optical maps
def writeMapDataAllMaps(opticalMapDataList, outFileName):
    fout = open(outFileName, 'w')
    firstMap = True
    for opMap in opticalMapDataList:
        # Write the separator before all maps except the first
        if not firstMap:
            fout.write('%.3f\t%.3f\n'%(0.001, 0.000))
        else:
            firstMap = False
        for fragL,fragSD in zip(opMap.fragLength, opMap.fragSD):
            fout.write('%.3f\t%.3f\n'%(fragL,fragSD))
    fout.close()

#############################################################################################
def main():
    mapFile = sys.argv[1]
    opMapList = readMapData(mapFile)
    
    # Remove all white space from restriction map names
    for opMap in opMapList:
        oldName = opMap.name
        opMap.name = ''.join(oldName.split())

    for opMap in opMapList:
        outFileName = opMap.name + '.opt'
        writeMapDataSingleMap(opMap, outFileName)


#############################################################################################
if __name__ == '__main__':
    main()
