########################################
# Filename: importSomaData.py
# Author: Lee Mendelowitz
# Date: 6/25/2012

# Tools to import data related to a SOMA alignment
############################################

import numpy as np
import cPickle
from MatchResult import MatchResult


##############################################
# Construct a Map from a line
class MapData(object):
    def __init__(self, line):
        fields = line.strip().split()
        self.mapId = fields[0]
        self.length = int(fields[1])
        self.numFrags = int(fields[2])
        self.frags = [int(f) for f in fields[3:]]
        assert(len(self.frags)==self.numFrags)

###############################################
def readMapFile(fname):
    fin = open(fname)
    maps = [MapData(l) for l in fin]
    fin.close()
    mapDict = dict((m.mapId, m) for m in maps)
    return mapDict


############################################
# Print a contig (type contigData)
def printContigFrags(contig):
    print contig.contigId
    print 'Size: %i numSites: %i'%(contig.length, contig.numSites)
    contigString = ['%i'%contig.frag[i] for i in range(contig.frag.shape[0])]
    contigString = '\n'.join(contigString)
    print contigString


############################################
# Read an optical map file
def readOpticalMap(mapFile):
    return readMapFile(mapFile)

############################################
# Read an insilico file and return a dictionary
# of contigIds to ContigData
def readSilicoFile(fileName):
    return readMapFile(fileName)

############################################
# Read pickle of list of MatchResults with
# unique matches (ie at most 1 MatchResult for contig)
def getContigToUniqueMatch(fileName):
    uniquePickleFile = open(fileName)
    ml = cPickle.load(uniquePickleFile)
    uniquePickleFile.close()
    contigToUniqueMatch = {}
    for mr in ml:
        contigToUniqueMatch[mr.contigId] = mr
    return contigToUniqueMatch

############################################
# Read pickle of list of MatchResults with
# (potentially) non-unique matches (ie allow for more than 1 MatchResult per contig)
def getContigToMatchList(fileName):
    pickleFile = open(fileName)
    ml = cPickle.load(pickleFile)
    pickleFile.close()
    contigToMatch = {}
    for mr in ml:
        if mr.contigId not in contigToMatch: contigToMatch[mr.contigId] = []
        contigToMatch[mr.contigId].append(mr)
    return contigToMatch
