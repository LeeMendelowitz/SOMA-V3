########################################
# Filename: importSomaData.py
# Author: Lee Mendelowitz
# Date: 6/25/2012

# Tools to import data related to a SOMA alignment
############################################

import numpy as np
import cPickle
from MatchResult import MatchResult

############################################
# Class to represent data about a contig
class ContigData:
    # Construct contigData from two lines of insilico file
    def __init__(self, lines):
        l1,l2 = lines

        # Read contig data from line 1
        self.contigId,length,numSites = l1.split()
        self.length = int(length); self.numSites = int(numSites)

        # Read contig data from line 2
        sites = [s.split(',')[0] for s in l2.split(';')]
        self.sites = np.array([int(s.strip()) for s in sites if len(s)])
        assert self.numSites == self.sites.shape[0]
    
        # Compute additional information about the contig
        self.frag = np.zeros(self.numSites + 1)
        self.frag[1:-1] = self.sites[1:] - self.sites[0:-1]
        self.frag[0] = self.sites[0] # first fragment
        self.frag[-1] = self.length-self.sites[-1] # last fragment
        self.numFrags = self.frag.shape[0]

############################################
# Print a contig (type contigData)
def printContigFrags(contig):
    print contig.contigId
    print 'Size: %i numSites: %i'%(contig.length, contig.numSites)
    contigString = ['%i'%contig.frag[i] for i in range(contig.frag.shape[0])]
    contigString = '\n'.join(contigString)
    print contigString


############################################
# Read a .opt optical map file
# Return a numpy array. Column 1 are fraglengths, col 2 are sd (in bp)
def readOpticalMap(mapFile):
    fin = open(mapFile)
    data = []
    for l in fin:
        fields = l.split()
        s = fields[0].strip()
        sd = fields[1].strip() if len(fields)>1 else -1.0
        data.append([float(s),float(sd)])
    fin.close()
    data = np.array(data)
    separators = data[:,0] == 0.001
    fragData = 1000*data[~separators,:] #Convert to bp
    return fragData

############################################
# Read an insilico file and return a dictionary
# of contigIds to ContigData
def readSilicoFile(fileName):
    fin = open(fileName)
    lines = fin.read().split('\n')
    fin.close()
    linePairs = [(lines[i-1], lines[i]) for i in range(1,len(lines),2)]
    contigDict = {}
    for lp in linePairs:
        cd = ContigData(lp)
        assert cd not in contigDict # Should be no duplicate contigs
        contigDict[cd.contigId] = cd
    return contigDict

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
