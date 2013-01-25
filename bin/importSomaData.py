########################################
# Filename: importSomaData.py
# Author: Lee Mendelowitz
# Date: 6/25/2012

# Tools to import data related to a SOMA alignment
############################################

import numpy as np
import cPickle
from MatchResult import MatchResult
import SOMAMap
from collections import defaultdict


############################################
# Print a contig (type contigData)
def printContigFrags(contig):
    print contig.contigId
    print 'Size: %i numSites: %i'%(contig.length, contig.numSites)
    contigString = ['%i'%contig.frag[i] for i in range(contig.frag.shape[0])]
    contigString = '\n'.join(contigString)
    print contigString

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
    contigToMatch = defaultdict(list)
    for mr in ml:
        contigToMatch[mr.contigId].append(mr)
    return contigToMatch
