##################################
# Author: Lee Mendelowitz (lmendelo@umiacs.umd.edu)
# Date: 1/22/2013
# File: significanceTest.py
#
# Description:

# Some code to generate random contig maps by sampling 
# with replacement from a distribution of contig fragments.
import numpy as np
import sys
import SOMAMap

##############################################################################
# Class to store fragments
class FragDatabase(object):

    def __init__(self, *args, **kwargs):
        self.interiorFrags = np.array([])
        self.allFrags = np.array([])
        if 'mapFile' in kwargs:
            sys.stdout.write('Adding frags from map %s...'%kwargs['mapFile'])
            self.addFragsFromMap(kwargs['mapFile'])
            sys.stdout.write('Added %i fragments\n'%len(self.allFrags))

    def addInteriorFrags(self, frags):
        self.interiorFrags = np.array(list(self.interiorFrags) + list(frags))

    def addFrags(self, frags):
        self.allFrags = np.array(list(self.allFrags) + list(frags))

    # Sample n restriction fragments (with replacement)
    def sampleFrags(self, n):
        ind = np.random.randint(0, self.allFrags.shape[0], n)
        return self.allFrags[ind]
        
    # Sample n interior restriction fragments (with replacement)
    def sampleInteriorFrags(self, n):
        ind = np.random.randint(0, self.interiorFrags.shape[0], n)
        return self.interiorFrags[ind]

    def addFragsFromMap(self, mapFileName):
        mapDict = SOMAMap.readMaps(open(mapFileName))
        for mapId, mapObj in mapDict.iteritems():
            frags = mapObj.frags
            if not frags:
                continue
            self.addFrags(frags)
            if len(frags) > 2:
                self.addInteriorFrags(frags[1:-1])
