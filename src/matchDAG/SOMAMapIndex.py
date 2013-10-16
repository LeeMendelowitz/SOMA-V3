"""
Index a SOMAMap. 

This allows the efficient query of all fragments/indices that are in a specified range,
given by bp position.
"""
import numpy as np

class SOMAMapIndex(object):

    def __init__(self, map):
        self.map = map
        self.fragEndPos = np.cumsum(map.frags)
        self.fragStartpos = np.array([0] + list(self.fragEndPos[:-1]))
    
    def searchPos(self, minDist, maxDist):
        if minDist > maxDist:
            raise RuntimeError('minDist must be less than maxDist!')
        left = np.searchsorted(self.fragStartPos, minDist, side='left')
        right = np.searchsorted(self.fragEndPos, maxDist, side='right')
        return (left, right)
