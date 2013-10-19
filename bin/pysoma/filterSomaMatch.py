"""
Filter a SOMA Match XML file
"""
# Author: Lee Mendelowitz
# Date: 10/19/2013
# LMendelo@umiacs.umd.edu

import xml.etree.cElementTree as ElementTree
ETree = ElementTree.ElementTree
from MatchResult import MatchResult
from FileWrapper import FileWrapper


#MAX_MISS_RATE = 0.10
MAX_SCORE_PER_CHUNK = 1.5


def makeFilterFunc(maxScorePerChunk = MAX_SCORE_PER_CHUNK):
    def filterPass(matchResult):
        numChunks = len(matchResult.alignment)
        scorePerChunk = abs(float(matchResult.score))/float(numChunks)
        if scorePerChunk <= maxScorePerChunk:
            return True
    return filterPass

def parse(fin):
    fin = FileWrapper(fin, 'r')
    for event, element in ElementTree.iterparse(fin):
        if element.tag == 'MatchResult':
            mr = MatchResult(xmlElement = element)
            yield mr
    fin.close()

def filter(fin, fout, filterFunc):
    def genMatches():
        for event, element in ElementTree.iterparse(fin):
            if element.tag == 'MatchResult':
                mr = MatchResult(xmlElement = element)
                if filterFunc(mr):
                    yield (element, mr)

    # Write matches to fout
    fout.write('<alignments>\n')
    for element, mr in genMatches():
        ETree(element).write(fout)
    fout.write('</alignments>\n')

def genMatchResultElements(fin):
    for event, element in ElementTree.iterparse(fin):
        if element.tag == 'MatchResult':
            yield(element)

def mergeXMLFiles(finList, fout):
    fout = FileWrapper(fout, 'w')
    fout.write('<?xml version="1.0" ?>\n')
    fout.write('<alignments>\n')
    for fin in finList:
        fin = FileWrapper(fin, 'r')
        for xmlElement in genMatchResultElements(fin):
            ETree(xmlElement).write(fout)
        fin.close()
    fout.write('</alignments>\n')
    fout.close()

def runFilter(fin, fout, maxScorePerChunk = MAX_SCORE_PER_CHUNK):
    fin = FileWrapper(fin, 'r')
    fout = FileWrapper(fout, 'w')

    filterFunc = makeFilterFunc(maxScorePerChunk)
    filter(fin, fout, filterFunc)

    fin.close()
    fout.close()
