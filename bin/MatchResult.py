import re
import numpy as np
from scipy.stats import chi2
from math import sqrt
import sys

# Count the number of instances of char c in string s
def countChar(s, c):
    return sum(1 for char in s if char == c)

def computeSd(c):
    sigma = 0.3*sqrt(1000)
    return sigma*sqrt(c)

def printAlignment(mr, fout=sys.stdout):
    mr.printAlignment(fout)

class MatchResult:
    def __init__(self, *args, **kwargs):
        assert ('xmlElement' in kwargs)
        xmlElement = kwargs['xmlElement']
        self.createFromXmlElement(xmlElement)

    def createFromXmlElement(self, xmlElement):
        # Extract contig information
        contigElement = xmlElement.find('contig');
        self.contigId = contigElement.find('name').text
        self.size = int(contigElement.find('size').text)

        # Chromosome name
        self.chromosomeName = xmlElement.find('chromosome').text

        # Extract alignment data
        self.forward = int(xmlElement.find('forward').text)
        self.score = float(xmlElement.find('score').text)
        self.pval = float(xmlElement.find('pval').text)
        self.opStart = int(xmlElement.find('opStartIndex').text)
        self.opEnd = int(xmlElement.find('opEndIndex').text)
        self.opStartBp = int(xmlElement.find('opStartBp').text)
        self.opEndBp = int(xmlElement.find('opEndBp').text)
        self.totalHits = int(xmlElement.find('totalHits').text)
        self.totalMisses = int(xmlElement.find('totalMisses').text)
        self.totalMissRate = float(xmlElement.find('totalMissRate').text)
        self.contigHits = int(xmlElement.find('contigHits').text)
        self.contigMisses = int(xmlElement.find('contigMisses').text)
        self.contigMissRate = float(xmlElement.find('contigMissRate').text)
        self.contigUnalignedBases = int(xmlElement.find('contigUnalignedBases').text)
        self.contigUnalignedBaseRatio = float(xmlElement.find('contigUnalignedBaseRatio').text)
        self.contigUnalignedFrags = int(xmlElement.find('contigUnalignedFrags').text)
        self.opticalHits = int(xmlElement.find('opticalHits').text)
        self.opticalMisses = int(xmlElement.find('opticalMisses').text)
        self.opticalMissRate = float(xmlElement.find('opticalMissRate').text)
        self.alignedLengthRatio = float(xmlElement.find('alignedLengthRatio').text)
        self.opticalMatchString = xmlElement.find('opticalMatchString').text
        self.contigMatchString = xmlElement.find('contigMatchString').text

        contigLostIndexText = xmlElement.find('contigLostIndex').text
        self.contigLostIndex = map(int, contigLostIndexText.split(',')) if contigLostIndexText else []

        self.alignedOpticalFrags = [np.array(sList.split(','),dtype=int) for sList in xmlElement.find('opticalAlignedIndex').text.split(';') if len(sList)]
        self.alignedContigFrags = [np.array(sList.split(','),dtype=int) for sList in xmlElement.find('contigAlignedIndex').text.split(';') if len(sList)]
        # Map the aligned fragments to numpy array indices
        offset = self.alignedOpticalFrags[0][0]
        self.alignedOpticalFragsInd = [np.array(fragList, dtype=int)-offset for fragList in self.alignedOpticalFrags]
        self.alignedContigFragsInd = [np.array(fragList, dtype=int) for fragList in self.alignedContigFrags] # Note: no offset here!

        self.unalignedContigFrags = self.contigLostIndex
        assert len(self.alignedOpticalFrags) == len(self.alignedContigFrags)

        # Compute additional information from the MatchStrings
        self.numAlignedFrags = len(self.alignedOpticalFrags)
        opticalFragPattern = re.compile(r'(\d+),-?(\d+)')
        fragDataList = opticalFragPattern.findall(self.opticalMatchString)
        self.opticalFragLengths = np.array([s[0] for s in fragDataList], dtype=int)
        self.opticalFragSD = np.array([s[1] for s in fragDataList], dtype=int)

        contigFragPattern = re.compile(r'(\d+)')
        contigDataList = contigFragPattern.findall(self.contigMatchString)
        self.contigFragLengths = np.array(contigDataList, dtype=int)
        self.contigFragSD = np.array([computeSd(c) for c in self.contigFragLengths])

        # Compute number of sites, hits, misses
        self.numContigSites = self.contigFragLengths.shape[0] - 1
        self.numOpticalSites = self.opticalFragLengths.shape[0] - 1
        self.numContigFragsUnaligned = len(self.unalignedContigFrags)

        # Check consistency of data:
        numContigHits = countChar(self.contigMatchString, ';')  # Number of matched restriction sites
        numContigMisses = sum([len(frags)-1 for frags in self.alignedContigFrags])
        numOpticalHits = countChar(self.opticalMatchString, ';')  # Number of matched restriction sites
        assert (numContigHits == self.contigHits)
        assert (numContigMisses == self.contigMisses)
        assert (numOpticalHits == self.opticalHits)
        #print self.contigId, self.chromosomeName, self.opStart, self.opEnd
        self.computeStats()

    # Compute statistics related to the alignment
    def computeStats(self):
        self.alignedOpticalLengths = np.array([np.sum(self.opticalFragLengths[frags]) for frags in self.alignedOpticalFragsInd])
        self.alignedOpticalVariance = np.array([np.sum(self.opticalFragSD[frags]**2) for frags in self.alignedOpticalFragsInd])
        self.alignedOpticalVariance = np.array([np.sum(self.opticalFragSD[frags]**2) for frags in self.alignedOpticalFragsInd])
        self.alignedContigVariance = np.array([np.sum(self.contigFragSD[frags]**2) for frags in self.alignedContigFragsInd])
        self.alignedContigLengths = np.array([np.sum(self.contigFragLengths[frags]) for frags in self.alignedContigFragsInd])
        # LMM 4/8: This data set does not have standard deviations reported for the optical fragments.
        # Avoid a divide by zero
        #self.alignedOpticalSD = (self.alignedContigLengths - self.alignedOpticalLengths)/np.sqrt(self.alignedOpticalVariance)
        self.alignedOpticalSD = np.zeros(self.numAlignedFrags)
        self.alignedContigSD = (self.alignedContigLengths - self.alignedOpticalLengths)/np.sqrt(self.alignedContigVariance)

        # Number of bases in a contig fragment that is a one-to-one alignment with optical fragment:
        oneToOneInd = np.array([ind[0] for ind in self.alignedContigFragsInd[1:-1] if ind.shape[0]==1], dtype=int)
        self.basesInOneToOne = np.sum(self.contigFragLengths[oneToOneInd])
        self.basesInOneToOneRatio = float(self.basesInOneToOne)/self.size
        self.medianAlignedSD = -1
        self.meanAlignedSD = -1
        self.minAlignedSD = -1
        self.maxAlignedSD = -1
        self.innerContigLength = -1
        self.innerOpticalLength = -1
        lengths = -1
        self.lengthRatio = -1
        self.alignedSD = self.alignedContigSD
        if self.numAlignedFrags > 2:
            self.medianAlignedSD = np.median(np.abs(self.alignedSD[1:-1]))
            self.meanAlignedSD = np.mean(np.abs(self.alignedSD[1:-1]))
            self.minAlignedSD = np.min(np.abs(self.alignedSD[1:-1]))
            self.maxAlignedSD = np.max(np.abs(self.alignedSD[1:-1]))
            self.innerContigLength = np.sum(self.alignedContigLengths[1:-1]) # Length of the inner aligned fragments
            self.innerOpticalLength = np.sum(self.alignedOpticalLengths[1:-1]) # Length of the inner aligned fragments
            lengths = np.array([self.innerContigLength, self.innerOpticalLength], dtype=float)
            self.lengthRatio = np.min(lengths)/np.max(lengths)
            assert len(self.alignedOpticalLengths) == len(self.alignedOpticalVariance) == len(self.alignedContigLengths)


    def printAlignment(self, fout=sys.stdout):
        mr = self
        fout.write('%s %i %f\n'%(mr.contigId, mr.size, mr.pval))
        fout.write('%s %i to %i\n'%(mr.chromosomeName, mr.opStart, mr.opEnd))
        numAligned = len(mr.alignedContigFragsInd)
        nextContigIndex = 0
        cStrings = []
        oStrings = []
        cStrings.append('Contig Frags')
        oStrings.append('Optical Frags')
        for i in range(numAligned):
            # Check if a gap needs to be printed
            firstContigIndex = mr.alignedContigFragsInd[i][0]
            if firstContigIndex != nextContigIndex:
                gapIndices = np.arange(nextContigIndex, firstContigIndex)
                gapContigFrags = mr.contigFragLengths[gapIndices]
                cfragString = ' + '.join(map(str, gapContigFrags))
                cfragString += ' = ' +  str(np.sum(gapContigFrags))
                ofragString = 'GAP'
                cStrings.append(cfragString)
                oStrings.append(ofragString)
                #fout.write('%25s | %25s\n'%(cfragString,ofragString))
            contigFrags = mr.contigFragLengths[mr.alignedContigFragsInd[i]]
            opticalFrags = mr.opticalFragLengths[mr.alignedOpticalFragsInd[i]]
            opticalFragSD = mr.opticalFragSD[mr.alignedOpticalFragsInd[i]]
            cfragString = ' + '.join(map(str, contigFrags))
            cfragString += ' = ' +  str(np.sum(contigFrags))
            ofragString = ' + '.join('%i'%(l) for l in opticalFrags)
            ofragString += ' = ' +  str(np.sum(opticalFrags))
            cStrings.append(cfragString)
            oStrings.append(ofragString)
            #fout.write('%25s | %25s\n'%(cfragString,ofragString))
            nextContigIndex = mr.alignedContigFragsInd[i][-1]+1

        cw = max(len(s) for s in cStrings)
        ow = max(len(s) for s in oStrings)
        # Fancy python string formatting:
        for c, o in zip(cStrings, oStrings):
            fout.write('{0:^{width1}} | {1:^{width2}}\n'.format(c, o, width1=cw, width2=ow))
