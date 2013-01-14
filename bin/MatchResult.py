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

###################################################################
# Define constructors of objects from XML Elements

# The constructor is callable object which processes an xml element
# into an object of some type

# Create an object from the xmlElement text field
class TextConstructor(object):
    def __init__(self, constructor):
        self.constructor = constructor
    def __call__(self, xmlElement):
        return self.constructor(xmlElement.text)

IntField = TextConstructor(int)
FloatField = TextConstructor(float)
StrField = TextConstructor(str)
BoolField = TextConstructor(lambda s: bool(int(s)))

# Constructs a list of elements under the current element with the matching tag
# For example:
# <alignments>
#    <alignment>......</alignment>
#    <alignment>......</alignment>
#    <alignment>......</alignment>
# </alignments>
# Assume ae is Element object for the root node.
# ac = AlignmentConstructor() # This builds an alignment objects
# alc = ListConstructor(ac, 'alignment')
# myAlignments = alc(ae) # This returns a list of 3 alignment objects
class ListConstructor(object):
    def __init__(self, constructor, tag):
        self.constructor = constructor
        self.tag = tag

    def __call__(self, xmlElement):
        children = (e for e in xmlElement.iter() if e.tag==self.tag)
        #return [self.constructor(e) for e in xmlElement.iter(tag=self.tag)]
        return [self.constructor(e) for e in children]

#####################################################################

# A chunk from an optical map
class MapChunk(object):
    def __init__(self, xmlElement):

        # Create list of (key, constructor) pairs
        objList = [ ('start', IntField),
                    ('end', IntField),
                    ('startBp', IntField),
                    ('endBp', IntField),
                    ('frags', ListConstructor(IntField, 'frag'))
                  ]

        # Add attributes to the current object
        self.__dict__.update((key, constructor(xmlElement.find(key))) for key,constructor in objList)
    def lengthBp(self):
        return self.endBp - self.startBp


# The Score of a MatchedChunk
class Score(object):
    def __init__(self, xmlElement):
        objList = [ ('contig', FloatField),
                    ('optical', FloatField),
                    ('sizing', FloatField)
                  ]
        # Add attributes to the current object
        self.__dict__.update((key, constructor(xmlElement.find(key))) for key,constructor in objList)

# A matched chunk
class MatchedChunk(object):
    def __init__(self, xmlElement):
        objList = [ ('contig', MapChunk),
                    ('optical', MapChunk),
                    ('isContigGap', BoolField),
                    ('isBoundaryChunk', BoolField),
                    ('score', Score) ]
        # Add attributes to the current object
        self.__dict__.update((key, constructor(xmlElement.find(key))) for key,constructor in objList)

class MatchResult:
    def __init__(self, *args, **kwargs):
        assert ('xmlElement' in kwargs)
        xmlElement = kwargs['xmlElement']
        self.createFromXmlElement(xmlElement)

    def createFromXmlElement(self, xmlElement):

        # Extract contig information
        contigElement = xmlElement.find('contig');
        self.contigId = contigElement.find('name').text
        self.contigLength = int(contigElement.find('length').text)

        # Define attributes to be read from the xmlElement
        strFields = ['chromosome',
                     'opticalMatchString',
                     'opticalAlignedIndex',
                     'contigMatchString',
                     'scoreString'
                     ]
        intFields = ['cStartIndex', 'cEndIndex', 'cStartBp', 'cEndBp',
                     'cAlignedBases',
                     'opStartIndex',
                     'opEndIndex',
                     'opStartBp',
                     'opEndBp',
                     'opAlignedBases',
                     'totalHits',
                     'totalMisses',
                     'contigHits',
                     'contigMisses',
                     'contigUnalignedBases',
                     'contigUnalignedFrags',
                     'opticalHits',
                     'opticalMisses'
                     ]
        boolFields = ['forward']                     
        floatFields = ['score',
                       'pval',
                       'chi2',
                       'totalMissRate',
                       'contigMissRate',
                       'contigUnalignedBaseRatio',
                       'opticalMissRate'
                       ]

        objList = []
        objList.extend((k, StrField) for k in strFields)
        objList.extend((k, IntField) for k in intFields)
        objList.extend((k, BoolField) for k in boolFields)
        objList.extend((k, FloatField) for k in floatFields)
        objList.append( ('alignment', ListConstructor(MatchedChunk, 'chunk')) )

        # Add attributes to the current object
        self.__dict__.update( (key, constructor(xmlElement.find(key))) for key,constructor in objList )
                       

    # Compute statistics related to the alignment
    # 1/4/2013: This function is now incompatible with the current
    # implementation of the MatchResult class
    def computeStats(self):
        pass
        '''
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
        '''

    def printAlignment(self, fout=sys.stdout):
        mr = self
        fout.write('%s %i %f\n'%(mr.contigId, mr.contigLength, mr.pval))
        fout.write('%s %i to %i\n'%(mr.chromosome, mr.opStartIndex, mr.opEndIndex))
        cStrings = []
        oStrings = []
        cStrings.append('Contig Frags')
        oStrings.append('Optical Frags')

        # Iterate over the MatchedChunks
        for mc in self.alignment:
            if mc.isContigGap:
                cfragString = ' + '.join(map(str, mc.contig.frags))
                cfragString += ' = ' +  str(np.sum(mc.contig.frags))
                ofragString = 'GAP'
            else:
                cfragString = ' + '.join(map(str, mc.contig.frags))
                cfragString += ' = ' +  str(np.sum(mc.contig.frags))
                ofragString = str(np.sum(mc.optical.frags))
                ofragString += ' = ' + ' + '.join(map(str,mc.optical.frags))
            cStrings.append(cfragString)
            oStrings.append(ofragString)

        cw = max(len(s) for s in cStrings)
        ow = max(len(s) for s in oStrings)
        # Fancy python string formatting:
        for c, o in zip(cStrings, oStrings):
            fout.write('{0:^{width1}} | {1:^{width2}}\n'.format(c, o, width1=cw, width2=ow))
