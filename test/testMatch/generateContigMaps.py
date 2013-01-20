import numpy as np

def parseOpticalMapFile(fname):
    fin = open(fname)
    mapList = []
    for l in fin:
        fields = l.split()
        mapId = fields[0]
        bp = int(fields[1])
        numFrags = int(fields[2])
        frags = [int(f) for f in fields[3:]]
        assert(len(frags)==numFrags)
        mapList.append((mapId, frags))
    fin.close()
    return mapList

# Define attributes of where a contig is sampled from:
#    - map, mapStart index, mapEnd index
#    - forward
class SampledContig(object):
    def __init__(self, frags, contigId, opMapId, opMapStartInd, opMapEndInd, contigIsForward):
        self.contigId = contigId 
        self.frags = frags
        self.opMapId = opMapId
        self.opMapStartInd = opMapStartInd
        self.opMapEndInd = opMapEndInd
        self.contigIsForward = contigIsForward

    def writeToMapFile(self, handle):
        numFrags = len(self.frags)
        bp = sum(self.frags)
        outS = [self.contigId, str(bp), str(numFrags)]
        outS.extend( str(f) for f in self.frags )
        outS = '\t'.join(outS)
        handle.write(outS + '\n')

    def writeToTruthFile(self, handle):
        outS = [self.contigId, self.opMapId, str(self.opMapStartInd), str(self.opMapEndInd), str(self.contigIsForward)]
        outS = '\t'.join(outS)
        handle.write(outS + '\n')


####################################################################################
NUM_CONTIGS = 20
MIN_FRAGS_CONTIG = 10
MAX_FRAGS_CONTIG = 20
ALLOW_REVERSE_CONTIG = True
ALLOW_FALSE_CUTS = True
FALSE_CUT_RATE = 0.1 # Number of false cuts to introduce per true site
ALLOW_MISSING_SITES = True
MISSING_SITE_RATE = 0.00 # Probability of missing a true site
MIN_CONTIG_FRAG = 500 # minimum frag length for randomly generated contig
MAX_CONTIG_FRAG = 70000 # max frag length for randomly generated contig

def makeContigFromOpticalMap(contigId, opMapId, opMapFrags):

    # Sample the contig frags
    contigNumFrags = np.random.randint(MIN_FRAGS_CONTIG, MAX_FRAGS_CONTIG)
    assert(contigNumFrags <= len(opMapFrags))
    contigStartLoc = np.random.randint(0, len(opMapFrags) - contigNumFrags)
    contigEndLoc = contigStartLoc + contigNumFrags
    contigFrags = opMapFrags[contigStartLoc:contigEndLoc]

    # Reverse contig frags, if allowed
    contigIsForward = True
    if ALLOW_REVERSE_CONTIG:
        if np.random.rand() > 0.5:
            contigFrags = contigFrags[::-1]
            contigIsForward = False

    # Convert the fragments to site locations.
    contigL = np.sum(contigFrags)
    sites = np.cumsum(contigFrags)[:-1]
    selectedSites = sites
    assert(len(sites) > 0)

    # Add false sites
    falseCutSites = []
    if (ALLOW_FALSE_CUTS):
        numFalseCuts = int(FALSE_CUT_RATE*len(sites))
        falseCutSites = list(np.random.randint(0, contigL, numFalseCuts))

    if (ALLOW_MISSING_SITES):
        probs = np.random.rand(len(sites))
        selectedSites = [s for s,p in zip(sites, probs) if p > MISSING_SITE_RATE]

    sites = sorted(selectedSites + falseCutSites) 
    contigFrags = list(np.diff(sites))
    contigFrags.append(contigL - sites[-1])

    sampledContig = SampledContig(contigFrags, contigId, opMapId, contigStartLoc, contigEndLoc, contigIsForward)
    return sampledContig


def makeRandomContig(contigId):
    contigNumFrags = np.random.randint(MIN_FRAGS_CONTIG, MAX_FRAGS_CONTIG)
    contigFrags = list(np.random.randint(MIN_CONTIG_FRAG, MAX_CONTIG_FRAG, contigNumFrags))
    contig = SampledContig(contigFrags, contigId, 'None', 'None', 'None', 'None')
    return contig



def run(opMapFileIn, contigMapFileOut, contigMapTruthOut, numTrueContigs, numRandomContigs):

    mapList = parseOpticalMapFile(opMapFileIn)

    contigList = []
    
    # Generate true contigs
    contigNum = 0
    for i in xrange(numTrueContigs):

        # Select an opMap at random
        opMapId, opMapFrags = mapList[np.random.randint(0,len(mapList))]
        contigId = 'contig-%i'%contigNum
        sampledContig = makeContigFromOpticalMap(contigId, opMapId, opMapFrags)
        contigList.append(sampledContig)
        contigNum += 1

    # Generate random contigs
    for i in xrange(numRandomContigs):
        contigId = 'contig-%i'%contigNum
        contigList.append(makeRandomContig(contigId))
        contigNum += 1

    # Write the output files
    fContigMap = open(contigMapFileOut, 'w')
    fContigMapTruth = open(contigMapTruthOut, 'w')

    for sampledContig in contigList:
        sampledContig.writeToMapFile(fContigMap)
        sampledContig.writeToTruthFile(fContigMapTruth)

    fContigMap.close()
    fContigMapTruth.close()


def main():
    opMapFile = 'map100.opt'
    contigFileOut = 'contigs.silico'
    contigTruthFileOut = 'contigs.truth'
    numTrueContigs = 10
    numRandomContigs = 10
    run(opMapFile, contigFileOut, contigTruthFileOut, numTrueContigs, numRandomContigs)

if __name__ == '__main__':
    main()
