import numpy as np

#######################################
class OpticalMap(object):

    # Frags is a list of fragment lengths
    def __init__(self, mapId, frags):
        self.mapId = mapId
        self.frags = frags
        self.length = sum(self.frags)

    def writeMap(self, handle):
        outS = [self.mapId, str(self.length), str(len(self.frags))]
        outS.extend(str(frag) for frag in self.frags)
        outS = '\t'.join(outS)
        handle.write(outS + '\n')

#######################################
def generateFrags(numFrags, minFragLength, maxFragLength):
    frags = list(np.random.randint(minFragLength, maxFragLength, numFrags))
    return frags


def main():
    outFile = 'map100.opt'
    fout = open(outFile, 'w')
    numFrags = 500
    minFragL = 500
    maxFragL = 70000
    frags = generateFrags(numFrags, minFragL, maxFragL)
    opMap = OpticalMap('opmap1', frags)
    opMap.writeMap(fout)

    fout.close()


if __name__ == '__main__':
    main()
