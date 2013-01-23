#############################
# File: SOMAMap.py
# Author: Lee Mendelowitz (lmendelo@umiacs.umd.edu)
#
# Description:
# Define the SOMA Map class for reading, writing, and creating restriction map
# objects.

class SOMAMap(object):
    def __init__(self, *args, **kwargs):
        if 'line' in kwargs:
            self.makeFromLine(kwargs['line'])
        else:
            self.makeFromAttributes(**kwargs)
        self.checkMap()

    # Create a SOMAMap from a line in maps file
    def makeFromLine(self, line):
        fields = line.strip().split()
        self.mapId = fields[0]
        self.length = int(fields[1])
        self.numFrags = int(fields[2])
        self.frags = [int(f) for f in fields[3:]]

    # Create from mapId and frags attribute
    def makeFromAttributes(self, **kwargs):
        self.frags = list(kwargs['frags'])
        self.mapId = kwargs['mapId']
        self.numFrags = len(self.frags)
        self.length = sum(self.frags)
       
        # Add any other attributes from kwargs 
        for attr in kwargs.iterkeys():
            if attr not in self.__dict__:
                self.__dict__[attr] = kwargs[attr]

    # write the SOMAMap to file
    def write(self, handle):
        f = handle
        fields = [self.mapId,
                  str(self.length),
                  str(len(self.frags))]
        fields.extend(str(frag) for frag in self.frags)
        outS = '\t'.join(fields) + '\n'
        f.write(outS)

    # Check the consistency of the object
    def checkMap(self):
        if len(self.frags) != self.numFrags:
            raise Exception('SOMAMap attributes are inconsistent!')


# Read map file in the SOMA map format
def readMaps(handle):
    fin = handle
    closeFile = False
    if type(fin) is str:
        fin = open(fin)
        closeFile = True
    maps = [SOMAMap(line=l) for l in fin]
    mapDict = dict((m.mapId, m) for m in maps)
    if closeFile:
        fin.close()
    return mapDict

# Wrate maps to a file in the SOMA map format
def writeMaps(mapList, handle):
    fout = handle
    closeFile = False
    if type(handle) is str:
        fout = open(handle, 'w')
        closeFile = True

    for map in mapList:
        map.write(fout)

    if closeFile:
        fout.close()
