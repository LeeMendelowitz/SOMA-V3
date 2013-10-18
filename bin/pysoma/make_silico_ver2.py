#!/usr/bin/env python
"""
Make an insilico map for the given fasta file.
"""
import sys, os, time
from SOMAMap import SOMAMap


def genSites(seq, forwardRecSeq, reverseRecSeq = None):
    """
    Find the recoginition sites in sequence seq.
    recSeqs are the recognition sequences.
    They should only consist of letters ACTG
    """

    if not seq:
        return

    charSet = set('ACTG')
    for recSeq in (forwardRecSeq, reverseRecSeq):
        if recSeq is None:
            continue
        for c in recSeq:
            if c not in charSet:
                raise RuntimeError('Recognition sequences must only have alphabet ACTG')

    if reverseRecSeq is not None and reverseRecSeq == forwardRecSeq:
        raise RuntimeError('forwardRecSeq and reverseRecSeq should not match.')

    L = len(seq)
    if reverseRecSeq:
        # Check every position for a cut site.
        lf = len(forwardRecSeq)
        lr = len(reverseRecSeq)
        yield (0, 's')
        for i in xrange(L):
            if seq[i:i+lf] == forwardRecSeq:
                yield (i, 'f')
            if seq[i:i+lr] == reverseRecSeq:
                yield (i, 'r')
        yield (L, 'e')
    else:
        # Check every position for a cut site.
        lf = len(forwardRecSeq)
        yield (0, 's')
        for i in xrange(L):
            if seq[i:i+lf] == forwardRecSeq:
                yield (i, 'f')
        yield (L, 'e')


def makeInsilicoMap(seq, forwardRecSeq, reverseRecSeq = None, mapHandle = sys.stdout, siteHandle = None, mapName = 'map'):
    siteGen = genSites(seq, forwardRecSeq, reverseRecSeq)
    sites = [s for s in siteGen]

    def genFrags():
        i = iter(sites)
        lastSite = i.next()
        for s in i:
            yield (s[0] - lastSite[0])
            lastSite = s
    frags = [f for f in genFrags()]

    if siteHandle:
        siteHandle.write('\n'.join(str(s) for s in sites))

    numFrags = len(frags)
    map = SOMAMap(frags=frags, mapId = mapName)
    map.write(mapHandle)

    sys.stderr.write('Found %i restriction fragments.\n'%numFrags)

def test():
    """
    Test code.
    """
    fasta = 'chr1.head.fa'
    forwardRecSeq = 'GCTCTTC'
    reverseRecSeq = 'GAAGAGC'
    runFasta(fasta, forwardRecSeq, reverseRecSeq, 'testMap')

def runFasta(fastaFile, forwardRecSeq, reverseRecSeq, mapName):
    """
    Compute restriction fragments on a fasta file.
    """

    startTime = time.clock()

    path, fn = os.path.split(fastaFile)
    bn, ext = os.path.splitext(fn)

    from Bio import SeqIO
    seqGen = SeqIO.parse(fastaFile, 'fasta')
    rec = seqGen.next()
    seq = str(rec.seq)
    sys.stderr.write('Read Fasta file %s\n'%fastaFile)
    sys.stderr.write('%i bases read.\n'%len(seq))

    mapHandle = open('%s.map'%mapName, 'w')
    siteHandle = open('%s.sites'%mapName, 'w')
    makeInsilicoMap(seq, forwardRecSeq, reverseRecSeq, mapHandle, siteHandle, mapName)
    mapHandle.close()
    siteHandle.close()

    endTime = time.clock()
    sys.stderr.write('Run time: %.3f seconds.\n'%(endTime - startTime))


if __name__ == '__main__':
    test()
