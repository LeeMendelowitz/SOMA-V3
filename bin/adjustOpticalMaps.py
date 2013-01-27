#!/usr/bin/env python
# adjustOpticalMaps.py
# Author: Lee Mendelowitz (lmendelo@umiacs.umd.edu)
# Date: 1/26/2013 
#
# Wrapper around adjustOpticalMaps.R
# Run the R script to adjust optical map fragments
# to remove the oversizing bias from a finished optical map.
# The bias is removed by using loess regression to predict the
# true restriction fragment size from the reported optical map
# restriction fragment size.

import argparse
import subprocess
import os
import sys

parser = argparse.ArgumentParser(description="adjustOpticalMaps: Remove oversizing bias from optical maps.")
parser.add_argument('inputMapFile', help='path to the input optical map file')
parser.add_argument('matchFragFile', help='path to the training data for loess regression')
parser.add_argument('outputMapFile', help='output optical map file')

R_SCRIPT = 'adjustOpticalMaps.R'

def MissingInputFileError(Exception):
    pass

def RFailure(Exception):
    pass


def checkFileExists(fname):
    if not os.path.exists(fname):
        raise MissingInputFileError('Input file does not exist: %s'%fname)

def run(mapFile, matchFragFile, outputMapFile):
    cmd = ['R', '--quiet', '--slave', 
            '--args',
            mapFile,
            matchFragFile,
            outputMapFile,
            ]
    #sys.stdout.write(' '.join(cmd) + '\n')
    #sys.stdout.flush()
    p = subprocess.Popen(cmd, stdin = subprocess.PIPE)
    script = open(R_SCRIPT).read()
    p.stdin.write(script)
    p.stdin.close()
    ret = p.wait()
    if (ret != 0):
        raise RFailure('R Script returned with error code: %i'%ret)

def makeMatchedChunkFile(matchList, outputFile):
    alns = [mc for mr in matchList for mc in mr.alignment if not (mc.isBoundaryChunk or mc.isContigGap)]
    contigChunks = [al.contig.lengthBp() for al in alns]
    opticalChunks = [al.optical.lengthBp() for al in alns]
    fout = open(outputFile, 'w')
    fout.write('contig\toptical\n')
    for c, o in zip(contigChunks, opticalChunks):
        fout.write('%i\t%i\n'%(c,o))
    fout.close()

def main():

    args = parser.parse_args()

    # Check that input files exist
    die = False
    checkFileExists(args.inputMapFile)
    checkFileExists(args.matchFragFile)
    checkFileExists(R_SCRIPT)

    # Run the R Script
    run(args.inputMapFile, args.matchFragFile, args.outputMapFile)


if __name__ == '__main__':
    main()
