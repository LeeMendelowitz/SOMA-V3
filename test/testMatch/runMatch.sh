#!/bin/bash -x
BIN_DIR=../../bin
MATCH_BIN_DEBUG=$BIN_DIR/match_d
MATCH_BIN_RELEASE=$BIN_DIR/match
# TEST RELEASE
$MATCH_BIN_RELEASE --output test100WithError --pvalue=0.05 --sdMin=4 --sdMax=4 --numPermutationTrials 1000 --numThreads 4 100_WithError.silico 100.opt
$MATCH_BIN_RELEASE --output test100NoError --pvalue=0.05 --sdMin=4 --sdMax=4 --numPermutationTrials 1000 --numThreads 4 100_NoError.silico 100.opt

# TEST DEBUG
#$MATCH_BIN_DEBUG --output test100NoError --pvalue=0.05 --sdMin=4 --sdMax=4 --numPermutationTrials 1000 --numThreads 1 100_NoError.silico 100.opt
#$MATCH_BIN_DEBUG --output contig8 --pvalue=0.05 --sdMin=4 --sdMax=4 --numPermutationTrials 0 --numThreads 1 contig8.silico 100.opt
