#!/bin/bash
/cbcb/project-scratch/lmendelo/soma/test/testMatch
BIN_DIR=../../bin
MATCH_BIN=$BIN_DIR/match
$MATCH_BIN --output test100WithError --local --pvalue=0.05 --sdMin=4 --sdMax=4 --numPermutationTrials 1000 --numThreads 4 100_WithError.silico 100.opt
$MATCH_BIN --output test100NoError --local --pvalue=0.05 --sdMin=4 --sdMax=4 --numPermutationTrials 1000 --numThreads 4 100_NoError.silico 100.opt
