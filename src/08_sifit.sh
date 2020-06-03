#!/bin/bash

set -e -x -u -o pipefail

SIFIT=~/bin/SiFit.jar

java -jar $SIFIT -m 1663 -n 157 -fp 0.002 -fn 0.2 -iter 10000 -df 1 -ipMat sifit/gistic.calls.im3.txt --cellNames sifit/cellID.txt


java -cp $SIFIT SiFit.algorithm.InferAncestralStates -fp 0.002 -fn 0.2 -df 0 -ipMat sifit/gistic.calls.im.txt -tree sifit/gistic.calls.im3_mlTree.newick 
