#!/bin/bash -l

# all this script does is call Heinz, which is also present in (or linked to) the same directory; Heinz uses the output of writeHeinzNodes() and writeHeinzEdges() in the R package Bionet

CMD="./heinz -n $1.Nodes.txt -e $1.Edges.txt -o $1.Output.txt"
eval $CMD

