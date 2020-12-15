#!/bin/bash

while read -a line
do
	CMD="pkurun-cns 1 8 bamCoverage --bam $line -o $line.BedGraph --binSize 100 --outFileFormat bedgraph   -p 8"
#	CMD="pkurun-cns 1 4 samtools index $line"
	echo $CMD
	eval $CMD
done < merged.bam.list
