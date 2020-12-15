#!/bin/bash

while read -a line
do
	CMD="pkurun-cns 1 4 samtools index -@ 4 $line"
	echo $CMD
	eval $CMD
done < merge.bam.list
