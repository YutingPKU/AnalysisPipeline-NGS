#!/bin/bash

while read -a line 
do
	echo $line
	CMD=`samtools idxstats $line | awk '{s+=$3}END{print s}'`
	echo $CMD
done < rmd.bam.list 
