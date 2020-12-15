#!/bin/bash

#if [ $# -lt 2 ];then
#	    echo "Need 2 parameters! <bedgraph> <chrom info>"
#		    exit
#fi

while read -a line
do
	F=$line
	G=/lustre1/lch3000_pkuhpc/liuyt/monkey-brain/ref_genome/rheMac8/rheMac8.brief.chrom.sizes
	CMD="pkurun-cns 1 1 bash bdg2bw.sh $F $G"
	echo $CMD 
	eval $CMD
done < bdg.list

