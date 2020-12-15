#!/bin/bash

if [ $# -lt 2 ];then
	    echo "Need 2 parameters! <bedgraph> <chrom info>"
		    exit
fi
		
F=$1
G=$2
		
bedtools slop -i ${F} -g ${G} -b 0 | bedClip stdin ${G} ${F}.clip
		
LC_COLLATE=C sort -k1,1 -k2,2n ${F}.clip > ${F}.sort.clip
		
bedGraphToBigWig ${F}.sort.clip ${G} ${F/bdg/bw}
		
rm -f ${F}.clip ${F}.sort.clip""

