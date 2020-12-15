#!/bin/bash

dir=`ls -d /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/RNA-Seq/rmdbamfiles/*`
for i in $dir
do
	cd $i
	sample=`basename $i`
	outbam=$sample.merge.bam
	bamlist=`ls *sort.rmd.bam`
	CMD="pkurun-cns 1 12 samtools merge $outbam $bamlist"
	echo $CMD
	eval $CMD
done
