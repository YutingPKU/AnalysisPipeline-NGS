#!/bin/bash

paste R1.clean.fastq.list R2.clean.fastq.list | while IFS="$(printf '\t')" read -r f1 f2
do
	path=`dirname $f1`
	dir=`basename ${path}`
	R1=`basename $f1`
	R2=`basename $f2`
	outdir=/lustre1/lch3000_pkuhpc/liuyt/monkey-brain/CTCF/aligned/$dir
	mkdir -p ${outdir}
	CMD="pkurun-cnlong 1 20 bowtie2 -p 20  --very-sensitive -x /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/ref_genome/rheMac8/Bowtie2Index/rheMac8  -1 $f1 -2 $f2 -S $outdir/$R1.sam"
	echo $CMD
	eval $CMD
	
	
done 
