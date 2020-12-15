#!/bin/bash

dir=`ls -d /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/DNase-Seq/Data_181022/aligned/merge/*bam`
for i in $dir
do
	cd $i
	sample=`basename $i`
	input=${sample}_R1.fq.gz.clean.fastq.sam.MAPQ30.sort.rmd.bam
	outdir="/lustre1/lch3000_pkuhpc/liuyt/monkey-brain/DNase-Seq/Data_181022/macs2E03/$sample"
	mkdir -p $outdir
	CMD="pkurun-cns 1 1 macs2 callpeak -B -p 1e-3 -g hs --nomodel --extsize 147  -n $sample -t $input --outdir $outdir"
	echo $CMD
	eval $CMD
done
