#!/bin/bash

dir=`ls  /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/DNase-Seq/Data_190119/aligned/DNase-07456A-GZ-3/*.rmd.bam`
for i in $dir
do
	#cd $i
	sample=`basename $i`
	#input=${sample}_R1.fq.gz.sam.MAPQ30.sort.rmd.bam
	outdir="/lustre1/lch3000_pkuhpc/liuyt/monkey-brain/DNase-Seq/Data_190119/macs2-replicate/$sample-E05"
	mkdir -p $outdir
	CMD="pkurun-cns 1 1 macs2 callpeak -B -p 1e-5 -g hs --nomodel  --extsize 147   -n $sample -t $i --outdir $outdir"
	echo $CMD
	eval $CMD
done
