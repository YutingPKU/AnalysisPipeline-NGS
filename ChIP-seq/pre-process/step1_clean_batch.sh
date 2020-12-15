#!/bin/bash


paste R1.raw.fq.list R2.raw.fq.list | while IFS="$(printf '\t')" read -r f1 f2
#while read -a line
do
	R1=`basename $f1`
	R2=`basename $f2`
	dir=`dirname $f1`
	sample=`basename $dir`
	outdir="/lustre1/lch3000_pkuhpc/liuyt/monkey-brain/CTCF/Cleandata/$sample"
	mkdir -p $outdir
	CMD="pkurun-cnlong 1 4 cutadapt  -q 20,15 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAA -o ${outdir}/${R1}.clean.fastq -p ${outdir}/${R2}.clean.fastq $f1 $f2"
#	R1=`basename $line`
#	out=../cleandata/$R1.clean.fastq
#	CMD="pkurun-cns 1 12 cutadapt -j 12 -m 30 -q 20,15 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o $out $line"
	echo $CMD
	eval $CMD
	
	
done 
