#!/bin/bash

#paste R1.SRA.list.txt R2.SRA.list.txt | while IFS="$(printf '\t')" read -r f1 f2
while read -a line
do
	ts=$line.sort.bam
	trm=$line.sort.rmd.bam
	trmm=$line.mertrics
	#CMD="pkurun-cns 1 12  samtools sort -@ 12 -f $line $ts"
	#echo $CMD
	#eval $CMD
	CMD="pkurun-cns 1 12 java -jar /home/lch3000_pkuhpc/lustre1/liuyt/software/picard-tools-1.118/MarkDuplicates.jar INPUT=$ts OUTPUT=$trm VALIDATION_STRINGENCY=SILENT  REMOVE_DUPLICATES=true METRICS_FILE=$trmm"
#	CMD="pkurun-cns 1 12 samtools index $trm"
	echo $CMD
	eval $CMD
	
	
done <  align.bam.list.2 
