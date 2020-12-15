#!/bin/bash


#paste R1.clean.fastq.list R2.clean.fastq.list | while IFS="$(printf '\t')" read -r f1 f2
while read -a line
do
	bam=$line.bam 
	ts=$line.sort.bam
	trm=$line.sort.rmd.bam
	trmm=$line.mertrics
#	CMD="pkurun-cns 1 20 samtools view -@ 20 -bS -q 30 -o $line.MAPQ30.bam $line "
#	CMD="pkurun-cns 1 20 samtools sort -@ 20 $line.MAPQ30.bam -o $line.MAPQ30.sort.bam "
	CMD="pkurun-cns 1 4 samtools index $line.MAPQ30.sort.rmd.bam"
#	CMD="pkurun-cns 1 4 samtools index $line"
	#CMD="pkurun-cns 1 4 samtools view -@ 4 -b -o $tsnomt $ts 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"
	#CMD="pkurun-cns 1 1 java -jar /home/lch3000_pkuhpc/lustre1/liuyt/software/picard-tools-1.118/MarkDuplicates.jar INPUT=$line.MAPQ30.sort.bam OUTPUT=$line.MAPQ30.sort.rmd.bam VALIDATION_STRINGENCY=SILENT  REMOVE_DUPLICATES=true METRICS_FILE=$trmm"
	echo $CMD
	eval $CMD
	
	
done < sam.list
