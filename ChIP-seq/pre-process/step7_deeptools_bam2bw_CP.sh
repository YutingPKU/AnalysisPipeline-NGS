#!/bin/bash

while read -a line
do 
	bam=$line
	CMD="pkurun-cns 1 20 bamCoverage --bam $bam -o $bam.CPM.bw --binSize 10 --normalizeUsing CPM --extendReads 147 -p 20"
	#CMD="pkurun-cns 1 4 computeMatrix reference-point --referencePoint TSS -p 4 -a 1000 -b 1000 --skipZeros -R /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/ref_genome/ensemble/Mmul_8.0.1/Macaca_mulatta.Mmul_8.0.1.89.chr.gtf -S $bam -o $bam.TSS.gz "
	#CMD="pkurun-cns 1 1 plotProfile -m $bam.TSS.gz -out $bam.TSSplot.png"
	echo $CMD
	eval $CMD
done < rmd.bam.07456BCP.list



#CMD="pkurun-cns 1 4 computeMatrix reference-point --referencePoint center -p 4 -a 1000 -b 1000 --skipZeros -R Rh_12Fpcw_H3K27ac_peaks.narrowPeak -o CP.H3K27ac.gz \
#-S  /home/lch3000_pkuhpc/lustre1/liuyt/monkey-brain/DNase-Seq/aligned/DNase-07456A-CP-100/DNase-07456A-CP-100_R1.fq.gz.sam.sort.rmd.bam.CPM.bw \
#/home/lch3000_pkuhpc/lustre1/liuyt/monkey-brain/DNase-Seq/aligned/DNase-07456A-CP-25/DNase-07456A-CP-25_R1.fq.gz.sam.sort.rmd.bam.CPM.bw \
#/home/lch3000_pkuhpc/lustre1/liuyt/monkey-brain/DNase-Seq/aligned/DNase-07456A-CP-50/DNase-07456A-CP-50_R1.fq.gz.sam.sort.rmd.bam.CPM.bw \
#/home/lch3000_pkuhpc/lustre1/liuyt/monkey-brain/DNase-Seq/aligned/DNase-07456A-CP-75/DNase-07456A-CP-75_R1.fq.gz.sam.sort.rmd.bam.CPM.bw"
#CMD="pkurun-cns 1 1 plotProfile -m ../plot/CP.H3K27ac.gz  -out ../plot/CP.H3K27acplot.png  --perGroup"
#echo $CMD
#eval $CMD

