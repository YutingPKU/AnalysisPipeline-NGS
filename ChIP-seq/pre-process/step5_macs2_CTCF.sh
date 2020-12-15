#!/bin/bash

input="/lustre1/lch3000_pkuhpc/liuyt/monkey-brain/CTCF/aligned/07456A-CP-MAPQ30.sort.rmd.merge.bam"
#input="../aligned/07456A-GZ-CTCF-1/07456A-GZ-CTCF-1_R1.fq.gz.clean.fastq.sam.MAPQ30.sort.rmd.bam"
ctrl="../aligned/Input-07456A-CP-1/Input-07456A-CP-1_R1.fq.gz.clean.fastq.sam.MAPQ30.sort.rmd.bam"
out="../macs2/07456A-CP-merge"
mkdir -p $out
CMD="pkurun-cns 1 1 macs2 callpeak -B -g hs -p 1e-5  --nomodel --extsize 200 -n 07456A-CP-merge-CTCF -t $input -c $ctrl --outdir $out"
echo $CMD
eval $CMD
