#!/bin/bash

ref=/lustre/user/liclab/liuyt/monkey-brain/ref_genome/rheMac3/rheMac3.ucsc.refgene.txt
i=/lustre/user/liclab/liuyt/monkey-brain/chip-seq/public3/tmpdata/bowtie/treat.sort.rmd.bam
outdir=/lustre/user/liclab/liuyt/monkey-brain/chip-seq/ChIP_plot/plot
mkdir $outdir"/public3"
nohup python ChIPseqPlot_V1.3.py -f $ref -b $i -o ${outdir}"/public3/" -u 4000 & 
#nohup python ChIPseqPlot_V1.3.py -f $ref -b $i -o $outdir -u 3000 &


