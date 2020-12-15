#!/bin/bash

ref=/lustre/user/liclab/liuyt/monkey-brain/ref_genome/rheMac8/rheMac8.ucsc.refgene.txt
#i=/lustre/user/liclab/liuyt/monkey-brain/chip-seq/public/tmpdata/bowtie/treat.sorted.bam
outdir=/lustre/user/liclab/liuyt/monkey-brain/chip-seq/ChIP_plot/plot

for i in `ls ../../testdata/tmpdata/bowtie2/*.bam`
do
{
	out=`basename $i`
	mkdir ${outdir}"/"${out}
	nohup python ChIPseqPlot_V1.3.py -f $ref -b $i -o ${outdir}"/"${out} -u 4000 & 
        
}
done
#nohup python ChIPseqPlot_V1.3.py -f $ref -b $i -o $outdir -u 3000 &


