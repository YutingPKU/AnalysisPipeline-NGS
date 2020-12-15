#!/bin/bash
#SBATCH -J rh-cp
#SBATCH -p cn-short
#SBATCH -N 1 
#SBATCH -o rh-cp_%j.out
#SBATCH -e rh-cp_%j.err
#SBATCH --no-requeue
#SBATCH -A lch3000_g1
#SBATCH --qos=lch3000cns
#SBATCH -c 1

samtools bedcov rheMac8.100bp.window.bed  /home/lch3000_pkuhpc/lustre1/liuyt/monkey-brain/RNA-Seq/rmdbamfiles/rheMac8.CP.merged.rmd.bam > rheMac8.CP.100bp.bedgraph 
