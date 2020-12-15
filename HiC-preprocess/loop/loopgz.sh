#!/bin/bash
#SBATCH -J gz
#SBATCH -p gpu
#SBATCH -N 1 
#SBATCH -o loop_gz_%j.out
#SBATCH -e loop_gz_%j.err
#SBATCH --no-requeue
#SBATCH -A lch3000_g1
#SBATCH --qos=lch3000gpu
#SBATCH -c 1
time java -Xmx8g -jar /lustre1/lch3000_pkuhpc/liuyt/software/bin/juicer_tools.1.7.6_jcuda.0.8.jar  hiccups -m 1024 -r 2000 -p 6 -i 20 -d 8000  --ignore_sparsity  /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/merged-hic/hicfiles/GZ_1kbcisValidPairs.hic /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/merged-hic/results/GZ-2k
 
