#!/bin/bash
#SBATCH -J loopcp
#SBATCH -p gpu
#SBATCH -N 1 
#SBATCH -o diffloop_%j.out
#SBATCH -e diffloop_%j.err
#SBATCH --no-requeue
#SBATCH -A lch3000_g1
#SBATCH --qos=lch3000gpu
#SBATCH -c 1
time java -Xmx8g -jar /lustre1/lch3000_pkuhpc/liuyt/software/bin/juicer_tools.1.7.6_jcuda.0.8.jar  hiccupsdiff --ignore_sparsity -m 1024 -r 10000 -k KR -f .1 -p 2 -i 5 -t 0.02,1.5,1.75,2 -d 20000   /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/merged-hic/hicfiles/CP_cis_allValidPairs.hic /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/merged-hic/hicfiles/GZ_cis_allValidPairs.hic /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/merged-hic/results/CP-10k/merged_loops.bedpe /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/merged-hic/results/GZ-10k/merged_loops.bedpe /lustre1/lch3000_pkuhpc/liuyt/monkey-brain/merged-hic/results/diffloop  
